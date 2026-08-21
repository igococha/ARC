/*
 * ARC: Additive uncorrelated Relaxed Clock models for BEAST 2.
 *
 * Copyright (C) 2019-2026 Xavier Didelot, Igor Siveroni and Erik Volz.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
package arc.evolution.branchratemodel;




import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel.Base;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;

/**
 * @author Igor Siveroni
 */

@Description("Additive uncorrelated relaxed clock (ARC). The rate of a branch of duration t "
        + "is Gamma distributed with mean rateMean and variance rateMean*rateOmega/t, so the "
        + "variance of the number of substitutions along a branch grows linearly with its "
        + "duration. That additivity makes the clock consistent under subdivision of a branch.")
@Citation(value =
        "Didelot X, Siveroni I, Volz EM (2021) Additive Uncorrelated Relaxed Clock Models\n" +
        "  for the Dating of Genomic Epidemiology Phylogenies.\n" +
        "  Molecular Biology and Evolution 38(1):307-317.",
        DOI = "10.1093/molbev/msaa193",
        year = 2021, firstAuthorSurname = "Didelot")

public class ARClockModel extends Base {

    final public Input<IntegerParameter> categoryInput = new Input<>("rateCategories", "the rate categories associated with nodes in the tree for sampling of individual rates among branches.");
    
    final public Input<RealParameter> rateProbsInput = new Input<>("rateProbs", "the accum probabilities of the rates associated with nodes in the tree for sampling of individual rates among branches.");   
    final public Input<Integer> numberOfDiscreteRates = new Input<>("numberOfDiscreteRates", "the number of discrete rate categories to approximate the rate distribution by. A value <= 0 will cause the number of categories to be set equal to the number of branches in the tree. (default = -1)", -1);
    final public Input<Tree> treeInput = new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);  
    // Gamma distribution parameters
    final public Input<RealParameter> ratesMeanInput = new Input<>("rateMean", "mean of rates Gamma probability distribution", Input.Validate.REQUIRED);
    final public Input<RealParameter> ratesOmegaInput = new Input<>("rateOmega", "omega factor of rates Gamma distribution", Input.Validate.REQUIRED);
    
    final public Input<Integer> numSitesInput = new Input<>("numSites", "Number of sites (sequence length)");
    final public Input<Alignment> alignmentInput = new Input<>("alignment", "Sequence alignment"); 
    
    boolean useCategories;
    // either categories or rateProbsParameter is used
    RealParameter rateProbs; //when mode=rates
    IntegerParameter categories; //when mode=categories
      
    //ContinuousDistribution gammaDist;
    RealParameter ratesMean;
    RealParameter ratesOmega;

    Function meanRate;
    Tree tree;
    private int branchCount;//the number of branches of the tree
    private int nodeCount;
    
    int numCategories = 100;
    int numSites;
    
    
    private boolean recompute = true;
    private boolean samplingOnly = false;
    
    private boolean treeOnly = false;
    private double[] branchLengths; // accessed with nr
    private double[] storedBranchLengths; 
    
    double[] rates;
    double[] storedRates;
   

    @Override
	public void initAndValidate() {
		tree = treeInput.get();
		nodeCount = tree.getNodeCount();
		branchCount = nodeCount - 1;						
		if (categoryInput.get() != null) { 
			categories = categoryInput.get();
			useCategories=true;
		} else if (rateProbsInput.get() != null) {
			rateProbs = rateProbsInput.get();
			useCategories = false;
		} else {
			throw new IllegalArgumentException("Input error: Must use rateCategories or rateProbs");
		}
		
		if (alignmentInput.get() != null) {
			numSites = alignmentInput.get().getSiteCount();
			if (numSitesInput.get() != null) {
				Log.warning.println("ARClockModel: both alignment and numSites given; using the alignment (" + numSites + " sites).");
			}
		} else if (numSitesInput.get() != null) {
			numSites = numSitesInput.get();
		} else {
			// Not a benign default. numSites enters the Gamma shape as
			// rateMean*l*numSites/rateOmega, so with numSites=1 the shape is
			// typically far below 1 and the rate distribution degenerates into a
			// spike at zero with a very long tail, which can drive the likelihood
			// to NaN.
			numSites = 1;
			Log.warning.println("ARClockModel: WARNING neither alignment nor numSites was given, so numSites "
					+ "defaults to 1. This is very probably not what you want: it changes the model, and can "
					+ "make the likelihood degenerate. Supply the alignment, or set numSites to the sequence length.");
		}
		
		
	    ratesMean = ratesMeanInput.get();
	    ratesOmega = ratesOmegaInput.get();
        
        if (categories==null) useCategories=false;
        else useCategories = true;
        
        if (useCategories) {
            numCategories = numberOfDiscreteRates.get();
            if (numCategories <= 0) numCategories = branchCount;
            Log.info.println("ARClockModel: using " + numCategories + " rate " +
                    "categories to approximate rate distribution across branches.");
        } else {
            if (numberOfDiscreteRates.get() != -1) {  // if it's not the default i.e. it has been set as input
                throw new IllegalArgumentException("Can't specify both numberOfDiscreteRates and  rateProbs inputs.");
            }
            Log.info.println("ARClockModel: using rate probs as parameter for sampling across branches.");          
        }
        Log.info.println("ARClockModel: " + numSites + " sites.");
        
        //initialize rates in both types of inputs
        if (useCategories) {       	
        	categories.setDimension(branchCount);
        	Integer[] initialCategories = new Integer[branchCount];
        	for (int i = 0; i < branchCount; i++) {
        		initialCategories[i] = Randomizer.nextInt(numCategories);
        	}
        	// set initial values of rate categories
        	IntegerParameter other = new IntegerParameter(initialCategories);
        	categories.assignFromWithoutID(other);
        	categories.setLower(0);
        	categories.setUpper(numCategories - 1);
        } else {        
        	if (rateProbs.getDimension() != branchCount) {
        		rateProbs.setDimension(branchCount);
        		//randomly draw rates from the lognormal distribution
        		Double [] initialProbs = new Double[branchCount];
        		for (int i = 0; i < branchCount; i++) {
        			initialProbs[i] = Randomizer.nextDouble();  // [0,1)
        		}
        		RealParameter other = new RealParameter(initialProbs);
        		rateProbs.assignFromWithoutID(other);
        	}
        	rateProbs.setLower(0.0);
        	rateProbs.setUpper(1.0);
        }
        
        
        // clock.rate, from the superclass. Only a value supplied in the XML counts:
        // the fallback below is a bare parameter that never enters the State, and it
        // reports estimate=true purely because that is the Input default.
        final Function suppliedMeanRate = meanRateInput.get();
        meanRate = (suppliedMeanRate == null) ? new RealParameter("1.0") : suppliedMeanRate;

        // clock.rate and rateMean are both mean clock rates: a branch rate has mean
        // clock.rate * rateMean. Estimating both is exactly non-identifiable -- the
        // chain wanders along the ridge where their product is constant -- so refuse
        // it rather than silently produce an unusable posterior. Fixing either one
        // is fine; clock.rate defaults to 1.0, which is the usual choice.
        if (isEstimated(suppliedMeanRate) && isEstimated(ratesMean)) {
            throw new IllegalArgumentException(
                    "ARClockModel: clock.rate and rateMean are both mean clock rates, "
                  + "so estimating both makes the model non-identifiable. "
                  + "Fix one of them (clock.rate defaults to 1.0).");
        }
       
        // new
        // unscaledBranchRates = new double[tree.getNodeCount()];
        rates = new double[branchCount];
        storedRates = new double[branchCount];
        branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];
        for(int i=0; i < nodeCount; i++) {
        	branchLengths[i]=-1;
        }
		//testGamma(0.020132723470917546, 1/ 49.67037874654845 );
      
	}
    
 
	

	/** Branch rate R = r * clock.rate, where r is the ARC-distributed rate. */
	@Override
	public double getRateForBranch(Node node) {
		
		if (node.isRoot()) {
			//endTime(false);
            return 1; // root has no rate
        }
		synchronized (this) {
    		if (recompute) {
    			calculateUnscaledBranchRates();
    			//if (normalize)
    			// recalculateScaleFactor();
    			recompute = false;
			}
        }
		double r = rates[getNr(node)];
		//endTime(false);
        return r * meanRate.getArrayValue();

	}

    /** True if f was supplied in the XML as a state node the MCMC is sampling. */
    private static boolean isEstimated(Function f) {
        return (f instanceof StateNode) && ((StateNode) f).isEstimatedInput.get();
    }
	
	
	private void calculateUnscaledBranchRates() {
		if (this.useCategories) {
			calculateRatesForCategories();
		} else {
			calculateRatesForProbs();
		}
	}
    
	
    
	
	// includes all optimisations in a single function
    private void calculateRatesForCategories() {   	
    	
    	double rate=1.0;
    	ContinuousDistribution gammaDist;
    	//tree.getNodesAsArray();
    	for(int i=0; i < tree.getNodeCount();i++) {
    		final Node node = tree.getNode(i);
    		if (! node.isRoot()) {   			
    			final int nr = getNr(node);
    			if (samplingOnly && !categories.isDirty(nr)) {
    				continue;
    			}
    			final double l = node.getLength();  
    			   			
    			if (treeOnly) {
    				// it seems that checking for node.isDirty() ==Tree.IS_CLEAN is not enough
    				final double diff = l - branchLengths[node.getNr()];
    				final boolean isDiff = Math.abs(diff) > 0.00000001;
    				if (!isDiff) {
    					continue;
    				}
    			}
    			
    			  			
    			final int category = categories.getValue(nr);
    			final double p = (category  + 0.5) / numCategories;
    			
    			final double theta = ratesOmega.getValue() /  l  ;
    			final double k = ratesMean.getValue() / theta;
    			
    			gammaDist = new GammaDistributionImpl(k*numSites, theta); 
    			
    			branchLengths[node.getNr()] = l;
    			
    			try {
    				rate = gammaDist.inverseCumulativeProbability( p ) / numSites;
    			} catch (MathException e) {
    				inverseCumError(e,k,theta,p);
    			}
    			rates[nr] = rate;  			
    		} else {  // node is root
    			branchLengths[node.getNr()] = node.getLength(); 
    		}
    	}
    	samplingOnly=false;
    	treeOnly = false;
    }
    
 
    private void calculateRatesForProbs() {
    	double rate=1.0;
    	ContinuousDistribution gammaDist;
    	//tree.getNodesAsArray();
    	for(int i=0; i < tree.getNodeCount();i++) {
    		final Node node = tree.getNode(i);
    		if (! node.isRoot()) {   			
    			final int nr = getNr(node);
    			if (samplingOnly && !rateProbs.isDirty(nr)) {
    				continue;
    			}
    			final double l = node.getLength();  
    			   			
    			if (treeOnly) {
    				// it seems that checking for node.isDirty() ==Tree.IS_CLEAN is not enough
    				final double diff = l - branchLengths[node.getNr()];
    				final boolean isDiff = Math.abs(diff) > 0.00000001;
    				if (!isDiff) {
    					continue;
    				}
    			}    			   			  			
    			
    			final double p = rateProbs.getValue(nr);
    			
    			final double theta = ratesOmega.getValue() /  l  ;
    			final double k = ratesMean.getValue() / theta;
    			
    			gammaDist = new GammaDistributionImpl(k*numSites, theta); 
    			
    			branchLengths[node.getNr()] = l;
    			
    			try {
    				rate = gammaDist.inverseCumulativeProbability( p ) / numSites;
    			} catch (MathException e) {
    				inverseCumError(e,k,theta,p);
    			}
    			rates[nr] = rate;  			
    		} else {  // node is root
    			branchLengths[node.getNr()] = node.getLength(); 
    		}
    	}
    	samplingOnly=false;
    	treeOnly = false;
    	
    	
    }
    
    
   private void inverseCumError(MathException e, double k, double theta, double p) {
	   throw new RuntimeException("ARClockModel: failed to compute the Gamma inverse cumulative probability"
			   + " (rateMean=" + ratesMean.getValue()
			   + ", rateOmega=" + ratesOmega.getValue()
			   + ", shape=" + k + ", scale=" + theta + ", p=" + p + ")", e);
   }
    
	
    /*  computes the idx used to access the arrays a[branchCount] 
     * needed in case of node swaps involving the root */
    private int getNr(Node node) {
        int nodeNr = node.getNr();
        if (nodeNr > tree.getRoot().getNr()) {
            // the root is not the highest numbered node, so shift down to keep
            // rates[] densely indexed over branches
            nodeNr--;
        }
        return nodeNr;
        
    }
    
    // Come back when parameters have settled - so far, recompute always
    @Override
    protected boolean requiresRecalculation() {
        recompute = false;
        samplingOnly=true;
        
        treeOnly= true;
        if (treeInput.get().somethingIsDirty()) {
        	samplingOnly=false;
        	recompute = true;
        	//return true;
        }
        
        // changed from isDirtyCalculation based on ratesOmegaInput feedback
        if (ratesMeanInput.get().somethingIsDirty()) {
        	samplingOnly=false;
        	treeOnly=false;
        	recompute = true;
        	//return true;
        }
     
        
        //if (ratesOmegaInput.get().isDirty(0)) {
        //}
        
        // isDirtyCalculation did not report changes ratesOmegaInput
        // somethingIsDirty seems to over-report them
        if (ratesOmegaInput.get().somethingIsDirty()) {
        	samplingOnly=false;
        	treeOnly=false;
        	recompute = true;
        	//return true;
        }
     
        if (useCategories) {
        	if ( categoryInput.get().somethingIsDirty() ) {
            	treeOnly=false;
            	recompute = true;  // come back to this
            	//return true;
        	}
        } else { // use probs
        	if (rateProbsInput.get().somethingIsDirty()) {
        		treeOnly=false;
            	recompute = true;
        	}
        }
        
        
        // still don't know the use of this rate
        if (InputUtil.isDirty(meanRateInput)) {
        	return true;
        }           
        return recompute;
    }
    
    @Override
    public void store() {
    	System.arraycopy(rates, 0, storedRates, 0, rates.length);
    	System.arraycopy(branchLengths, 0, storedBranchLengths, 0, branchLengths.length);
        super.store();
    }

    @Override
    public void restore() {
    	double[] tmp = rates;
        rates = storedRates;
        storedRates = tmp;
        // restore branchLengths
        tmp = branchLengths;
        branchLengths = storedBranchLengths;
        storedBranchLengths = tmp;
        super.restore();
    }
      
    
  	
	
	

}
