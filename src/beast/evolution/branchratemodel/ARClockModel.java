package beast.evolution.branchratemodel;



import java.util.Arrays;
import java.util.List;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import com.sun.xml.internal.ws.message.MimeAttachmentSet;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.branchratemodel.BranchRateModel.Base;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.LogNormalDistributionModel;
import beast.math.distributions.LogNormalDistributionModel.LogNormalImpl;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;
import util.DistTest;

/**
 * @author Igor Siveroni
 */

@Description("Implements ARC: th Additive Uncorrelated Ralexed Clock Model")
@Citation(value =
        "Didelot X and Volz EM (2019) Additive uncorrelated  relaxed clock models\n" +
                "  for the dating of epidemiological models.", DOI = "TBC",
        year = 2019, firstAuthorSurname = "Didelot")

public class ARClockModel extends Base {

    final public Input<IntegerParameter> categoryInput = new Input<>("rateCategories", "the rate categories associated with nodes in the tree for sampling of individual rates among branches.");
    final public Input<RealParameter> rateProbsInput = new Input<>("rateProbs", "the accum probabilities of the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.XOR, categoryInput);   
    final public Input<Integer> numberOfDiscreteRates = new Input<>("numberOfDiscreteRates", "the number of discrete rate categories to approximate the rate distribution by. A value <= 0 will cause the number of categories to be set equal to the number of branches in the tree. (default = -1)", -1);
    final public Input<Tree> treeInput = new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);  
    // Gamma distribution parameters
    final public Input<RealParameter> ratesMeanInput = new Input<>("rateMean", "mean of rates Gamma probability distribution", Input.Validate.REQUIRED);
    final public Input<RealParameter> ratesOmegaInput = new Input<>("rateOmega", "omega factor of rates Gamma distribution", Input.Validate.REQUIRED);
    final public Input<Integer> numSitesInput = new Input<>("numSites", "Number of sites (sequence length)");
    
    boolean useCategories;
    // either categories or rateProbsParameter is used
    RealParameter rateProbs; //when mode=rates
    IntegerParameter categories; //when mode=categories
      
    //ContinuousDistribution gammaDist;
    RealParameter ratesMean;
    RealParameter ratesOmega;

    RealParameter meanRate;
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
        
				
		//TaxonSet taxa = tree.getTaxonset();	
		//List<Taxon> lt = taxa.taxonsetInput.get();
		//if (lt==null) System.out.println("lsyt if taxon null");
		//else System.out.println("list lengh = "+lt.size());
		
		//Alignment alignment = taxa.alignmentInput.get();
		//if (alignment==null) {
		//	System.out.println("Alignment is NULL");
		//} else {
		//int count = alignment.getSiteCount();
		//System.out.println("------------------ Counts sites = "+count );
		//}
		
		if (numSitesInput.get()==null) 
			numSites = 1;
		else
			numSites = numSitesInput.get();
		System.out.println("ARC model. numSites ="+numSites);
		
	    ratesMean = ratesMeanInput.get();
	    ratesOmega = ratesOmegaInput.get();
        
        if (categories==null) useCategories=false;
        else useCategories = true;
        
        if (useCategories) {
            numCategories = numberOfDiscreteRates.get();
            if (numCategories <= 0) numCategories = branchCount;
            Log.info.println("  ARClockModel: using " + numCategories + " rate " +
                    "categories to approximate rate distribution across branches.");
        } else {
            if (numberOfDiscreteRates.get() != -1) {  // if it's not the default i.e. it has been set as input
                throw new IllegalArgumentException("Can't specify both numberOfDiscreteRates and  rateProbs inputs.");
            }
            Log.info.println("  ARClockModel: using rate probs as parameter for sampling across branches.");          
        }
        
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
        			//System.out.println("p_i="+initialProbs[i]);
        		}
        		RealParameter other = new RealParameter(initialProbs);
        		rateProbs.assignFromWithoutID(other);
        	}
        	rateProbs.setLower(0.0);
        	rateProbs.setUpper(1.0);
        }
        
        
        meanRate = meanRateInput.get();  // from superclass
        if (meanRate == null) {
            meanRate = new RealParameter("1.0");
        }
        
        if (ratesMean == null) {
        	Log.warning.println("rateMean is null");
        } else {
        	Log.warning.println("rateMean :"+ratesMean.getValue());
        }
        
        //if (Math.abs(ratesMean.getValue()-1.0) > 1e-6 ) { // igor: should this be ONE in ARC's case?
        //	Log.warning.println("WARNING: mean of distribution for additive relaxed clock model is not 1.0.");
        //}
       
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
        //DistTest.timeDistributions();
      
	}
    
 
	

	//get the rate for node: R = r*scale*meanRate
	@Override
	public double getRateForBranch(Node node) {
		//startTime();
		
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
		//System.out.println("rate = "+r+" nr="+node.getNr()+"  newNr="+getNr(node));
		//endTime(false);
        return r;	
     
	}
	
	
	private void calculateUnscaledBranchRates() {
		//System.out.println("recalculate branch rates");
		//debugState();
		//debugTree();
		if (this.useCategories) {
			calculateRatesForCategories();
		} else {
			calculateRatesForProbs();
		}
	}
    
	void debugState() {
		//System.out.println(tree.getRoot().getNr()+ " =? " + (this.nodeCount-1));
		if (tree.getRoot().getNr() != (nodeCount-1)) {
			System.out.println("Root Swap!!");
		}
		if (this.useCategories) {
			Integer[] vs = categories.getValues();
			//System.out.println("array size = "+vs.length+ " num branches="+branchCount);
			System.out.print("Cs = ");
			for(int i=0; i < vs.length; i++) {
				System.out.print(vs[i]+" ");
			}
			System.out.println(" ");
		} else {
			Double[] vs = rateProbs.getValues();
			System.out.print("Ps = ");
			for(int i=0; i < vs.length; i++) {
				System.out.print(vs[i]+" ");
			}
			System.out.println(" ");
		}
	}
	
	void debugTree() {
		for(int i=0; i < tree.getNodeCount();i++) {
    		final Node node = tree.getNode(i);
    		System.out.print(node.getNr()+"/"+String.format("%.3f", node.getLength())+" ");
		}
		System.out.println("");
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
    				//System.out.println("nr = "+nr+" l ="+l+" "+branchLengths[node.getNr()]);   				
    			}
    			
    			  			
    			final int category = categories.getValue(nr);
    			final double p = (category  + 0.5) / numCategories;
    			
    			final double theta = ratesOmega.getValue() /  l  ;
    			final double k = ratesMean.getValue() / theta;
    			
    			gammaDist = new GammaDistributionImpl(k*numSites, theta); 
    			
    			branchLengths[node.getNr()] = l;
    			
    			//System.out.println("k = "+k+ " theta="+theta+"  p="+p);
    			try {
    				rate = gammaDist.inverseCumulativeProbability( p ) / numSites;
    			} catch (MathException e) {
    				inverseCumError(e,k,theta,p);
    			}
    			rates[nr] = rate;  			
    		} else {  // node is root
    			branchLengths[node.getNr()] = node.getLength(); 
    		}
    		// System.out.println(" rate = "+rate);
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
    				//System.out.println("nr = "+nr+" l ="+l+" "+branchLengths[node.getNr()]);   				
    			}    			   			  			
    			
    			final double p = rateProbs.getValue(nr);
    			
    			final double theta = ratesOmega.getValue() /  l  ;
    			final double k = ratesMean.getValue() / theta;
    			
    			gammaDist = new GammaDistributionImpl(k*numSites, theta); 
    			
    			branchLengths[node.getNr()] = l;
    			
    			//System.out.println("k = "+k+ " theta="+theta+"  p="+p);
    			try {
    				rate = gammaDist.inverseCumulativeProbability( p ) / numSites;
    			} catch (MathException e) {
    				inverseCumError(e,k,theta,p);
    			}
    			rates[nr] = rate;  			
    		} else {  // node is root
    			branchLengths[node.getNr()] = node.getLength(); 
    		}
    		// System.out.println(" rate = "+rate);
    	}
    	samplingOnly=false;
    	treeOnly = false;
    	
    	
    }
    
    
   private  void inverseCumError(MathException e, double k, double theta, double p) {
	   System.out.println("Problems computing Gamma Inverse Cumm. Prob. :");
	   System.out.println("mu = "+ratesMean.getValue());
	   System.out.println("omega = "+ratesOmega.getValue());
	   System.out.println("kappa ="+k+"  theta = "+theta);
	   System.out.println("p = "+p);
	   throw new RuntimeException("Failed to compute inverse cumulative probability");
   }
    
	
    /*  computes the idx used to access the arrays a[branchCount] 
     * needed in case of node swaps involving the root */
    private int getNr(Node node) {
        int nodeNr = node.getNr();
        if (nodeNr > tree.getRoot().getNr()) {
        	System.out.println("-- root swap");
            nodeNr--;
        }
        return nodeNr;
        
    }

    // indexing used by UCRelaxedClock
    private int getNr2(Node node) {
    	int nodeNr = node.getNr();
        if (nodeNr == branchCount) { // this should be root if no swap has happened
            // root node has nr less than #categories, so use that nr
            nodeNr = node.getTree().getRoot().getNr();
        }
        return nodeNr;
    }
    
    // Come back when parameters have settled - so far, recompute always
    @Override
    protected boolean requiresRecalculation() {
        recompute = false;
        samplingOnly=true;
        
        treeOnly= true;
        //System.out.println("calling requires recalculation");
        if (treeInput.get().somethingIsDirty()) {
        	//System.out.println("tree changes");
        	samplingOnly=false;
        	recompute = true;
        	//return true;
        }
        
        // changed from isDirtyCalculation based on ratesOmegaInput feedback
        if (ratesMeanInput.get().somethingIsDirty()) {
        	//System.out.println("Gamma mu changes");
        	samplingOnly=false;
        	treeOnly=false;
        	recompute = true;
        	//return true;
        }
     
        //System.out.println("omega="+ratesOmegaInput.get().getArrayValue());
        
        //if (ratesOmegaInput.get().isDirty(0)) {
        //	System.out.println("------ omega isDIRTY!!!");
        //}
        
        // isDirtyCalculation did not report changes ratesOmegaInput
        // somethingIsDirty seems to over-report them
        if (ratesOmegaInput.get().somethingIsDirty()) {
        	//System.out.println("Gamma omega changes");
        	samplingOnly=false;
        	treeOnly=false;
        	recompute = true;
        	//return true;
        }
     
        if (useCategories) {
        	if ( categoryInput.get().somethingIsDirty() ) {
        		//System.out.println("category changes");
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
        if (meanRate.somethingIsDirty()) {
        	return true;
        }           
        //System.out.println("recompute:"+recompute);
        return recompute;
    }
    
    @Override
    public void store() {
    	System.arraycopy(rates, 0, storedRates, 0, rates.length);
    	System.arraycopy(branchLengths, 0, storedBranchLengths, 0, branchLengths.length);
        super.store();
        // System.out.println("num calls="+numCalls+"  duration="+duration);
    }

    @Override
    public void restore() {
    	//System.out.println("--- restore ARC");
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
