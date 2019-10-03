package beast.evolution.substitutionmodel;


import java.util.Arrays;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.branchratemodel.BranchRateModel.Base;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.LogNormalDistributionModel;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

/**
 * @author Igor Siveroni
 */

@Description("Implements ARC: th Additive Uncorrelated Ralezed Clock Model")
@Citation(value =
        "Didelot X and Volz EM (2019) Additive uncorrelated  relaxed clock models\n" +
                "  for the dating of epidemiological models.", DOI = "TBC",
        year = 2019, firstAuthorSurname = "Didelot")

public class ARClockModel extends Base {
    final public Input<IntegerParameter> categoryInput = new Input<>("rateCategories", "the rate categories associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    final public Input<Integer> numberOfDiscreteRates = new Input<>("numberOfDiscreteRates", "the number of discrete rate categories to approximate the rate distribution by. A value <= 0 will cause the number of categories to be set equal to the number of branches in the tree. (default = -1)", -1);
    final public Input<RealParameter> rateProbsInput = new Input<>("rateProbs", "the accum probabilities of the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.XOR, categoryInput);
    final public Input<Tree> treeInput = new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
    final public Input<Boolean> normalizeInput = new Input<>("normalize", "Whether to normalize the average rate (default false).", false);
   
    // Gamma distribution parameters
    final public Input<RealParameter> ratesMeanInput = new Input<>("rateMean", "mean of rates Gamma probability distribution", Input.Validate.REQUIRED);
    final public Input<RealParameter> ratesOmegaInput = new Input<>("rateOmega", "omega factor of rates Gamma distribution", Input.Validate.REQUIRED);

	
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
    private boolean normalize = false;
    
    int numCategories = 100;
    
    private boolean recomputeScaleFactor = true;
    
    private boolean recompute = true;
    
    // state that we might want to store/restore
    
    double[] unscaledBranchRates;
    double[] storedUnscaledBranchRates;
    private double scaleFactor = 1.0; 
    private double storedScaleFactor = 1.0;
    

    @Override
	public void initAndValidate() {
		tree = treeInput.get();
		nodeCount = tree.getNodeCount();
		branchCount = nodeCount - 1;
		categories = categoryInput.get();
		rateProbs = rateProbsInput.get();
		
	    ratesMean = ratesMeanInput.get();
	    ratesOmega = ratesOmegaInput.get();
        
        if (categories==null) useCategories=true;
        else useCategories = true;
        
        if (useCategories) {
            numCategories = numberOfDiscreteRates.get();
            if (numCategories <= 0) numCategories = branchCount;
            Log.info.println("  ARClockModel: using " + numCategories + " rate " +
                    "categories to approximate rate distribution across branches.");
        } else {
            if (numberOfDiscreteRates.get() != -1) {  // if it's not the default i.e. it has been set as input
                throw new IllegalArgumentException("Can't specify both numberOfDiscreteRates and  rates inputs.");
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
        			initialProbs[i] =Randomizer.nextDouble();  // [0,1)
        		}
        		RealParameter other = new RealParameter(initialProbs);
        		rateProbs.assignFromWithoutID(other);
        	}
        	rateProbs.setLower(0.0);
        	rateProbs.setUpper(1.0);
        }
        
       
        normalize = normalizeInput.get();
        
        meanRate = meanRateInput.get();  // from superclass
        if (meanRate == null) {
            meanRate = new RealParameter("1.0");
        }
        
        if (ratesMean == null) {
        	Log.warning.println("rateMean is null");
        } else {
        	Log.warning.println("rateMean :"+ratesMean.getValue());
        }
        
        if (Math.abs(ratesMean.getValue()-1.0) > 1e-6 ) { // igor: should this be ONE in ARC's case?
        	Log.warning.println("WARNING: mean of distribution for additive relaxed clock model is not 1.0.");
        }
       
        // new
        unscaledBranchRates = new double[tree.getNodeCount()];
		
	}
	
    private void initialize() {
        // code to be executed whenever we start a fresh MCMC run (and no restore has been called)
    }
	

	//get the rate for node: R = r*scale*meanRate
	@Override
	public double getRateForBranch(Node node) {
		if (node.isRoot()) {
            return 1; // root has no rate
        }
		synchronized (this) {
    		if (recompute) {
    			calculateUnscaledBranchRates();
                recalculateScaleFactor();
                recompute = false;
			}
        }
		double r = unscaledBranchRates[getNr(node)] * scaleFactor;
		System.out.println("rate = "+r+" nr="+node.getNr()+"  newNr="+getNr(node));
        return 1;	
     
	}
	
	private void recalculateScaleFactor() {
        if (normalize) {
            double timeTotal = 0.0;
            double branchTotal = 0.0;

            for (int i = 0; i < tree.getNodeCount(); i++) {
                Node node = tree.getNode(i);
                if (!node.isRoot()) {
                    double branchInTime = node.getLength();
                    double branchLength = branchInTime * unscaledBranchRates[node.getNr()];
                    timeTotal += branchInTime;
                    branchTotal += branchLength;
                }
            }

            scaleFactor = timeTotal / branchTotal;

            scaleFactor *= meanRate.getValue();
        } else {
            scaleFactor = 1.0;
        }
        // scaleFactor *= meanRate.getValue();
    }
	
    
	private void calculateUnscaledBranchRates() {
		if (this.useCategories)
			calculateUnscaledRatesForCategories();
		else
			calculateUnscaledRatesForRates();
	}
    
    
    
    private void calculateUnscaledRatesForCategories() {
    	double rate=1.0;
    	ContinuousDistribution gammaDist;
    	for(int i=0; i < tree.getNodeCount();i++) {
    		final Node node = tree.getNode(i);
    		final int nr = getNr(node);
    		if (! node.isRoot()) {   		
    			int category = categories.getValue(nr);
    			System.out.println("cat = "+category+ " node length="+node.getLength());
    			System.out.println("mu = "+ratesMean.getValue()+ " omega="+ratesOmega.getValue());
    			
    			final double beta = node.getLength() / ratesOmega.getValue();
    			final double alpha = ratesMean.getValue() * beta;
    			gammaDist = new GammaDistributionImpl(alpha, beta); 
    			final double p = (category + 0.5) / numCategories;
    			System.out.println("alpha = "+alpha+ " beta="+beta+"  p="+p);
    			try {
    				rate = gammaDist.inverseCumulativeProbability( p );
    			} catch (MathException e) {
    				throw new RuntimeException("Failed to compute inverse cumulative probability");
    			}
    		}
    		System.out.println(" rate = "+rate);
    		unscaledBranchRates[nr] = rate;
    	}   	
    }
    
    private void calculateUnscaledRatesForRates() {
    	double rate=1.0;
    	ContinuousDistribution gammaDist;
    	for(int i=0; i < tree.getNodeCount();i++) {
    		final Node node = tree.getNode(i);
    		final int nr = getNr(node);
    		if (! node.isRoot()) {   		
    			double p = rateProbs.getValue(nr);
    	        final double beta = node.getLength() / ratesOmega.getValue();
    	        final double alpha = ratesMean.getValue() * beta;
    	        gammaDist = new GammaDistributionImpl(alpha, beta); 
    	        try {
    	        	rate = gammaDist.inverseCumulativeProbability( p );
    	        } catch (MathException e) {
    	        	throw new RuntimeException("Failed to compute inverse cumulative probability");
    	        }
    		}
    		unscaledBranchRates[nr] = rate;
    	}   	
    }
    
    
	
    /*  computes the idx used to access the arrays a[branchCount] 
     * needed in case of node swaps involving the root */
    private int getNr(Node node) {
        int nodeNr = node.getNr();
        if (nodeNr > tree.getRoot().getNr()) {
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
        recompute = true;
        recomputeScaleFactor = true;

        return recompute;
    }
    
    @Override
    public void store() {
    	recompute = true;
        super.store();
    }

    @Override
    public void restore() {
    	recompute = true;
        super.restore();
    }
   
    
	
}
