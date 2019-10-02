package beast.evolution.substitutionmodel;

import java.util.Arrays;

import org.apache.commons.math.MathException;
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

@Description("Implements ARC: Additive Uncorrelated Ralezed Clock Model")
@Citation(value =
        "Didelot X and Volz EM (2019) Additive uncorrelated  relaxed clock models\n" +
                "  for the dating of epidemiological models.", DOI = "TBC",
        year = 2019, firstAuthorSurname = "Didelot")

public class ARClockModelTemplate extends Base {
	final public Input<ParametricDistribution> rateDistInput = new Input<>("distr", "the distribution governing the rates among branches. Must have mean of 1. The clock.rate parameter can be used to change the mean rate.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> categoryInput = new Input<>("rateCategories", "the rate categories associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    final public Input<Integer> numberOfDiscreteRates = new Input<>("numberOfDiscreteRates", "the number of discrete rate categories to approximate the rate distribution by. A value <= 0 will cause the number of categories to be set equal to the number of branches in the tree. (default = -1)", -1);
    final public Input<RealParameter> rateProbsInput = new Input<>("rateprobs", "the accum probabilities of the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.XOR, categoryInput);
    final public Input<Tree> treeInput = new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
    final public Input<Boolean> normalizeInput = new Input<>("normalize", "Whether to normalize the average rate (default false).", false);
   
	
    boolean useCategories;
    // either categories or rateProbsParameter is used
    RealParameter rateProbs; //when mode=rates
    IntegerParameter categories; //when mode=categories
    
    ParametricDistribution distribution; //the distribution of the rates
    // i.e. LogNormal(M,S,MeanInRealSpace), then get the three parameters
    double M;double S;boolean MeanInRealSpace;
    
    GammaDistribution m_dist = new GammaDistributionImpl(1, 1);

    RealParameter meanRate;
    Tree tree;
    private int branchCount;//the number of branches of the tree
    private boolean normalize = false;//   
    private boolean renormalize = true;
    
    private boolean recompute = true;
    private double[] ratesxCategory; //the output rates
    private double[] storedRatesxCategory; //
    private double scaleFactor = 1.0; //initial
    private double storedScaleFactor = 1.0; //initial
    int numCategories = 100;//


    @Override
	public void initAndValidate() {
		tree = treeInput.get();
		branchCount = tree.getNodeCount() - 1;
		categories = categoryInput.get();
		rateProbs = rateProbsInput.get();
		distribution = rateDistInput.get();
		
		//get the mean and standard deviation of lognormal distribution
        if(distribution instanceof LogNormalDistributionModel){
            LogNormalDistributionModel mylognormal=(LogNormalDistributionModel)distribution;
             M=mylognormal.MParameterInput.get().getValue();
             S=mylognormal.SParameterInput.get().getValue();
             MeanInRealSpace=mylognormal.hasMeanInRealSpaceInput.get();
        } else {
        	throw new IllegalArgumentException("ARClock model only implemented for LogNormal distribution");
        }
        
        if (categories==null) useCategories=true;
        else useCategories = true;
        
        if (useCategories) {
            numCategories = numberOfDiscreteRates.get();
            if (numCategories <= 0) numCategories = branchCount;
            Log.info.println("  ARCClockModel: using " + numCategories + " rate " +
                    "categories to approximate rate distribution across branches.");
        } else {
            if (numberOfDiscreteRates.get() != -1) {  // if it's not the default i.e. it has been set as input
                throw new IllegalArgumentException("Can't specify both numberOfDiscreteRates and  rates inputs.");
            }
            Log.info.println("  ARClockModel: using rate probs as parameter for sampling  across branches.");          
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
        
        if (useCategories) {
            // rates are initially zero and are computed by getRawRate(int i) as needed
            ratesxCategory = new double[numCategories];
            storedRatesxCategory = new double[numCategories];
        }
        normalize = normalizeInput.get();
        meanRate = meanRateInput.get();
        if (meanRate == null) {
            meanRate = new RealParameter("1.0");
        }
        if (Math.abs(distribution.getMean()-1.0) > 1e-6 ) { // igor: should this be ONE in ARC's case?
        	Log.warning.println("WARNING: mean of distribution for additive relaxed clock model is not 1.0.");
        }
       
        m_dist = new GammaDistributionImpl(1, 1);

		
	}
	
    private void initialize() {
        if (useCategories) {
           Arrays.fill(ratesxCategory, 0.0);
        }
    }
	
	/**
     * Computes a scale factor for normalization. Only called if normalize=true.
     */
    private void computeScaleFactor() {
        //scale mean rate to 1.0 or separate parameter
        double treeRate = 0.0;
        double treeTime = 0.0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);
            if (!node.isRoot()) {
                treeRate += getRawRate(node) * node.getLength();
                treeTime += node.getLength();
            }
        }
        scaleFactor = 1.0 / (treeRate / treeTime);
    }

	//get the rate for node: R = r*scale*meanRate
	@Override
	public double getRateForBranch(Node node) {
		if (node.isRoot()) {
            return 1; // root has no rate
        }
        if (recompute) {
           // avoid thread clash
            synchronized (this) {
                initialize();
                recompute = false;
            }
        }
        if (renormalize) {
            if (normalize) {
                synchronized (this) {
                    computeScaleFactor();
                }
            }
            renormalize = false;
        }
        return getRawRate(node) * scaleFactor * meanRate.getValue();
	}
	
    private double getRawRate(Node node) {
    	if (useCategories) {
    		return getRawRateForCategories(node);
    	} else {
    		return getRawRateForRates(node);
    	}
    }
    
    private double getRawRateForCategories(Node node) {
        int nodeNumber = node.getNr();
        if (nodeNumber == branchCount) {
            // root node has nr less than #categories, so use that nr
            nodeNumber = node.getTree().getRoot().getNr();
        }
        int category = categories.getValue(nodeNumber);
        if (ratesxCategory[category] == 0.0) {
            try {
                ratesxCategory[category] = distribution.inverseCumulativeProbability((category + 0.5) / numCategories);
            } catch (MathException e) {
                throw new RuntimeException("Failed to compute inverse cumulative probability");
            }
        }
        return ratesxCategory[category];
    }
    
    private double getRawRateForRates(Node node) {
    	int nodeNumber = node.getNr();
        if (nodeNumber == branchCount) {
            // root node has nr less than #categories, so use that nr
            nodeNumber = node.getTree().getRoot().getNr();
        }
        double p = rateProbs.getValue(nodeNumber);
        double rate;
        try {
            rate = distribution.inverseCumulativeProbability(p);
        } catch (MathException e) {
            throw new RuntimeException("Failed to compute inverse cumulative probability");
        }
    	return rate;
    }
	
    
    // Need to review this
    @Override
    protected boolean requiresRecalculation() {
        recompute = false;
        renormalize = true;

//        if (treeInput.get().somethingIsDirty()) {
//        	recompute = true;
//            return true;
//        }
        
        // rateDistInput cannot be dirty?!?
        if (rateDistInput.get().isDirtyCalculation()) {
            recompute = true;
            return true;
        }
        // NOT processed as trait on the tree, so DO mark as dirty
        if (categoryInput.get() != null && categoryInput.get().somethingIsDirty()) {
            //recompute = true;
            return true;
        }

       
        if (meanRate.somethingIsDirty()) {
            return true;
        }

        return recompute;
    }
    
    @Override
    public void store() {
    	if (useCategories) {
           System.arraycopy(ratesxCategory, 0, storedRatesxCategory, 0, numCategories);
            
        }
        storedScaleFactor = scaleFactor;
        super.store();
    }

    @Override
    public void restore() {
        if(useCategories){
            double[] tmp = ratesxCategory;
            ratesxCategory = storedRatesxCategory;
            storedRatesxCategory = tmp;
        }
        scaleFactor = storedScaleFactor;
        super.restore();
    }
   
    
	
}
