package beast.evolution.substitutionmodel;

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

public class ARClockModel extends Base {
	final public Input<ParametricDistribution> rateDistInput = new Input<>("distr", "the distribution governing the rates among branches. Must have mean of 1. The clock.rate parameter can be used to change the mean rate.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> categoryInput = new Input<>("rateCategories", "the rate categories associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    final public Input<Integer> numberOfDiscreteRates = new Input<>("numberOfDiscreteRates", "the number of discrete rate categories to approximate the rate distribution by. A value <= 0 will cause the number of categories to be set equal to the number of branches in the tree. (default = -1)", -1);
    final public Input<RealParameter> rateProbsInput = new Input<>("rateprobs", "the accum probabilities of the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.XOR, categoryInput);
    final public Input<Tree> treeInput = new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
    final public Input<Boolean> normalizeInput = new Input<>("normalize", "Whether to normalize the average rate (default false).", false);
    // there are three modes to represent the rates on the branches
    enum Mode {
        categories,
        rates
    }
	
    boolean useCategories;
    // either categories or rateProbsParameter is used
    RealParameter rateProbsParameter; //when mode=rates
    IntegerParameter categories; //when mode=categories
    
    ParametricDistribution distribution; //the distribution of the rates
    // i.e. LogNormal(M,S,MeanInRealSpace), then get the three parameters
    double M;double S;boolean MeanInRealSpace;

    RealParameter meanRate;
    Tree tree;
    private int branchCount;//the number of branches of the tree
    private boolean normalize = false;//
    private boolean recompute = true;//distribution.inverseCumulativeProbability((category + 0.5) / rates.length);
    private boolean renormalize = true;//
    private double[] rates; //the output rates
    private double[] storedRates; //
    private double scaleFactor = 1.0; //initial
    private double storedScaleFactor = 1.0; //initial
    int numCategories = 100;//


	public ARClockModel() {
		tree = treeInput.get();
		branchCount = tree.getNodeCount() - 1;
		categories = categoryInput.get();
		rateProbsParameter = rateProbsInput.get();
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
        	if (rateProbsParameter.getDimension() != branchCount) {
        		rateProbsParameter.setDimension(branchCount);
        		//randomly draw rates from the lognormal distribution
        		Double [] initialProbs = new Double[branchCount];
        		for (int i = 0; i < branchCount; i++) {
        			initialProbs[i] =Randomizer.nextDouble();  // [0,1)
        		}
        		RealParameter other = new RealParameter(initialProbs);
        		rateProbsParameter.assignFromWithoutID(other);
        	}
        	rateProbsParameter.setLower(0.0);
        	rateProbsParameter.setUpper(1.0);
        }
        
        if (useCategories) {
            // rates are initially zero and are computed by getRawRate(int i) as needed
            rates = new double[numCategories];
            storedRates = new double[numCategories];
        }
        normalize = normalizeInput.get();
        meanRate = meanRateInput.get();
        if (meanRate == null) {
            meanRate = new RealParameter("1.0");
        }
        try { // igor: should this be ONE in ARC's case?
            double mean = rateDistInput.get().getMean();
            if (Math.abs(mean - 1.0) > 1e-6) {
                Log.warning.println("WARNING: mean of distribution for additive relaxed clock model is not 1.0.");
            }
        } catch (RuntimeException e) {
            // ignore
        }
        
		
	}

	@Override
	public double getRateForBranch(Node node) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub

	}

}
