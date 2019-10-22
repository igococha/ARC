package beast.evolution.branchratemodel;


import java.util.Arrays;

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
import beast.evolution.branchratemodel.BranchRateModel.Base;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.LogNormalDistributionModel;
import beast.math.distributions.LogNormalDistributionModel.LogNormalImpl;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

/**
 * @author Igor Siveroni
 */

@Description("Implements ARC: th Additive Uncorrelated Ralexed Clock Model")
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
    
    
    private boolean recompute = true;
    private boolean categoriesOnly = false;
    
    
    private boolean treeOnly = false;
    private double[] branchLengths; // accessed with nr
    private double[] storedBranchLengths; 
    
    double[] rates;
    double[] storedRates;
    
    private boolean recomputeScaleFactor = true;
    private double scaleFactor = 1.0; 
    private double storedScaleFactor = 1.0;
    
    private int numCalls=0;
    private long startTime, endTime, duration;
    

    @Override
	public void initAndValidate() {
		tree = treeInput.get();
		nodeCount = tree.getNodeCount();
		branchCount = nodeCount - 1;
		categories = categoryInput.get();
		rateProbs = rateProbsInput.get();
		
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
        // unscaledBranchRates = new double[tree.getNodeCount()];
        rates = new double[branchCount];
        storedRates = new double[branchCount];
        branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];
        for(int i=0; i < nodeCount; i++) {
        	branchLengths[i]=-1;
        }
		// testGamma(9.0, 2.0);
        //timeDistributions();
	}
    
 
	

	//get the rate for node: R = r*scale*meanRate
	@Override
	public double getRateForBranch(Node node) {
		startTime();
		//System.out.println("Getting branch");
		if (node.isRoot()) {
			//System.out.println("root");
            return 1; // root has no rate
        }
		synchronized (this) {
    		if (recompute) {
    			calculateUnscaledBranchRates();
    			if (normalize)
    				recalculateScaleFactor();
    			recompute = false;
			}
        }
		double r = rates[getNr(node)] * scaleFactor;
		//System.out.println("rate = "+r+" nr="+node.getNr()+"  newNr="+getNr(node));
		endTime(false);
        return 1;	
     
	}
	
	private void recalculateScaleFactor() {
		System.out.println("recalculate scale factor");
        if (normalize) {
            double timeTotal = 0.0;
            double branchTotal = 0.0;

            for (int i = 0; i < tree.getNodeCount(); i++) {
                Node node = tree.getNode(i);
                if (!node.isRoot()) {
                    double branchInTime = node.getLength();
                    double branchLength = branchInTime * rates[node.getNr()];
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
		//System.out.println("recalculate branch rates");
		//debugState();
		//debugTree();
		if (this.useCategories) {
			calculateRatesForCategories();
		} else {
			calculateRatesForRates();
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
    
    
    private void calculateRatesForCategories() {
    	
    	double rate=1.0;
    	ContinuousDistribution gammaDist;
    	for(int i=0; i < tree.getNodeCount();i++) {
    		final Node node = tree.getNode(i);
    		if (! node.isRoot()) {   			
    			final int nr = getNr(node);
    			if (categoriesOnly && !categories.isDirty(nr)) {
    				continue;
    			}
    			final double l = node.getLength();  
    			if (treeOnly) {
    				// it seems that checking for node.isDirty==Tree.IS_CLEAN is not enough
    				final double diff = l - branchLengths[node.getNr()];
    				final boolean isDiff = Math.abs(diff) > 0.00001;
    				if (!isDiff) {
    					continue;
    					//System.out.println("skip node");
    				}
    				//System.out.println("nr = "+nr+" l ="+l+" "+branchLengths[node.getNr()]);   				
    			}
    			final int category = categories.getValue(nr);			
    			final double theta = ratesOmega.getValue() /  l  ;
    			final double k = ratesMean.getValue() / theta;
    			gammaDist = new GammaDistributionImpl(k, theta); 
    			final double p = (categorr  + 0.5) / numCategories;
    			branchLengths[node.getNr()] = l;
    			numCalls++;
    			//System.out.println("k = "+k+ " theta="+theta+"  p="+p);
    			try {
    				rate = gammaDist.inverseCumulativeProbability( p );
    			} catch (MathException e) {
    				throw new RuntimeException("Failed to compute inverse cumulative probability");
    			}
    			rates[nr] = rate;  			
    		} else {  // node is root
    			branchLengths[node.getNr()] = node.getLength(); 
    		}
    		// System.out.println(" rate = "+rate);
    	}
    	categoriesOnly=false;
    	treeOnly = false;
    }
    
    private void calculateRatesForRates() {
    	double rate=1.0;
    	ContinuousDistribution gammaDist;
    	for(int i=0; i < tree.getNodeCount();i++) {
    		final Node node = tree.getNode(i);
    		if (! node.isRoot()) {   
    			final int nr = getNr(node);
    			double p = rateProbs.getValue(nr);
    			final double theta = ratesOmega.getValue() /  node.getLength()   ;
    			final double k = ratesMean.getValue() / theta;
    			gammaDist = new GammaDistributionImpl(k, theta);   			  	       
    	        try {
    	        	rate = gammaDist.inverseCumulativeProbability( p );
    	        } catch (MathException e) {
    	        	throw new RuntimeException("Failed to compute inverse cumulative probability");
    	        }
    	        rates[nr] = rate;
    		}
    	}   	
    }
    
    
	
    /*  computes the idx used to access the arrays a[branchCount] 
     * needed in case of node swaps involving the root */
    private int getNr(Node node) {
        int nodeNr = node.getNr();
        if (nodeNr > tree.getRoot().getNr()) {
        	System.out.println("-- if root swap");
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
        recomputeScaleFactor = true;
        categoriesOnly=true;
        treeOnly= true;
        //System.out.println("calling requires recalculation");
        if (treeInput.get().somethingIsDirty()) {
        	//System.out.println("tree changes");
        	categoriesOnly=false;
        	recompute = true;
        	//return true;
        }
        
        if (ratesMeanInput.get().isDirtyCalculation()) {
        	//System.out.println("Gamma mu changes");
        	categoriesOnly=false;
        	treeOnly=false;
        	recompute = true;
        	//return true;
        }
     
        if (ratesOmegaInput.get().isDirtyCalculation()) {
        	//System.out.println("Gamma omega changes");
        	categoriesOnly=false;
        	treeOnly=false;
        	recompute = true;
        	//return true;
        }
     
        if (categoryInput.get() != null && categoryInput.get().somethingIsDirty()) {
        	//System.out.println("category changes");
        	treeOnly=false;
        	recompute = true;  // come back to this
        	//return true;
        }
    
        if (rateProbsInput.get() != null && rateProbsInput.get().somethingIsDirty()) {
        	//System.out.println("rate probs changes");
        	categoriesOnly=false;
        	treeOnly=false;
        	recompute = true;  // come back to this
        	//return true;
        }
    

        // still don't know the use of this rate
        if (meanRate.somethingIsDirty()) {
        	return true;
        }     
        
       
        
        return recompute;
    }
    
    @Override
    public void store() {
    	System.arraycopy(rates, 0, storedRates, 0, rates.length);
    	//System.arraycopy(branchLengths, 0, branchLengths, 0, branchLengths.length);
    	storedScaleFactor = scaleFactor;
        super.store();
        System.out.println("num calls="+numCalls+"  duration="+duration);
    }

    @Override
    public void restore() {
    	//System.out.println("restore branches");
    	double[] tmp = rates;
        rates = storedRates;
        storedRates = tmp;
        //tmp = branchLengths;
        //branchLengths = storedBranchLengths;
        //storedBranchLengths = tmp;
        scaleFactor = storedScaleFactor;
        super.restore();
    }
   
    /* piggybacking some tests */
    
    private void testGamma(double alpha,double beta) {
    	double k = alpha; 
    	double theta = 1/beta;
    	// apache commons implementations corresponds to Gamma(k,theta)
    	ContinuousDistribution gammaDist = new GammaDistributionImpl(k, theta); 
    	for (double p=0; p < 1.0; p += 0.1) {
    		double pinv;
    		try {
    			pinv = gammaDist.inverseCumulativeProbability( p );
    		} catch (Exception e) {
    			pinv = Double.NaN;
    		}
    		System.out.println(p + " --> "+pinv);
    	}
    	  	
    }
    
    private void timeDistributions() {
    	final int numIter=1000;
    	final double start=0.1;
    	final double end=0.9;
    	final double delta=(end-start)/numIter;
    	
    	timeLogNormal(start,end,delta,numIter);
    	
    	timeGamma(start,end,delta,numIter);
    	
    	timeLogNormal(start,end,delta,numIter);
    	
    	timeGamma(start,end,delta,numIter);
    	
    	timeLogNormal(start,end,delta,numIter);
    	
    	timeGamma(start,end,delta,numIter);
    	
    	timeLogNormal(start,end,delta,numIter);
    	
    	timeGamma(start,end,delta,numIter);
    	
    	timeLogNormal(start,end,delta,numIter);
    	
    	timeGamma(start,end,delta,numIter);
    	
    	/* time ln = 3/4 time gamma */
    	
    }
    
    private void timeLogNormal(double start, double end, double delta, int numIter) {
    	
    	LogNormalDistributionModel distLN = new LogNormalDistributionModel();
    	Double[] M = new Double[1]; M[0] = 1.0;
    	RealParameter pM = new RealParameter(M);
    	Double[] S = new Double[1]; S[0] = 0.5;
    	RealParameter pS = new RealParameter(S);
    	
    	distLN.MParameterInput.setValue(pM, distLN);
    	distLN.SParameterInput.setValue(pS, distLN);
    	   	
    	/* LogNormal */	
    	
    	startTime();
		  	
    	double p = start;
    	for(int i=0; i < numIter; i++) {
    		try {
    			final double x = distLN.inverseCumulativeProbability(p);
    		} catch (Exception e) {
    			throw new RuntimeException("problem with inverse cum prob");
    		}
    		p += delta;
    	}
    	
    	endTime(true);
		
		
		System.out.println("Time log normal:"+duration);
    }
    
    private void timeGamma(double start, double end, double delta, int numIter) {
    	
    	startTime();
		
    	GammaDistribution g;
    	ContinuousDistribution gammaDist = new GammaDistributionImpl(2, 0.5); 
    	
    	double p = start;
    	for(int i=0; i < numIter; i++) {
    		try {
    			final double x = gammaDist.inverseCumulativeProbability(p);
    		} catch (Exception e) {
    			throw new RuntimeException("problem with inverse cum prob");
    		}
    		p += delta;
    	}
    	
    	endTime(true);
		
		System.out.println("Time gamma:"+duration);
    	
    }
    
    private void startTime() {
    	startTime = System.nanoTime();
    }
    
    private void endTime(boolean reset) {
    	endTime = System.nanoTime();
    	if (reset)
    		duration = (endTime - startTime)/100000;
    	else
    		duration += (endTime - startTime)/100000;
    }
    
    
    
    
	
}
