package util;

import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.core.parameter.RealParameter;
import beast.math.distributions.LogNormalDistributionModel;

public class DistTest {

	public static void timeDistributions() {
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
    
    public static void timeLogNormal(double start, double end, double delta, int numIter) {
    	
    	LogNormalDistributionModel distLN = new LogNormalDistributionModel();
    	Double[] M = new Double[1]; M[0] = 1.0;
    	RealParameter pM = new RealParameter(M);
    	Double[] S = new Double[1]; S[0] = 0.5;
    	RealParameter pS = new RealParameter(S);
    	
    	distLN.MParameterInput.setValue(pM, distLN);
    	distLN.SParameterInput.setValue(pS, distLN);
    	   	
    	/* LogNormal */	
    	Timer timer = new Timer();
    	timer.startTime();
		  	
    	double p = start;
    	for(int i=0; i < numIter; i++) {
    		try {
    			final double x = distLN.inverseCumulativeProbability(p);
    		} catch (Exception e) {
    			throw new RuntimeException("problem with inverse cum prob");
    		}
    		p += delta;
    	}
    	
    	timer.endTime(true);
				
		System.out.println("Time log normal:"+timer.duration);
    }
    
    public static void timeGamma(double start, double end, double delta, int numIter) {
    	Timer timer = new Timer();
    	timer.startTime();
		
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
    	
    	timer.endTime(true);
		
		System.out.println("Time gamma:"+ timer.duration);
    	
    }
	
	

}
