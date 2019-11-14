package util;

public class Timer {
	
	public long startTime, endTime, duration;
	
	public void startTime() {
    	startTime = System.nanoTime();
    }
    
    public void endTime(boolean reset) {
    	endTime = System.nanoTime();
    	if (reset)
    		duration = (endTime - startTime)/100000;
    	else
    		duration += (endTime - startTime)/100000;
    }
    

}
