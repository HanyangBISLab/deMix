package kr.ac.hanyang.bislab.demix.distribution;

public class Binomial {
	
	public static double logsum1(int n, int k){
		double value1=0;
		for (int i=1;i<=k; i++) value1 += Math.log(n+1-i);
		return value1;
    }
	
    public static double logsum2(int k){
    	double value2 = 0;
    	for (int i=1;i<=k;i++) value2 += Math.log(i);
    	return value2;
    }
    
    public static double pdf(int n, double p, int k){
    	if( p == 0 ){
    		if( k == 0 ) return 1;
    		else return 0;
    	}
    	if( p == 1 ){
    		if( k == n ) return 1;
    		else return 0;
    	}
    	
    	double logP= logsum1(n,k) - logsum2(k) + k*Math.log(p) + (n-k)*Math.log(1-p) ;
    	return Math.exp(logP);
    }
    
    
}







