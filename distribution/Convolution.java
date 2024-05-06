package kr.ac.hanyang.bislab.demix.distribution;

import java.util.ArrayList;
//convolution two arrays

public class Convolution {
	
	private ArrayList<Double> convResult;
	private double convIntSum;
	private double maxInt;
	private int maxIndex;
	
	public Convolution(ArrayList<Double>x, ArrayList<Double>y){//convolution {1,2,3), {4,5,6}
		
		convResult 	= new ArrayList<Double>();
		convIntSum 	= 0;
		maxInt 		= 0;
		maxIndex 	= 0;
		
		int len_x = x.size();
		int len_y = y.size();
		int len_x_y = 2*len_x + len_y - 2;
		double[] a = new double[len_x_y];			
		double[] b = new double[len_x_y];			
		double[] c = new double[len_x+len_y-1];
		for(int i=0 ; i<len_x ; i++)	a[i] = x.get(len_x-i-1);//3,2,1,0,0,0,0
		for(int i=len_x ; i<len_x+len_y-1 ; i++) a[i] = 0;
		for(int i=0 ; i<len_x-1 ; i++) 	b[i] = 0;			//0,0,4,5,6,0,0
		for(int i=0 ; i<len_y ; i++) 	b[len_x-1+i] = y.get(i);
		for(int i=0 ; i<len_x-1 ; i++) 	b[len_x-1+len_y+i] = 0;
		for(int i=0 ; i<len_x+len_y-1; i++){				//convolution
			for(int j=0 ; j<len_x_y ; j++)  c[i] += a[j] * b[j];
			a = shiftArray(a);
		}
		
		for(int i=0;i<c.length;i++){
		//	convResult.add(c[i]);
			convIntSum +=c[i];
			if(c[i] > maxInt) {
				maxInt = c[i];
				maxIndex = i;
			}
		}
		
		for(int i=0; i<c.length; i++) {
			convResult.add(c[i]/convIntSum);
		}
		
	}
	
	private double[] shiftArray(double[]a){//right shift array and fill with 0 
		double[] c = new double[a.length]; //{3,2,1,0,0,0,0} ->{0,3,2,1,0,0,0}
		for(int i=0 ; i<a.length-1 ; i++) c[i+1] = a[i];
		return c;
	}
	
	public ArrayList<Double> getConvResult() { return convResult; }
	public double getConvIntSum() { return convIntSum; }
	public double getMaxInt() { return maxInt; }
	public int getMaxIndex() { return maxIndex; }
	
}

