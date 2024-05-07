package kr.ac.hanyang.bislab.demix.distribution;

import java.util.ArrayList;

public class DDMatch implements Comparable<DDMatch> {

	int 	deuterated;
	int 	deu2nd;
	
	int 	mpc;
	double 	mpc_err;
	double  scailing_factor;
	
	double 	mainPortion;
	double 	interperted;
	double 	error;
	
	int 	side;//0:leftmatched, 1:rightmatched
	
	ArrayList<Double> predictedCluster;
	
	public DDMatch() {
	}
	
	public DDMatch(int deu, int m, double merr, double sc, double inted, double err, int fs) {
		
		deuterated= deu;
		deu2nd= -1;
		
		mpc= m;
		mpc_err= merr;

		scailing_factor= sc;
		interperted= inted;
		error= err;
		
		side= fs;
	}
	
	public int getDeuteratedNum() { return deuterated; }
	public double getMPC() { return mpc; }
	public double getMPCError() { return mpc_err; }
	public String toString() { 
		StringBuffer buf= new StringBuffer();
		if( deu2nd == -1 ) {
			buf.append(String.format("%d/1.00\t-", deuterated));
		}
		else buf.append(String.format("%d/%.2f\t%d/%.2f", deuterated, mainPortion, deu2nd, 1-mainPortion));
		
		buf.append( String.format("\t%.2f\t%d", interperted, mpc) );
		
		return buf.toString();
	}
	
	public String toOutput(int hdxSites) { 
		StringBuffer buf= new StringBuffer(getStrHDXRate(hdxSites));
		if( deu2nd == -1 ) {
			buf.append(String.format("\t%d\t1.00\t-\t-\t", deuterated));
		}
		else buf.append(String.format("\t%d\t%.2f\t%d\t%.2f\t", deuterated, mainPortion, deu2nd, 1-mainPortion));
				
		if( predictedCluster.size() == 0 ) buf.append('-');
		else {
			
			if( predictedCluster.get(0) == 0. ) buf.append("0.0");
			else buf.append(String.format("%.3f", scailing_factor*predictedCluster.get(0)));
			
			for(int i=1; i<predictedCluster.size(); i++ ) {
				if( predictedCluster.get(i) == 0. ) buf.append(";0.0");
				else buf.append(String.format(";%.3f", scailing_factor*predictedCluster.get(i)));
			}
		}
		
		return buf.toString();
	}
	public ArrayList<Double> getPedictedCluster() { 
		return predictedCluster;
	}
	public String getStrMatchQuality() { 
		return String.format("%.2f", interperted);
	}
	
	public String getStrHDXNumer() { 
		
		if( deu2nd == -1 ) 
			return String.valueOf(deuterated);
		else 
			return String.format("%.1f", (deuterated*mainPortion+deu2nd*(1-mainPortion)));

	}
	
	public String getStrHDXRate(int numSites) { 
		
		if( deu2nd == -1 ) 
			return String.format("%.1f", 100.*deuterated/numSites);
		else 
			return String.format("%.1f", 100.*(deuterated*mainPortion+deu2nd*(1-mainPortion))/numSites);
	}
	
	public void setPredictedCluster(ArrayList<Double> pc) {
		predictedCluster= pc;
	}
	
	public void setDeuterated(int deu, double por, ArrayList<Double> pc) {
		deuterated= deu;
		mainPortion= por;
		deu2nd= -1;
		predictedCluster= pc;
	}
	
	public void setDeuterated(int deu, double por, int d2, ArrayList<Double> pc) {
		deuterated= deu;
		mainPortion= por;
		deu2nd= d2;
		predictedCluster= pc;
	}
	
	public boolean isExplainedEnough(double cut) {
		return cut <= interperted;
	}
	
	
	public int compareTo(DDMatch x){
		
		if( mpc < x.mpc ) return 1;
		else if( mpc > x.mpc ) return -1;
		
		if( interperted < x.interperted ) return 1;
		else if( interperted > x.interperted ) return -1;
		
		if( mpc_err > x.mpc_err ) return 1;
		else if( mpc_err < x.mpc_err ) return -1;
		return 0;	
		
	}

	
	
	
	
	
	
	
}
