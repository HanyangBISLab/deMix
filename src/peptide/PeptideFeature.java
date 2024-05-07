package kr.ac.hanyang.bislab.demix.peptide;

import java.util.ArrayList;

public class PeptideFeature {

	int startScan;
	int endScan;
	int weightedCenter;
	
	double startRT;
	double endRT;
	double centerRT;
	
	private int bestFittingScan;
	private double bestFittingRT;
	
	double weightedAverMZ;
	
	ArrayList<Double> aggregatedCluster;
	
	public PeptideFeature(int st, double stR, int ed, double edR, int ct, double ctR, int best, double beR, double weiMz) {
		
		startScan		= st;
		endScan			= ed;
		weightedCenter	= ct;
		
		startRT			= stR;
		endRT			= edR;
		centerRT		= ctR;

		bestFittingScan	= best;
		bestFittingRT	= beR;
		
		weightedAverMZ	= weiMz;
	}
	
	public void setIsotopeCluster(double[] peaks) {
		
		aggregatedCluster= new ArrayList<Double>();
		
		double sum= 0.;
		for(int i=0; i<peaks.length; i++) {
		//	aggregatedCluster.add(peaks[i]);
			sum += peaks[i];
		}
		
		if( sum != 0 ) {
			for(int i=0; i<peaks.length; i++) {
				aggregatedCluster.add(peaks[i]/sum);
			/*	double norm= peaks[i]/sum;
				if( norm < 0.01 ) aggregatedCluster.add(0.);
				else aggregatedCluster.add(norm);//*/
			}
		}
		else {
			for(int i=0; i<peaks.length; i++) {
				aggregatedCluster.add(0.);
			}
		}
	}
	
	public String toString(double pmz) {
		return String.format("%.3f\t%.3f\t%d\t%d\t%d\t%.2f", weightedAverMZ, weightedAverMZ-pmz, startScan, endScan, bestFittingScan, bestFittingRT/60);
	}
	
	public String toOutput() {
		StringBuffer buf= new StringBuffer( String.format("%d\t%d\t%d\t%.2f\t", startScan, endScan, bestFittingScan, bestFittingRT/60) );
		
		if( aggregatedCluster.size() == 0 ) buf.append('-');
		else {
			
			if( aggregatedCluster.get(0) == 0. ) buf.append("0.0");
			else buf.append(String.format("%.3f", aggregatedCluster.get(0)));
			
			for(int i=1; i<aggregatedCluster.size(); i++ ) {
				if( aggregatedCluster.get(i) == 0. ) buf.append(";0.0");
				else buf.append(String.format(";%.3f", aggregatedCluster.get(i)));
			}
		}
		
		return buf.toString();
	}
	
	public int getStartScan() { return startScan; }
	public int getEndScan() { return endScan; }
	
	public double getStartRT() { return startRT; }
	public double getEndRT() { return endRT; }
	public double getCenterRT() { return centerRT; }
	public double getBestRT() { return bestFittingRT; }
	
	public double getWeightedAverMZ() { return weightedAverMZ; }

	public ArrayList<Double> getIsotopeCluster() { return aggregatedCluster; }
	

}












