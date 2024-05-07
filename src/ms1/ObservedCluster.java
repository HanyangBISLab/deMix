package kr.ac.hanyang.bislab.demix.ms1;

public class ObservedCluster {
	
	int	   scanNo;
	double monoMz;
	double amount;
	double shapeScore;
	double pif;
	
	public ObservedCluster(int sn, double m, double a, double s, double p) {
		scanNo= sn;
		monoMz= m;;
		amount= a;
		shapeScore= s;
		pif= p;
	}
	public int getScanNo() {
		return scanNo;
	}
	public double getWeightedScore() {
		return shapeScore*pif;
	}
	public double getWeightedMZ() {
		return monoMz*amount;
	}
	
	
	
}
