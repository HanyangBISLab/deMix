package kr.ac.hanyang.bislab.demix.peptide;

import java.util.ArrayList;

import kr.ac.hanyang.bislab.demix.protein.Match2Protein;
import kr.ac.hanyang.bislab.demix.Params;

public class Peptide implements Comparable<Peptide> {
	
	int		index;
	String 	peptide;
	String 	stripSeq;
	double 	moleWeight;

	int 	charge;
	double 	pmz;
	double  peakWidth;
	double  theo_pmz;
	
	double 	rt;
	int 	scanNo;
	
	ArrayList<Double> 	theoriticalIsotopes;
	ArrayList<Double> 	trimmedIsotopes;
	PeptideFeature 		ms1Feature;
		
	int					numHDXSite; 
	
	Match2Protein 		protMatch;
	
	public Peptide(int num, String inPept, int cs, double p, double w) {
		
		index= num;
		
		peptide= inPept;	
		charge	= cs;
		pmz		= p;
		peakWidth= w;
		
		theo_pmz = p;
		
		moleWeight= (pmz - Params.Proton)*charge;

		int startP= 0, endP= peptide.length();
		if( peptide.charAt(1)=='.' && peptide.charAt(peptide.length()-2)=='.' ) {
			startP= 2;
			endP -= 2;
		}
		
		StringBuffer buf= new StringBuffer();
		for(int i=startP; i<endP; i++) {
			if( Character.isLetter(peptide.charAt(i)) ) 
				buf.append(peptide.charAt(i));
		}
		stripSeq= buf.toString();
		
		numHDXSite = stripSeq.length()-1;
		for(int i=0 ; i<stripSeq.length() ; i++)
			if(stripSeq.charAt(i)=='P') numHDXSite--;

		theoriticalIsotopes = IsotopicCluster.get(stripSeq.toString(), charge);
		
		double natArea = 0.;
		for(int k=0; k<theoriticalIsotopes.size(); k++){
			natArea += theoriticalIsotopes.get(k);
		}
		
		double tpArea = 0;
		trimmedIsotopes = new ArrayList<Double>();								
		for(int k=0; k<theoriticalIsotopes.size(); k++){
			tpArea += theoriticalIsotopes.get(k);
			trimmedIsotopes.add( theoriticalIsotopes.get(k) );
			if( tpArea/natArea > 0.9 ) break;
		}

		protMatch= null;
		
	}

	public void     correctSystemBias(double ppmshift, double error) { 
		pmz += (pmz/1000000*ppmshift);
		peakWidth = error;
		ms1Feature = null;
	}
	
	public void     setIndex(int i) { index = i; }
	public void 	setMS1Feature(PeptideFeature pf) { ms1Feature= pf; }
	public void 	setProtMatch(Match2Protein pf)   { protMatch= pf; }
	
	public int 		getIndex() { return index; }
	public String 	getInputPeptide() { return peptide; }
	public String 	getStripSequence() { return stripSeq; }

	public double 	getPMZ() { return pmz; }
	public double 	getTheoreticalPMZ() { return theo_pmz; }
	public double 	getPeakWidth() { return peakWidth; }
	public int 		getCharge() { return charge; }
	public double 	getMoleWeight() { return moleWeight; }
	
	public int 		getHDXSiteNum() { return numHDXSite; }
	public PeptideFeature getMS1Feature() { return ms1Feature; }
	public Match2Protein 	getMatch2Protein() { return protMatch; }
	
	public ArrayList<Double> getIsotopeCluster() { return theoriticalIsotopes; }
	
	public String toStringOfIsotopeCluster() {
		StringBuffer buf= new StringBuffer();
		
		buf.append(String.format("%.3f", theoriticalIsotopes.get(0)/100.));		
		for(int i=1; i<theoriticalIsotopes.size(); i++ ) {
			buf.append(String.format(";%.3f", theoriticalIsotopes.get(i)/100.));
		}
		return buf.toString();
	}
	
	public ArrayList<Double> getTrimmedIsotopeCluster() { return trimmedIsotopes; }

	public int compareTo(Peptide x) {
		
		if( pmz > x.pmz ) return 1;
		else if( pmz < x.pmz ) return -1;
		
		return 0;
	}
	
	
	private static final double	Electron = 0.000549;
	private static final double	Hydrogen = 1.007825035;
	private static final double	Oxygen = 15.99491463;
	private static final double	Proton = Hydrogen-Electron;
	private static final double	H2O = Hydrogen*2 + Oxygen;	
	
	private static double aaMass[]={
		71.03711, 0, 103.00919, 115.02694, 129.04259, 
		147.06841, 57.02146, 137.05891, 113.08406, 0, 
		128.09496, 113.08406, 131.04049, 114.04293, 0, 
		97.05276, 128.05858, 156.10111, 87.03203, 101.04768,
		150.953636, 99.06841, 186.07931, 0, 163.06333, 0};

	public static double getPeptMass(String pept){
		double mass= 0.;
		for(int i=0; i<pept.length(); i++){			
			if( aaMass[pept.charAt(i)-'A'] < 1. )
				return -1;
			mass += aaMass[pept.charAt(i)-'A'];
		}
		return mass;
	}
	public static double getPeptMW(String pept){
		return getPeptMass(pept)+H2O;
	}
	public static double getPeptMZ(String pept, int cs){
		return (getPeptMW(pept)+cs*Proton)/cs;
	}
	public static String getPeptideSequence(String peptide){
		int dot = 0;
		int len= peptide.length();
		
		if( peptide.charAt(1) == '.' &&  peptide.charAt(len-2) == '.' ) dot = 2;

		StringBuffer sequence = new StringBuffer();
		for(int i=dot; i<len-dot; i++) {
			if( Character.isLetter(peptide.charAt(i)) )
				sequence.append(peptide.charAt(i));
		}
		return sequence.toString();
	}
	
}





