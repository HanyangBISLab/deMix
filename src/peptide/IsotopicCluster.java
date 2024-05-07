package kr.ac.hanyang.bislab.demix.peptide;

import java.util.ArrayList;

import kr.ac.hanyang.bislab.demix.distribution.Binomial;
import kr.ac.hanyang.bislab.demix.distribution.Convolution;

public class IsotopicCluster {

		
	public static ArrayList<Double> get(String pept, int charge){
		return get( getComponentElements(pept, charge) );
	}
	
	
	public static ArrayList<Double> get(int[] elements){

	/*	System.out.print("C"+elements[C]+" | ");
		System.out.print("H"+elements[H]+" | ");
		System.out.print("N"+elements[N]+" | ");
		System.out.print("O"+elements[O]+" | ");
		System.out.print("S"+elements[S]+" | ");
		System.out.println();//*/
		
		ArrayList<Double> cluster = new ArrayList<Double>();

		double T1 = getIthT(elements, 1);
		double T2 = getIthT(elements, 2);
		double T4 = getIthT(elements, 4);
		
		double I0 = getI0(elements);
		cluster.add(I0);
		
	//	System.out.println(T1 + " " + T2 + " " + T4);

		double area = I0;
		for(int Kth=1; ; Kth++){
			double Ik = getKthIntensity(Kth, I0, T1, T2, T4);

			area += Ik;
			if( Ik/area < 0.001 && Ik < I0 ) {
				area -= Ik;
				break;
			}
			cluster.add( Ik );
		}
		
		for(int i=0; i<cluster.size(); i++) {
			cluster.set(i, cluster.get(i)*100/area);
		}
		
		return cluster;
	}
	
	
	private static double getI0(int[] ele){
		double I0 = 1;
		for(int i=0; i<NoElements; i++){
			I0 *= Math.pow(IsotopeComposition[i][0], ele[i]);
		}
		return I0;
	}
	
	
	private static double getKthIntensity(int Kth, double I0, double T1, double T2, double T4){
		
		double Ik = 0.;
		int pTwo = Kth/2+1, pFour = Kth/4+1;

		for(int i=0; i<pFour; i++){
			for(int j=0; j<pTwo; j++){
				int k = Kth - i*4 - j*2;
				if( k < 0 ) break;
		
			//	Ik += (Math.pow(T4, i)*Math.pow(T2, j)*Math.pow(T1, k))/(getFactorial(i)*getFactorial(j)*getFactorial(k));
				
				double denorm = Math.pow(T4, i)*Math.pow(T2, j)*Math.pow(T1, k);				
				Ik += Math.exp(Math.log(denorm) - Math.log(getFactorial(i)) - Math.log(getFactorial(j)) - Math.log(getFactorial(k))) ;
			}
		}
	
		return I0*Ik;
	}
	
	private static double  getFactorial(int index){
		double  fac = 1.;
		for(int i=2; i<=index; i++)
			fac *= i;
		return fac;
	}
	
	private static double getIthT(int[] ele, int index){
		double T = 0;
		for(int i=0; i<NoElements; i++){
			T += ele[i]*IsotopeComposition[i][index]/IsotopeComposition[i][0];
		}
		return T;
	}
	
	public static int[] getComponentElements(String pept, int charge){
	//	System.out.println("@@@"+pept+"@@@");
		int[] elements = new int[NoElements];		
		for(int i=0; i<pept.length(); i++){
			for(int j=0; j<NoElements; j++){
				elements[j] += AminoAcid[pept.charAt(i)-'A'][j];
			}
		}
		elements[H] += 2+charge;
		elements[O] += 1;
		
		return elements;
	}	
	public static int[] getAminoAcidChemComponenta(char aa){
		return AminoAcid[aa-'A'];
	}
	
	public static ArrayList<Double> getEnvelope(String chemStr){
		
		String[] cnum= chemStr.split("_");
		int[] elements = new int[NoElements];	
		for(int i=0; i<elements.length; i++) {
			elements[i]= Integer.parseInt(cnum[i]);
		//	System.out.println(elements[i]);
		}
		
		double mass= 0.;
		for(int i=0; i<elements.length; i++) {
			mass += elements[i]*monoMass[i];
		}
		System.out.println(mass);
		
		return get(elements);
	}
	
	public static double getMass(String chemStr){
		
		double m= 0.;
		
		String[] cnum= chemStr.split("_");
		
		for(int i=0; i<monoMass.length; i++) {
			m += monoMass[i]*Integer.parseInt(cnum[i]);
		}
		
		return m;
	}
	
	private static int NoElements = 5; 
	private static int C=0, H=1, N=2, O=3, S=4;
	private static double[] monoMass = {12.0, 1.0078250321, 14.0030740052, 15.9949146221, 31.97207069}; 
	private static double[][] IsotopeComposition={

		//mono    , 1       , 2      , 3, 4
		//internet
		{0.9893  , 0.0107  , 0      , 0, 0     }, //C
		{0.999885, 0.000115, 0      , 0, 0     }, //H
		{0.99632 , 0.00368 , 0      , 0, 0     }, //N
		{0.99757 , 0.00038 , 0.00205, 0, 0     }, //O
		{0.9493  , 0.0076  , 0.0429 , 0, 0.0002}, //S//*/

	/*	//emass
		{0.988930, 0.011070, 0       , 0, 0     }, //C
		{0.99985 , 0.00015 , 0       , 0, 0     }, //H
		{0.996337, 0.003663, 0       , 0, 0     }, //N
		{0.997590, 0.000374, 0.002036, 0, 0     }, //O
		{0.9502  , 0.0075  , 0.0421  , 0, 0.0002}, //S//*/
		
		//???
	/*	{0.9893  ,0.0107  , 0      , 0, 0     }, //C
		{0.99985 , 0.00015 , 0      , 0, 0     }, //H
		{0.99632 , 0.00368 , 0      , 0, 0     }, //N
		{0.99757 , 0.00038 , 0.00205, 0, 0     }, //O
		{0.9493  , 0.0076  , 0.0429 , 0, 0.0002}, //S//*/
	};
	private static int[][] AminoAcid={
		//C, H, N, O, S
		{ 3, 5, 1, 1, 0}, //A
		{ 0, 0, 0, 0, 0}, //B
		{ 3, 5, 1, 1, 1}, //C
		{ 4, 5, 1, 3, 0}, //D
		{ 5, 7, 1, 3, 0}, //E
		{ 9, 9, 1, 1, 0}, //F
		{ 2, 3, 1, 1, 0}, //G
		{ 6, 7, 3, 1, 0}, //H
		{ 6,11, 1, 1, 0}, //I
		{ 0, 0, 0, 0, 0}, //J
		{ 6,12, 2, 1, 0}, //K
		{ 6,11, 1, 1, 0}, //L
		{ 5, 9, 1, 1, 1}, //M
		{ 4, 6, 2, 2, 0}, //N
		{ 0, 0, 0, 0, 0}, //O
		{ 5, 7, 1, 1, 0}, //P
		{ 5, 8, 2, 2, 0}, //Q
		{ 6,12, 4, 1, 0}, //R
		{ 3, 5, 1, 2, 0}, //S
		{ 4, 7, 1, 2, 0}, //T
		{ 0, 0, 0, 0, 0}, //U
		{ 5, 9, 1, 1, 0}, //V
		{11,10, 2, 1, 0}, //W
		{ 0, 0, 0, 0, 0}, //X
		{ 9, 9, 1, 2, 0}, //Y
		{ 0, 0, 0, 0, 0}, //Z	
	};
	
	
	//C 4.9384 H 7.7583 N 1.3577 O 1.4773 S 0.0417 (Senko et al, J Am Soc MS 1995 pp. 229-233).
	private static double[] averagine = {4.9384, 7.7583, 1.3577, 1.4773, 0.0417}; 
	private static double averagineWeight = monoMass[C]*averagine[C] +
											monoMass[H]*averagine[H] +
											monoMass[N]*averagine[N] +
											monoMass[O]*averagine[O] +
											monoMass[S]*averagine[S];
	
	private static int[] getComponentElements(double avgMass){
		
		double factor= avgMass/averagineWeight;
		
		int[] elements = new int[NoElements];				
		elements[C] = (int)Math.round(factor*averagine[C]);
		elements[H] = (int)Math.round(factor*averagine[H]);
		elements[N] = (int)Math.round(factor*averagine[N]);
		elements[O] = (int)Math.round(factor*averagine[O]);
		elements[S] = (int)Math.round(factor*averagine[S]);
		
		double sumup =  monoMass[C]*elements[C] +
						monoMass[H]*elements[H] +
						monoMass[N]*elements[N] +
						monoMass[O]*elements[O] +
						monoMass[S]*elements[S];
		
		elements[H] += (int)Math.round((avgMass-sumup)/monoMass[H]);
		
		return elements;
	}
	
	
	public static ArrayList<Double> get(double avgMass){
		return get( getComponentElements( avgMass ) );
	}
	
	
	public static double[][] getAveragineIsotopeTable(double minMass, double maxMass, double interval){
		
		int tSize= (int)Math.round((maxMass-minMass)/interval) + 1;		
		double[][] table= new double[tSize][];
		
		for(int m=1; m<tSize; m++) {
		//	System.out.println(m*interval);
		//	int[] avgChemComp = getComponentElements( m*interval );		
		//	ArrayList<Double> isotopePattern = get(avgChemComp);
			ArrayList<Double> isotopePattern= get(m*interval);
			
			table[m] = new double[isotopePattern.size()];
			for(int i=0; i<isotopePattern.size(); i++)
				table[m][i]= isotopePattern.get(i);
			
		}

		return table;
		
	}
	
	
	
	public static double[] getPredictedDdeu(String peptide, String strDnat, int hdxnum){
		
		int exchangeableSites = peptide.length()-1;
		for(int i=0 ; i<peptide.length() ; i++)
			if( peptide.charAt(i)=='P' ) exchangeableSites--;
		
		ArrayList<Double> natDist= new ArrayList<Double>();
		String[] nTok= strDnat.split(";");
		for( String s: nTok ) natDist.add(Double.parseDouble(s));
		
		ArrayList<Double> deuLevelDist= new ArrayList<Double>();
		double exProb 	  = (double)hdxnum/exchangeableSites;	//exSiteNum: number of exchangeable sites (id sequence length -p site# -1)
		for(int j=0 ; j<=exchangeableSites ; j++) {		//p: probability for the binomial distribution
			deuLevelDist.add( Binomial.pdf(exchangeableSites, exProb, j) );						//binomial distribution
		}
		
		Convolution conv= new Convolution(natDist, deuLevelDist);
		ArrayList<Double> cvRes= conv.getConvResult();
		double[] ddeu= new double[cvRes.size()];
		for(int i=0; i<ddeu.length; i++) {
			ddeu[i]= cvRes.get(i);
		}
		
		return ddeu;
	}
	
	
	
}









































