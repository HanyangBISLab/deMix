package kr.ac.hanyang.bislab.demix;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;

import kr.ac.hanyang.bislab.demix.distribution.BimodalFitting;
import kr.ac.hanyang.bislab.demix.distribution.Binomial;
import kr.ac.hanyang.bislab.demix.distribution.Convolution;
import kr.ac.hanyang.bislab.demix.distribution.DDMatch;
import kr.ac.hanyang.bislab.demix.distribution.UnimodalFitting;
import kr.ac.hanyang.bislab.demix.ms1.Experiment;
import kr.ac.hanyang.bislab.demix.ms1.HDXExperiment;
import kr.ac.hanyang.bislab.demix.peptide.Peptide;
import kr.ac.hanyang.bislab.demix.peptide.PeptideFeature;
import kr.ac.hanyang.bislab.demix.protein.ProxDB;

public class HDXAnalysis {

	static double grand_ms1_tolerance = 0.1;
	
	ArrayList<Peptide> peptList;
	
	String[] 			labeledTag;
	PeptideFeature[][] 	aggreatedDeuDist;//condition, peptidenx
	DDMatch[][] 		analysisResults;
	boolean 			hasProtein;
	
	String str_peptide 	="peptide";
	String str_charge	="charge";
	String str_mz		="mz";
	
	int col_peptide;
	int col_charge;
	int col_mz;	 
		
	public HDXAnalysis(deMix hdx, int groupIndex) throws Exception {
				
		peptList= new ArrayList<Peptide>();

		HashSet<String> seenPept= new HashSet<String>();
		
		String buf;
		BufferedReader in = new BufferedReader( new FileReader(hdx.pept_input_file) );	
		setTSVColumnIndex(in.readLine(), hdx.pept_input_file);
		
		while( (buf = in.readLine()) != null ){
			
			String[] tok= buf.split("\t");
			if( tok.length < 3 ) continue;
			
			int charge= Integer.parseInt(tok[col_charge]);
			String stripSeq= Peptide.getPeptideSequence(tok[col_peptide]);
			double pmz= Peptide.getPeptMZ(stripSeq, charge);
			
			if( !seenPept.add(tok[col_peptide]+"_"+charge+"_"+(int)pmz) ) continue;

			double peakW= (hdx.PPMTolerance == 0)? hdx.DATolerance : (pmz/1000000*hdx.PPMTolerance);
			double wideW= 0;//pmz/10000.;//i.e., 100ppm
			
			peptList.add( new Peptide(peptList.size(), tok[col_peptide], charge, pmz, Math.max(peakW, wideW)) );
		}
		in.close();
		System.out.println("#Input Peptides:: " + peptList.size());
				
		//Get RT time of peptides from ctrl experiment **********************************************************************************************
		String ctrl_ms_file= hdx.groupList.get(groupIndex).ctrl_ms_file;
		Experiment ctrlE= new Experiment(ctrl_ms_file, "CTRL");
		ctrlE.getElutionProfile(peptList);
		
	/*	//Correcting systematic error bias***********************************************************************************************************
		ArrayList<Double> errorAtGlance= new ArrayList<Double>();
		for( Peptide pept : peptList ) {
			if( pept.getMS1Feature() != null ) {
				errorAtGlance.add( (pept.getMS1Feature().getWeightedAverMZ()-pept.getPMZ()) * 1000000 / pept.getPMZ() );
			}
		}
		
		double ppmshift= 0;
		if( errorAtGlance.size() > 2 ) {
			Collections.sort(errorAtGlance);
			ppmshift= errorAtGlance.get(errorAtGlance.size()/2);
			System.out.println( String.format("Systematic m/z correction:: %.2f ppm", ppmshift) );
		}
		for( Peptide pept : peptList ) {
			pept.correctSystemBias(0, (hdx.PPMTolerance == 0)? hdx.DATolerance : (pept.getPMZ()/1000000*hdx.PPMTolerance));
		}
		ctrlE.getElutionProfile(peptList);//*/

		
		//align to protein sequence *****************************************************************************************************************
		if( hdx.fasta_input_file != null ) {
			hasProtein= true;
			System.out.print("#Input Proteins:: ");
			
			ProxDB pdb= new ProxDB();
			pdb.readFasta(hdx.fasta_input_file);
			
			for( Peptide pept : peptList ) {
				pept.setProtMatch( pdb.match2Protein(pept.getStripSequence()) );
			}
			
			Collections.sort(peptList, new Comparator<Peptide>() {
				public int compare(Peptide aa, Peptide bb) {
					if( aa.getMatch2Protein().compareTo(bb.getMatch2Protein()) < 0 ) return -1;
					else if( aa.getMatch2Protein().compareTo(bb.getMatch2Protein()) > 0 ) return 1;
					return 0;
				}
			});
			
			for(int i=0;i<peptList.size(); i++) {
				peptList.get(i).setIndex(i);
			}
		}
		else {
			hasProtein= false;
			
			Collections.sort(peptList, new Comparator<Peptide>() {
				public int compare(Peptide aa, Peptide bb) {
					if( aa.getIndex() < bb.getIndex() ) return -1;
					else if( aa.getIndex() > bb.getIndex() ) return 1;
					return 0;
				}
			});
		}
		System.out.println();

		
		//Get deuterated peptides from HDX experiment ***********************************************************************************************
	
		LinkedHashMap<String, String> hdx_ms_files= hdx.groupList.get(groupIndex).hdx_ms_files;
		
		labeledTag= new String[hdx_ms_files.size()];
		aggreatedDeuDist= new PeptideFeature[hdx_ms_files.size()][];		
		
		int label= 0;
		for( Map.Entry<String, String> entry : hdx_ms_files.entrySet() ) {
			System.out.println("Analysing : " + entry.getKey() +", "+entry.getValue());
			HDXExperiment hdxexp= new HDXExperiment(entry.getValue(), entry.getKey());		
			labeledTag[label]= entry.getKey();
			aggreatedDeuDist[label++]= hdxexp.getAggreatedCluster(peptList);
		}
		
		analysisResults= new DDMatch[labeledTag.length][peptList.size()];
	}
	
	
	
	
	
	public void decode() {

		int peptListSize= peptList.size();

		ArrayList<Double>[][] calcDeuDist= new ArrayList[peptListSize][];
		int[][] calc_max_index= new int[peptListSize][];
				
		for(int pti=0; pti<peptListSize; pti++ ) {
			
			Peptide pept= peptList.get(pti);
			
			int exchangeableSites= pept.getHDXSiteNum();
			calcDeuDist[pti]= new ArrayList[exchangeableSites+1];
			calc_max_index[pti]= new int[exchangeableSites+1];
			
			ArrayList<Double> natDist= pept.getIsotopeCluster();			
			for(int hlevel=0 ; hlevel<=exchangeableSites; hlevel++){
				
				//preparation////////////////////////////////////////////////////////////////////////////////////////////////////
				ArrayList<Double> deuLevelDist= new ArrayList<Double>();
				double exProb 	  = (double)hlevel/exchangeableSites;	//exSiteNum: number of exchangeable sites (id sequence length -p site# -1)
				for(int j=0 ; j<=exchangeableSites ; j++) {		//p: probability for the binomial distribution
					deuLevelDist.add( Binomial.pdf(exchangeableSites, exProb, j) );						//binomial distribution
				}
				
				Convolution conv= new Convolution(natDist, deuLevelDist);
				calcDeuDist[pti][hlevel]= conv.getConvResult();
				calc_max_index[pti][hlevel]= conv.getMaxIndex();//////////////////////////////////////////////////////////////////////////////////////////////////				
			}
		}
		
		for(int pti=0; pti<peptListSize; pti++ ) {
			
			for(int label=0; label<labeledTag.length; label++) {
				
				if( aggreatedDeuDist[label][pti] == null ) continue;
				
				ArrayList<Double> aggDeuDist= aggreatedDeuDist[label][pti].getIsotopeCluster();
				
				//score possible single hdx
				ArrayList<DDMatch> unimodalAnalysis= UnimodalFitting.run(aggDeuDist, calcDeuDist[pti], calc_max_index[pti]);
				if( unimodalAnalysis.size() == 0 ) continue;
				
				if( unimodalAnalysis.get(0).isExplainedEnough(0.9) ) {
					analysisResults[label][pti]= unimodalAnalysis.get(0);
					continue;
				}
				
				//bimodal hdx
				DDMatch bimodalAnalysis = BimodalFitting.run(unimodalAnalysis, aggDeuDist, calcDeuDist[pti]);
				analysisResults[label][pti]= bimodalAnalysis;
			}
		}
	}
	
	
	public void printHDXProfile(String outfile) throws Exception {
		PrintWriter out = new PrintWriter( new BufferedWriter( new FileWriter(outfile) ) );
		printHDXProfile(out);
		out.close();
		System.out.println("The results were written to " + outfile);
	}
	
	public void printHDXProfile(PrintWriter out) throws Exception {
		
		int peptListSize= peptList.size();
		
		out.print("Id\tm/z\tCharge\tPeptide\t");		
		if( hasProtein ) out.print("Protein\tPosFrom\tPosTo\t");
		out.print("ExpM/z\tm/zShift\tStartScan\tEndScan\tApexScan\tApexRT(m)\tPredictedDnat");
		
		for(int label=0; label<labeledTag.length; label++) {
			out.printf("\t%s (#HDX)", labeledTag[label]);
		}
		for(int label=0; label<labeledTag.length; label++) {
			out.printf("\t%s (%cHDX)", labeledTag[label], '%');
		}
		out.println();
		
		
		for(int pti=0; pti<peptListSize; pti++ ) {
			
			Peptide pept= peptList.get(pti);
			out.printf("%d\t%.3f\t%d\t%s\t", pept.getIndex()+1, pept.getPMZ(), pept.getCharge(), pept.getInputPeptide());
			
			if( hasProtein ) out.printf("%s\t", pept.getMatch2Protein().toString()); 
			
			if( pept.getMS1Feature() == null ) out.print("-\t-\t-\t-\t-\t-");
			else out.print(pept.getMS1Feature().toString(pept.getPMZ()));
			
			out.print("\t"+pept.toStringOfIsotopeCluster());
			
			for(int label=0; label<labeledTag.length; label++) {
				if( analysisResults[label][pti] == null ) out.printf("\t-");
				else out.printf("\t%s", analysisResults[label][pti].getStrHDXNumer() );
			}
			
			for(int label=0; label<labeledTag.length; label++) {
				if( analysisResults[label][pti] == null ) out.printf("\t-");
				else out.printf("\t%s", analysisResults[label][pti].getStrHDXRate(pept.getHDXSiteNum()) );
			}
			
			out.println();
		}
	}	
	
	public void printDdeuAnalysis(String outfile) throws Exception {
		PrintWriter out = new PrintWriter( new BufferedWriter( new FileWriter(outfile) ) );
		printDdeuAnalysis(out);
		out.close();
		System.out.println("The results were written to " + outfile);
	}
	public void printDdeuAnalysis(PrintWriter out) throws Exception {
		
		out.printf("Id\tm/z\tCharge\tPeptide\tD2OLabel\tDeuRate(%c)\t1stDeuNum\t1stDeuWeight\t2ndDeuNum\t2ndDeuWeight\tPredictedDdeu\t", '%');
		out.println("StartScan\tEndScan\tApexScan\tApexRT(m)\tObservedDdeu(aggregated)\tMatchedScore");

		int peptListSize= peptList.size();
		for(int pti=0; pti<peptListSize; pti++ ) {
			for(int label=0; label<labeledTag.length; label++) {
				
				Peptide pept= peptList.get(pti);
				out.printf("%d\t%.3f\t%d\t%s\t%s\t", pept.getIndex()+1, pept.getPMZ(), pept.getCharge(), pept.getInputPeptide(), labeledTag[label]);		
				if( analysisResults[label][pti] == null ) {
					for(int i=0; i<12; i++)
						out.printf("-\t");
				}
				else {
					out.printf("%s\t%s\t%s", analysisResults[label][pti].toOutput(pept.getHDXSiteNum()), aggreatedDeuDist[label][pti].toOutput(),
						analysisResults[label][pti].getStrMatchQuality());
				}
				out.println();
			}
		}
	}
	
	
	
	
	
	
	
	
	private void setTSVColumnIndex(String header, String file) throws Exception {
		
		String[] tok= header.split("\t");
		
		int required= 0;
		for(int i=0; i<tok.length; i++){
			if( tok[i].compareToIgnoreCase(str_peptide) == 0 ){
				col_peptide= i;
				required++;
			}
			else if( tok[i].compareToIgnoreCase(str_charge) == 0 ){
				col_charge= i;
				required++;
			}
			else if( tok[i].compareToIgnoreCase(str_mz) == 0 ){
				col_mz= i;
				required++;
			}
		}
		
		if( required < 3 ) {
			System.out.println(header);
			System.out.println("[deMix] Unmatched column names! Required columns: peptide, charge, mz");
			System.exit(1);
		}	
		
		System.out.println("Matched input columns in " + file);
		System.out.println("Peptide is in the " + (col_peptide+1) + "-th cloumn." );
		System.out.println("Charge  is in the " + (col_charge+1) + "-th cloumn." );
		System.out.println("m/z     is in the " + (col_mz+1) + "-th cloumn." );
		
	}

	
	
	
	
	public void XXX_print(String outfile) throws Exception {
		
		PrintWriter out = new PrintWriter( new BufferedWriter( new FileWriter(outfile) ) );
		out.print("index\tpmz\tcharge\tpeptide\tstartScan\tendScan\tapexScan\tapexRT(m)");
		
		for(int label=0; label<labeledTag.length; label++) {
			out.printf("\t%s_hdx1\t%s_hdx2\t%s_explained\t%s_MPC", labeledTag[label], labeledTag[label], labeledTag[label], labeledTag[label] );
		}
		out.println();
		
		int peptListSize= peptList.size();
		for(int pti=0; pti<peptListSize; pti++ ) {
			
			Peptide pept= peptList.get(pti);
			out.printf("%d\t%.3f\t%d\t%s\t%s", pept.getIndex()+1, pept.getPMZ(), pept.getCharge(), pept.getInputPeptide(), pept.getMS1Feature().toString());
			
			for(int label=0; label<labeledTag.length; label++) {
				if( analysisResults[label][pti] == null ) out.printf("\t-\t-\t-\t-");
				else out.printf("\t%s", analysisResults[label][pti].toString() );
			}
			out.println();
		}
		
		out.close();
		
		System.out.println("The results were written to " + outfile);
	}
	
	private double calcStdDeviation(double aver, ArrayList<Double> mzError){
		double sd= 0;
		for( Double err : mzError)
			sd +=  Math.pow(aver-err, 2);
		sd = sd / (mzError.size()-1);
		return Math.sqrt(sd);	
	}
	
}







