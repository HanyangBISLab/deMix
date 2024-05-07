package kr.ac.hanyang.bislab.demix.ms1;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Set;

import kr.ac.hanyang.bislab.demix.mzfile.MZFormatParser;
import kr.ac.hanyang.bislab.demix.Params;
import kr.ac.hanyang.bislab.demix.peptide.Peptide;
import kr.ac.hanyang.bislab.demix.peptide.PeptideFeature;

public class Experiment {
	
	static double min_PIF_cluster = 0.5;
	static double min_similarity_cluster = 0.95;
	
	static double min_weight_similarity = 0.5;
	static double max_minus_oneDa_ratio = 0.5;
	
	String label;
	String rawFile;
	String baseName;
	
	MZFormatParser mzParser;
		
	public Experiment(String fileName, String l) throws Exception {
		
		label= l;
		rawFile= fileName;
		
		baseName= rawFile.replace("\\", "/");
		baseName= baseName.substring(baseName.lastIndexOf('/')+1, baseName.lastIndexOf('.'));	
		
		mzParser= MZFormatParser.get(rawFile);
	}
	
	static double neighboring_peak_max_ratio= 0.2;
	
	public void getElutionProfile(ArrayList<Peptide> peptList) throws Exception {

		//make sure sorted pmz 
		Collections.sort(peptList);
		
		int peptListSize= peptList.size();
		ArrayList<ObservedCluster>[] elutionChromatogram= new ArrayList[peptListSize];
		for(int i=0; i<peptListSize; i++) {
			elutionChromatogram[i]= new ArrayList<ObservedCluster>();
		}
	
		Set<Integer> scanSet= mzParser.getScanSet();
		for( Integer sn : scanSet ) {
			double[][] peakList = mzParser.getPeakList(sn);//0:mass, 1:inten
			int peaksCnt  = peakList[0].length;
			
			int start= 0;
			ArrayList<Integer>[] pmzIndex= new ArrayList[peptListSize];
			for(int pti=0; pti<peptListSize; pti++ ) {
				
				pmzIndex[pti]= new ArrayList<Integer>();
				
				double massTol= 	peptList.get(pti).getPeakWidth();
				double lowerBound = peptList.get(pti).getPMZ() - massTol;
				double upperBound = peptList.get(pti).getPMZ() + massTol;
				
				for(int k=start; k<peaksCnt; k++) {
					
					if( peakList[0][k] < lowerBound ) continue;
					else if( peakList[0][k] > upperBound ) {
						start= k;
						break;
					}
					if( peakList[1][k] > 0 ) pmzIndex[pti].add(k);
				}
				
				if( pmzIndex[pti].size() != 0 ) {
					start= pmzIndex[pti].get(0);
				}
			}//mono indices found
			
			
			for(int pti=0; pti<peptListSize; pti++ ) {
				elutionChromatogram[pti].add( getCluster(sn, peakList, pmzIndex[pti], peptList.get(pti).getCharge(), 
													peptList.get(pti).getTrimmedIsotopeCluster(), peptList.get(pti).getPeakWidth()) );
			}
		}//for each scan, MS1 spectrum
		
		for(int pti=0; pti<peptListSize; pti++ ) {
			
			ArrayList<ObservedCluster> choromatogram= elutionChromatogram[pti];
	
			int elutionApex= getScoreApex(choromatogram);
			if( elutionApex < 1 ) {//indices in terms of MS1, not scan number
				elutionApex = 1;
				continue;
			}
			
					
			int[] rtRange = getElutionRange(choromatogram, elutionApex);//indices in terms of MS1, not scan number
			int sttNat 		= rtRange[0];
			int endNat 		= rtRange[1];
			int maxNat		= elutionApex;
						
			double weightedAverScan = 0., weightedAverMZ = 0., amountSum = 0., maxWeight= 0.;
			for(int on=sttNat; on<=endNat; on++){
				
				ObservedCluster curClst= choromatogram.get(on);
				if( curClst == null ) continue;
				
				weightedAverScan = (weightedAverScan*amountSum + curClst.scanNo*curClst.amount) / (amountSum+curClst.amount);	
				weightedAverMZ   = (weightedAverMZ*amountSum + curClst.getWeightedMZ()) / (amountSum+curClst.amount);				
				amountSum += curClst.amount;
				if( maxWeight < curClst.amount ) {
					maxWeight= curClst.amount;
					maxNat= on;
				}
			}
			int elutionWeightedAver = (int)Math.round(weightedAverScan);
			if( elutionWeightedAver == 0 ) elutionWeightedAver= 1;
			
			int startScanNo= choromatogram.get(sttNat).scanNo;
			int endScanNo= choromatogram.get(endNat).scanNo;
			int apexScanNo= choromatogram.get(maxNat).scanNo;
			
			peptList.get(pti).setMS1Feature( new PeptideFeature(startScanNo, 			mzParser.getRetentionTime(startScanNo),
																endScanNo, 				mzParser.getRetentionTime(endScanNo),
																elutionWeightedAver, 	mzParser.getRetentionTime(elutionWeightedAver),
																apexScanNo, 			mzParser.getRetentionTime(apexScanNo),
																weightedAverMZ) );
		}
	}
	
	static int continuity= 3;
	private int getScoreApex(ArrayList<ObservedCluster> choromatogram){		
				
		int apex = 0;
		double maxScore = 0, maxAmount = 0;
		
		for(int i=continuity-1; i<choromatogram.size(); i++){			
			// at least 3 observed
			int observed= 0;
			for(int c=0; c<continuity; c++) {
				if( choromatogram.get(i-c) != null ) {
					observed++;
				}
			}
			if( observed != continuity ) continue;
			
			double quality = 0., amount= 0.;
			for(int c=0; c<continuity; c++) {
				quality += choromatogram.get(i-c).getWeightedScore();
				amount += choromatogram.get(i-c).amount;
			}
			
			quality /= continuity;
			if( quality < min_weight_similarity ) continue;
			
			amount /= continuity;
			
			if( quality > 0.90 ){
				if( maxAmount < amount ) {
					maxScore = quality;
					maxAmount = amount;
					apex = i;
				}
			}
			else {
				if( maxScore < quality ) {
					maxScore = quality;
					maxAmount = amount;
					apex = i;
				}
			}
		}	
		
		return apex-1;
	}
	
	
	private int[] getElutionRange( ArrayList<ObservedCluster> choromatogram, int apexIndex ){
		
		int sttNat= apexIndex, endNat= apexIndex, obConfirmed= 1;//Observed by apex;
		
		double minAmount= choromatogram.get(apexIndex).amount*0.1;
		int hole = 0;
		
		for(int on=apexIndex-1; on>=1; on--){
			if( choromatogram.get(on) != null && 
					minAmount < choromatogram.get(on).amount ) {
				hole = 0;
				obConfirmed++;
			}
			else hole++;
			
			sttNat = on;
			if( hole == Params.maxElutionHoleScans ){
				sttNat = on + hole;
				break;
			}
		}
		
		hole = 0;
		for(int on=apexIndex+1; on<choromatogram.size(); on++){
			if( choromatogram.get(on) != null && 
					minAmount < choromatogram.get(on).amount ) {
				hole = 0;
				obConfirmed++;
			}
			else hole++;
			
			endNat = on;
			if( hole == Params.maxElutionHoleScans ){
				endNat = on-hole;
				break;
			}
		}
		
		int[] range = new int[3];
		range[0] = sttNat;
		range[1] = endNat;
		range[2] = obConfirmed;
		return range;	
	}
		
	private ObservedCluster getCluster(int scanNo, double[][] ms1spectrum, ArrayList<Integer> monoIndex, 
											int charge, ArrayList<Double> distNat, double tolerance) {
		
		double isotopeDiff= Params.IsotopeSpace/charge;

		int peaksCnt= ms1spectrum[0].length;
		double mono_mass= -1., lsf_score= 0., lsf_sum= 0., lsf_pif= 1.;										
		for(int mix=0; mix<monoIndex.size(); mix++){	

			int mono_index 	 = monoIndex.get(mix);						
			double isotopeMZ = ms1spectrum[0][mono_index];
			int search_index = mono_index+1;
			
			ArrayList<Double> distObserved = new ArrayList<Double>();//observed nat dist
			distObserved.add(ms1spectrum[1][mono_index]);

			double obNatSum= ms1spectrum[1][mono_index];
			
			for(int k=1; k<distNat.size(); k++){
				
				isotopeMZ += isotopeDiff;
				
				double tpdelta =tolerance, tpInten= distObserved.get(distObserved.size()-1)*neighboring_peak_max_ratio; // matched isotope peak's error
				int tpindex = -1; // matched peak index
				
				double lowerBound = isotopeMZ - tpdelta;
				double upperBound = isotopeMZ + tpdelta;
				
				for(int m=search_index ; m<peaksCnt; m++){			
					
					if( ms1spectrum[0][m] < lowerBound ) continue;
					else if( ms1spectrum[0][m] > upperBound ) break;
					
					double merr = Math.abs(isotopeMZ-ms1spectrum[0][m]);
					if( merr < tpdelta && ms1spectrum[1][m] > tpInten ) {								
						tpdelta = merr;
						tpindex = m;
					}
				}

				if( tpindex < 0 ) break;
				else {
					obNatSum += ms1spectrum[1][tpindex];
					distObserved.add( ms1spectrum[1][tpindex] );
					
					search_index = tpindex+1;
					isotopeMZ = ms1spectrum[0][tpindex];	
				}
			}
			if( distObserved.size() < distNat.size() ) continue;
			
			double minus1DaMz 		= ms1spectrum[0][mono_index]- isotopeDiff;
			double minus1DaHeight 	= 0.;
			double regionalSum= 0.;

			double lowerBound = ms1spectrum[0][mono_index]- isotopeDiff- tolerance; // previous.
			for(int m=mono_index-1; m>-1; m--){
				if( ms1spectrum[0][m] < lowerBound ) break;
				regionalSum += ms1spectrum[1][m];
				if( Math.abs(ms1spectrum[0][m]-minus1DaMz) <= tolerance/2. ) {
					if( minus1DaHeight < ms1spectrum[1][m] ) minus1DaHeight= ms1spectrum[1][m];
				}
			}
			if( ms1spectrum[1][mono_index]*max_minus_oneDa_ratio < minus1DaHeight ) continue;
			
			double upperBound = ms1spectrum[0][search_index-1]+ isotopeDiff+ tolerance;
			for(int m=mono_index; m<peaksCnt; m++){
				if( upperBound < ms1spectrum[0][m] ) break;
				regionalSum += ms1spectrum[1][m];
			}
			if( obNatSum/regionalSum < min_PIF_cluster ) continue; 
		
			double similarityScore= dotProduct(distNat, distObserved);			
			if( similarityScore < min_similarity_cluster ) continue;
			
			if( lsf_score < similarityScore ){
				lsf_score = similarityScore;
				lsf_sum = obNatSum;
				lsf_pif = obNatSum/regionalSum;
				mono_mass = ms1spectrum[0][mono_index];
			}
		}
		return new ObservedCluster(scanNo, mono_mass, lsf_sum, lsf_score, lsf_pif);
	}
	
	private double dotProduct(ArrayList<Double> one, ArrayList<Double> two){
		double dotP = 0;
		for(int i=0; i<two.size(); i++){
			dotP += one.get(i)*two.get(i);
		}
		double oneLen = 0;
		for(int i=0; i<one.size(); i++){
			oneLen += one.get(i)*one.get(i);
		}
		double twoLen = 0;
		for(int i=0; i<two.size(); i++){
			twoLen += two.get(i)*two.get(i);
		}
		return dotP/(Math.sqrt(oneLen)*Math.sqrt(twoLen));
	}

	
	
}












