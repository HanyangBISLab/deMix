package kr.ac.hanyang.bislab.demix.ms1;

import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;

import kr.ac.hanyang.bislab.demix.mzfile.MZFormatParser;
import kr.ac.hanyang.bislab.demix.Params;
import kr.ac.hanyang.bislab.demix.peptide.Peptide;
import kr.ac.hanyang.bislab.demix.peptide.PeptideFeature;

public class HDXExperiment {

	String label;
	String rawFile;
	String baseName;
	
//	private MSXMLParser parser;
	
	private MZFormatParser mzParser;
	private int maxScanNo;
	private TreeMap<Double, Integer> rt2scan;
	
	public HDXExperiment(String fileName, String l) throws Exception {
		
		label= l;
		rawFile= fileName;
		
		baseName= rawFile.replace("\\", "/");
		baseName= baseName.substring(baseName.lastIndexOf('/')+1, baseName.lastIndexOf('.'));	
		
	//	parser = new MSXMLParser(rawFile);	
	//	maxScanNo = parser.getMaxScanNumber();
		
		mzParser= MZFormatParser.get(rawFile);
		maxScanNo = mzParser.getMaxScanNumber();
		rt2scan = mzParser.getRt2ScanMap();
		
	}
	

	public PeptideFeature[] getAggreatedCluster(ArrayList<Peptide> peptList) throws Exception {
		
		int peptListSize= peptList.size();
		
		PeptideFeature[] aggreatedDeuDist= new PeptideFeature[peptListSize];	
		
		for(int pti=0; pti<peptListSize; pti++ ) {
			
			Peptide pept= peptList.get(pti);
			if( pept.getMS1Feature() == null ) continue;
			
			double targetedMZ= pept.getPMZ();
			double isotopeDiff= Params.DeutriumIsotopeSpace/pept.getCharge();
			double massTol= pept.getPeakWidth();
			
			PeptideFeature elutionProfile= pept.getMS1Feature();

			int maxClusterLength= pept.getIsotopeCluster().size()+pept.getHDXSiteNum()-1;//(int)(pept.getHDXSiteNum()*1.5);	
			int requiredObservedPeaks = (int)(pept.getMoleWeight()/1000)+3;

			//we have to convert scanNo. to rt time.
			Map.Entry<Double, Integer> alignScan= rt2scan.floorEntry(elutionProfile.getStartRT());
			int startScan = (alignScan != null)? alignScan.getValue() : rt2scan.firstEntry().getValue();
			
			alignScan= rt2scan.ceilingEntry(elutionProfile.getEndRT());
			int endScan = (alignScan != null)? alignScan.getValue() : rt2scan.lastEntry().getValue();

			alignScan= rt2scan.floorEntry(elutionProfile.getBestRT());
			int centerScan = (alignScan != null)? alignScan.getValue() : (startScan+endScan)/2;
		//	System.out.println(pept.getInputPeptide() + " " +pept.getCharge() + " " + startScan + " " +endScan);
			
		/*	int centerScan	= elutionProfile.getBestScan();
			int startScan	= elutionProfile.getStartScan();
			int endScan		= elutionProfile.getEndScan();//*/

			int initHDXNum= 0;
			int xstart 	= Math.max(startScan-5, 1);	//prevent startIndex scan <0
			for(int nh=centerScan; nh>=xstart; nh--){
				
				double[][] peakList= mzParser.getPeakList(nh);
				if( peakList == null ) continue;
				
				initHDXNum = getInitHdxNum(peakList, targetedMZ, isotopeDiff, massTol, requiredObservedPeaks, maxClusterLength);
				if( initHDXNum != 0 ) {
					centerScan = nh;
					break;
				}
			}
			if( initHDXNum == 0 ){
				int xend   	= Math.min(endScan+5, maxScanNo);//prevent e
				for(int nh=centerScan+1; nh<=xend; nh++){
					
					double[][] peakList= mzParser.getPeakList(nh);
					if( peakList == null ) continue;
					
					initHDXNum = getInitHdxNum(peakList, targetedMZ, isotopeDiff, massTol, requiredObservedPeaks, maxClusterLength);
					if( initHDXNum != 0 ) {
						centerScan = nh;
						break;
					}
				}
			}
			if( initHDXNum == 0 ) {
				continue;
			}//*/
			
		//	System.out.println(pept.getInputPeptide() + " $$ " +initHDXNum + " " + startScan + " " +endScan);

			startScan= Math.max(1, centerScan-Params.extendibleScanRange);	//prevent startIndex scan <0
			endScan  = Math.min(centerScan+Params.extendibleScanRange, maxScanNo);//prevent endIndex scan > end scan of control_mzxml
			
			int checkRange = (int)(pept.getMoleWeight()/1000 + 1);
			int possibleHDXStart = Math.max(initHDXNum-checkRange, 0); 
			int possibleHDXEnd = Math.min(initHDXNum+checkRange, maxClusterLength);

			double[] aggregatedDistribution= new double[maxClusterLength];		
			int[] aggregatedMembers= new int[maxClusterLength];
			
			int minAggedScan= centerScan, maxAggedScan= centerScan,  apexAggedScan= centerScan, conScanHole = 0;
		//	System.out.println(startScan + " " + centerScan + " " +endScan);
			
			double observedAmount= 0., maxAmount= 0.;
			for(int j=centerScan; j<endScan ; j++) {
			
				double[][] peakList= mzParser.getPeakList(j);
				if( peakList == null ) continue;

				observedAmount= aggreatedScan(aggregatedDistribution, aggregatedMembers, peakList, targetedMZ, isotopeDiff, massTol, 
						possibleHDXStart, initHDXNum, possibleHDXEnd, requiredObservedPeaks+1, maxClusterLength);
				
				if( observedAmount < 1 ){
					conScanHole++;
				}
				
				else {						
					conScanHole = 0;
					if( maxAggedScan < j ) maxAggedScan= j;
					if( maxAmount < observedAmount ) {
						maxAmount= observedAmount;
						apexAggedScan= j;
					}
				}
				
				if( conScanHole > 2 ) break;
			}
			
			observedAmount= 0;
			conScanHole = 0;
			for(int j=centerScan-1; j>=startScan ; j--) {
				
				double[][] peakList= mzParser.getPeakList(j);
				if( peakList == null ) continue;
				
				observedAmount= aggreatedScan(aggregatedDistribution, aggregatedMembers, peakList, targetedMZ, isotopeDiff, massTol, 
						possibleHDXStart, initHDXNum, possibleHDXEnd, requiredObservedPeaks+1, maxClusterLength);
				
				if(  observedAmount < 1 ) conScanHole++;
				else {
					conScanHole = 0;
					if( j < minAggedScan ) minAggedScan= j;
					if( maxAmount < observedAmount ) {
						maxAmount= observedAmount;
						apexAggedScan= j;
					}
				}
				if( conScanHole > 2 ) break;
			}
						
			PeptideFeature pf= new PeptideFeature(	minAggedScan,  mzParser.getRetentionTime(minAggedScan),
													maxAggedScan,  mzParser.getRetentionTime(maxAggedScan),
													apexAggedScan, mzParser.getRetentionTime(apexAggedScan),
													apexAggedScan, mzParser.getRetentionTime(apexAggedScan),
													0. );
			
		//	for(int i=0; i<aggregatedDistribution.length; i++) System.out.println(aggregatedDistribution[i]);
			pf.setIsotopeCluster(aggregatedDistribution);
			aggreatedDeuDist[pti]= pf;
			
		//	for(int i=0;i<maxClusterLength; i++)
		//		System.out.println(aggregatedMembers[i]);			
		}

		return aggreatedDeuDist;
	}
	
	
//	private int n_GetScanSum(int idindex, int scanNum, int n_startHDXNum, int n_initHDXNum, int n_endHDXNum, int prePeakNum) throws IOException{
	
	private double aggreatedScan(double[] distAggDeu, int[] aggregatedMembers, double[][] peakList, double targetedMZ, double isotopeDiff, double massTol, 
									int possibleHDXStart, int possibleHDXNum, int possibleHDXEnd, int requiredObservedPeaks, int maxClusterLength){
		
	//	double[][] peakList= scan.getMassIntensityList();//0:mass, 1:inten
		int peaksCnt= peakList[0].length;
		
		double[] observedDist = new double[maxClusterLength];
//		double[] observedError = new double[maxClusterLength];
		
		int k= 0;		
		for(int pi=0 ; pi<maxClusterLength ; pi++){		//deuterium exchanged peaks max#: maxPeakNum
			
			int tpindex= -1;
			double tpdelta= massTol+1.;
			
			double lowerBound = targetedMZ - massTol;
			double upperBound = targetedMZ + massTol;
			
		//	double matchSum = 0.;
			for( ; k<peaksCnt; k++) {
				
				if( peakList[0][k] < lowerBound ) continue;
				else if( peakList[0][k] > upperBound ) {
					break;
				}
			/*	if( scanNum == 286 ) {
					System.out.println(peakList[0][k] + "\t" + peakList[1][k] + "\t" + (targetedMZ-peakList[0][k]) + "\t" + massTol);
				}//*/
			//	observedDist[pi] += peakList[1][k];
			//	observedError[pi] += peakList[1][k]*peakList[0][k];
				
				double merr = Math.abs(targetedMZ-peakList[0][k]);
				if( merr < tpdelta ) {								
					tpdelta = merr;
					tpindex = k;
				}//*/
			}
			
			if( tpindex != -1 ) {
				observedDist[pi] = peakList[1][tpindex];
			//	observedError[pi] = peakList[0][tpindex]-targetedMZ;
				k= tpindex+1;
			}//*/
		//	if( observedError[pi] != 0. )
		//		observedError[pi] = observedError[pi]/observedDist[pi]-targetedMZ;
			targetedMZ += isotopeDiff;
		}
	/*	if( scanNum == 286 ) {
			for(int i=0;i<maxClusterLength; i++)
				System.out.println(observedDist[i]+"\t"+observedError[i]);
		}//*/
		double minHdxInt= Params.minPeaksIntensiy;
		
		int distPCnt = 0;					
		int local_series = 0, start_dist = 0, end_dist = 0;
		boolean n_criteria = true;
		
		for(int n=0; n<maxClusterLength; n++){			//check the number of valid peaks 
			
			if( observedDist[n] > minHdxInt ){
				local_series++;
			}
			else {
				if( distPCnt < local_series ){
					distPCnt = local_series;
					start_dist = n - local_series;
					end_dist = n - 1;
				}
				local_series = 0;
			}
			
			if( possibleHDXStart <= n && n<=possibleHDXEnd ){
				if( observedDist[n] < minHdxInt ) n_criteria = false;
			}
		
		}
		if( distPCnt < local_series ){
			distPCnt = local_series;
			start_dist = maxClusterLength - local_series;
			end_dist = maxClusterLength - 1;
		}
		
		if( !n_criteria ) return 0;		
		if( possibleHDXNum < start_dist || end_dist < possibleHDXNum ) return 0;	
				
	/*	for(int n=0; n<maxPeakNum; n++){
			if( n < start_dist || n > end_dist ) observedDist[n] = 0;
		}//*/
	
		double amount = 0.;
		if( requiredObservedPeaks < distPCnt ){//if number of valid peaks(peaks which has intensity>minInt) exceeds minPeakNum
		
			for(int n=0; n<maxClusterLength; n++){ 				//and apex intensity of deuteriun exchanged distribution> max intensity of previous peaks,
				if( start_dist <= n && n <= end_dist ) {
					distAggDeu[n] += observedDist[n];				//sum up the sub distributions to make a distribution sum
					amount += observedDist[n];
					aggregatedMembers[n]++;
				}
			}
			return amount;
		} 
		else return 0;

	}
	
	
	
	
	private int getInitHdxNum(double[][] peakList,  double targetedMZ, double isotopeDiff, double massTol, 
															int requiredObservedPeaks, int maxClusterLength){
		
	//	Scan scan = parser.rap(scanNum);	
	//	double[][] peakList= scan.getMassIntensityList();//0:mass, 1:inten
		int peaksCnt= peakList[0].length;

		double[] observedDist = new double[maxClusterLength];
		
		int k= 0;		
		for(int pi=0 ; pi<maxClusterLength ; pi++){		//deuterium exchanged peaks max#: maxPeakNum
			
			double lowerBound = targetedMZ - massTol;
			double upperBound = targetedMZ + massTol;
			
			for( ; k<peaksCnt; k++) {
				
				if( peakList[0][k] < lowerBound ) continue;
				else if( peakList[0][k] > upperBound ) {
					break;
				}
				observedDist[pi] += peakList[1][k];
			}
			
			targetedMZ += isotopeDiff;
		}
	//	System.out.println(scanNum + "  " +maxClusterLength + " // " );
	//	for(int i=0;i<maxClusterLength;i++)
	//		System.out.println(observedDist[i]);

		int local_series = 0;	
		int maxindex = 0, serindex = 0;
		double localmax = 0, series_max = 0;
		
		for(int n=0; n<maxClusterLength; n++){//check the number of valid peaks 
			
			if( 0 < observedDist[n] ){
				local_series++;
				if( series_max < observedDist[n] ){
					series_max = observedDist[n];
					serindex = n;
				}
			}
			else {
				if( requiredObservedPeaks < local_series ){
					if( localmax < series_max ){
						localmax = series_max;
						maxindex = serindex;
					}
				}
				local_series = 0;
				series_max = 0;
			}		
		}
		
		if( requiredObservedPeaks < local_series ){
			if( localmax < series_max ){
				localmax = series_max;
				maxindex = serindex;
			}
		}
		
		return maxindex;
	}
	
	
	
	
	
	
	
	
}















