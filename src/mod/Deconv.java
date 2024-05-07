package mod;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import dst.IsoCluster;
import ppk.Peak;
import ppk.PeakData;
import ppk.PeakProcessing;
import rds.RawData;
import rds.ReaderFactory;
import utils.Constants;
import utils.Constants.PEAK_FIT_TYPE;

public class Deconv {
	public DeconvPept dp;
	public PeakProcessing peakPs;
	public Deconv() {
		dp = new DeconvPept();
		peakPs = new PeakProcessing();
	}
	
	public static long stime;
	
	
	public static void main(String[] args) {
		Deconv dev = new Deconv();
		int argc = args.length;
		int start_idx = -1, end_idx = -1;
		if (argc == 3 ) {
			start_idx = Integer.parseInt(args[1]);
			end_idx = Integer.parseInt(args[2]);
		}
		else if( !dev.isUsageRight(args) ){
			return;
		}
		
		Param param = new Param();
		param.readFile(args[0]);
		
		param.dump();	
		
		// load param
		String DataType;
		String PrintPk;
		DataType = param.get("File::DataType");
		if (DataType != null) {
			if (DataType.compareTo("MZXML") == 0) Constants.DataFileType = 6;
		}

		//System.out.println("%d",DataType);

		Constants.maxCharge			= Integer.parseInt(param.get("DeconvPep::MaxCharge"));
		Constants.thScore 			= Double.parseDouble(param.get("DeconvPep::ThScore"));

		Constants.maxAbundancePeak	= Integer.parseInt(param.get("AdvDeconv::MaxAbundancePeak"));
		Constants.scanNoModifier	= Integer.parseInt(param.get("AdvDeconv::ScanNoModifier"));		// 1.02
		Constants.maxMissPeak		= Integer.parseInt(param.get("AdvDeconv::MaxMissPeak"));		// 1.02
		Constants.massErr			= Double.parseDouble(param.get("AdvDeconv::MassErr"));
		Constants.thClustExt		= Double.parseDouble(param.get("AdvDeconv::ThClustExt"));
		Constants.Partition			= Integer.parseInt(param.get("PeakPicking::Partition"));		// 1.09
		
		double snt = 0.0, bgr = 0.0, bgr2 = 0.0;
		
		snt		= Double.parseDouble(param.get("PeakPicking::SNRThreshold"));
		bgr		= Double.parseDouble(param.get("PeakPicking::BackgroundRatio"));
		bgr2	= Double.parseDouble(param.get("PeakPicking::BackgroundRatio2"));
		PrintPk = param.get("PeakPicking::PrintPeak");
		if (PrintPk != null) {
			if (PrintPk.compareTo("YES") == 0) Constants.PrintPeak = true;
		}
		Constants.ResultOrder = param.get("DeconvPep::ResultOrder");	// 1.01
		Constants.OutputFormat = param.get("DeconvPep::OutputFormat");	// 1.02
		Constants.MSMS = param.get("DeconvPep::MakeMSMS");				// 1.1
		if (Constants.MSMS != null) {
			if (Constants.MSMS.compareTo("MASSONLY") == 0) Constants.MakeMSMS = 1;
			else if (Constants.MSMS.compareTo("DEISOTOPING") == 0) Constants.MakeMSMS = 2;
		}


		/*	Target = CParam::get("DeconvPep::Target");				// 1.03
		if (Target == NULL) {
		Target = new char[3];
		strcpy(Target, "MS") ;
		}*/
		/*	Truncated = CParam::get("DeconvPep::Truncated");			// 1.07
		if (Truncated == NULL) {
		Truncated = new char[3];
		strcpy(Truncated, "NO") ;
		}*/
		if (Constants.ResultOrder == null) {
			Constants.ResultOrder = new String();
			Constants.ResultOrder = "ABUNDANCE";
		}
		if (Constants.OutputFormat == null) {
			Constants.OutputFormat = new String();
			Constants.OutputFormat = "PEK";
		}

		if (Constants.maxAbundancePeak < 1) {
			System.out.println("ERROR: MaxAbundancePeak must be higher than 1\n");
			return ;
		}
		if (snt < 0.0) {
			System.out.println("ERROR: SNRThreshold must be positive\n");
			return ;
		}
		if (bgr < 0.0) {
			System.out.println("ERROR: BackgroundRatio must be positive\n");
			return ;
		}
		if (Constants.Partition <= 0) {
			System.out.println("ERROR: Partition must be positive\n");
			return ;
		}
		if (Constants.maxCharge <= 0) {
			System.out.println("ERROR: MaxCharge must be positive\n");
			return ;
		}
		if (Constants.maxMissPeak <= 0) {
			System.out.println("ERROR: MaxMissPeak must be positive\n");
			return ;
		}
		if (Constants.massErr <= 0.0) {
			System.out.println("ERROR: MassErr must be positive\n");
			return ;
		}
		/*	if (CIntsRangeFunc::FctRange > 0.9 || CIntsRangeFunc::FctRange < 0.0) {
		System.out.println("ERROR: IntsRangeErr must be 0.0 to 0.9\n");
		return -1;
		}*/
		// file input/output

		List<String> vif = param.getList("<File::DataList>");
		System.out.println("Read " + vif.size() + " files...");
		if (vif.size() > 0) {
			for (int i = 0; i < vif.size(); ++i) {
				System.out.println("Processing "+vif.get(i)+"...\n");
				try {
					dev.process_rawfile(vif.get(i), Constants.OutputFormat, param, start_idx, end_idx);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		else {
			System.out.println("ERROR: Num of DataFiles != Num of ResultFiles\n");
			return ;
		}
		
		System.out.println("[jRAPID] BYE!!");
	}
	
	double get_TIC(double min_mz, double max_mz, List<Peak> vec_pk, float minIntensity,
			double[] bpi, double[] bp_mz)
		{
			int numPts = vec_pk.size();
			if (numPts == 0)
				return 0;

			double sum = 0;
			//int numPtsUsed = 0;
			for (int i = 0; i < numPts; i++)
			{
				if (vec_pk.get(i).mz() >= min_mz && vec_pk.get(i).mz() < max_mz && vec_pk.get(i).in() >= minIntensity)
				{
//					System.out.println(vec_pk.get(i).in() + "\t" + vec_pk.get(i).mz());
					sum += vec_pk.get(i).in();
					if (vec_pk.get(i).in() > bpi[0])
					{
						bpi[0] = vec_pk.get(i).in();
						bp_mz[0] = vec_pk.get(i).mz();
					}
				}
			}
			return sum;
		}
	
	public List<IsoCluster> process(List<Peak> pkl, List<Peak> pks) {
		return dp.process(pkl, pks);
	}
	
	/**
	 *  sets the options for this instance.
	 * @param s2n sets the threshold signal to noise value.
	 * @param thresh sets the peak intensity threshold.
	 * @param thresholded this is true if FileType == FINNIGAN or MzXMLRAWDATA
	 * @param type sets the type of peak fitting algorithm used.
	 */
	public void setOptions(double s2n, double thresh, boolean thresholded, PEAK_FIT_TYPE type){
		peakPs.setOptions(s2n, thresh, thresholded, type);
	}
	
//	public DeconvolutedPair determineMono(double monomz, List<IsoCluster> plist) {
//		List<Peak> res_peak = new ArrayList<>();
//		
//		int cs = 0;
//		double m = 0, in = 0.0;
//		
//		//determine precursor m/z
//		for (int j = 0; j < plist.size(); ++j) {
//			int check_flag = 0;
//			for (int k = 0; k<plist.get(j).plist.size(); k++) {
////				monomz = plist.get(j).plist.get(k).mdbl_mz;
//				if (plist.get(j).plist.get(k).mz() > monomz*(1 - Constants.massErr/*1.0E-05*/) && plist.get(j).plist.get(k).mz() < monomz*(1 + Constants.massErr) && in < plist.get(j).sum) {
//					if (in == 0) {
//						in = plist.get(j).plist.get(plist.get(j).abundant).in();
//						m = plist.get(j).mass;
//						cs = plist.get(j).charge;
//						check_flag = 1;
//					}
//					else {
//						if (in < plist.get(j).plist.get(plist.get(j).abundant).in()) {
//							in = plist.get(j).plist.get(plist.get(j).abundant).in();
//							m = plist.get(j).mass;
//							cs = plist.get(j).charge;
//							check_flag = 1;
//						}
//					}
//				}
//			}
//
//			if (plist.get(j).peakNum > 0 && check_flag == 0) {
//				double misspeak_mz[] = new double[3];
//				int cnt_misspeak;
//
//				misspeak_mz[0] = (plist.get(j).mass + plist.get(j).charge*1.00235) / plist.get(j).charge;
//				misspeak_mz[1] = misspeak_mz[0] + 1 / plist.get(j).charge;
//				misspeak_mz[2] = misspeak_mz[1] + 1 / plist.get(j).charge;
//
//				if (misspeak_mz[0] > monomz*(1 - Constants.massErr) && misspeak_mz[0] < monomz*(1 + Constants.massErr)) {
//					if (in == 0) {
//						in = plist.get(j).plist.get(plist.get(j).abundant).in();
//						m = plist.get(j).mass;
//						cs = plist.get(j).charge;
//					}
//					else {
//						if (in < plist.get(j).plist.get(plist.get(j).abundant).in()) {
//							in = plist.get(j).plist.get(plist.get(j).abundant).in();
//							m = plist.get(j).mass;
//							cs = plist.get(j).charge;
//						}
//					}
//				}
//				else if (misspeak_mz[1] > monomz*(1 - Constants.massErr) && misspeak_mz[1] < monomz*(1 + Constants.massErr)) {
//					if (plist.get(j).peakNum > 1) {
//						if (in == 0) {
//							in = plist.get(j).plist.get(plist.get(j).abundant).in();
//							m = plist.get(j).mass;
//							cs = plist.get(j).charge;
//						}
//						else {
//							if (in < plist.get(j).plist.get(plist.get(j).abundant).in()) {
//								in = plist.get(j).plist.get(plist.get(j).abundant).in();
//								m = plist.get(j).mass;
//								cs = plist.get(j).charge;
//							}
//						}
//					}
//				}
//				else if (misspeak_mz[2] > monomz*(1 - Constants.massErr) && misspeak_mz[2] < monomz*(1 + Constants.massErr)) {
//					if (plist.get(j).peakNum == 3) {
//						if (in == 0) {
//							in = plist.get(j).plist.get(plist.get(j).abundant).in();
//							m = plist.get(j).mass;
//							cs = plist.get(j).charge;
//						}
//						else {
//							if (in < plist.get(j).plist.get(plist.get(j).abundant).in()) {
//								in = plist.get(j).plist.get(plist.get(j).abundant).in();
//								m = plist.get(j).mass;
//								cs = plist.get(j).charge;
//							}
//						}
//					}
//				}
//			}
//			Peak rp = new Peak();
//			rp.mdbl_mz = plist.get(j).mass;
//			rp.mdbl_intensity = plist.get(j).sum;
//			rp.mdbl_charge = plist.get(j).charge;
//			res_peak.add(rp);
//		}
//		DeconvolutedPair dp = new DeconvolutedPair();
//		dp.m = m;
//		dp.cs = cs;
//		return dp;
//	}

	public double background_intensity(List<Peak> vpk, double snt) {
		double t = 0, s = 0;
		int n = 0;

		for (int i = 0; i < (int)vpk.size(); ++i) {
			if (vpk.get(i).mdbl_intensity <= Double.MAX_VALUE && vpk.get(i).mdbl_intensity != 0) {
//				System.out.println("vpk.get(i).mdbl_intensity : " + vpk.get(i).mdbl_intensity);
//				System.out.println("vpk.get(i).in() : " + vpk.get(i).in());
				t += vpk.get(i).mdbl_intensity;
				++n;
			}
		}
		if (n > 0) {
			t /= n;	//	average
			n = 0;
			for (int i = 0; i < (int)vpk.size(); ++i) {
				if (vpk.get(i).mdbl_intensity < 5.0/*snt*/*t && vpk.get(i).mdbl_intensity != 0) {
					s += vpk.get(i).mdbl_intensity;
					++n;
				}
			}
			if (n > 0) {
				return s / n;
			}
		}

		return 0;
	}
	
	/**
	 * filter noise peaks for input of cluster. 
	 * @param profileScan is  a boolean, true if profile data, false if centroided
	 * @param pl	input peak list
	 * @param snt	signal to noise threshold
	 * @param PeakBackgroundRatio 
	 * @param thresholded this is true if FileType == FINNIGAN or MzXMLRAWDATA
	 * @param fit sets the type of peak fitting algorithm used.
	 * @param numOfGroups
	 * @return
	 */
	public List<Peak> peakPicking(boolean profileScan, List<Peak> pl,
			double snt, double PeakBackgroundRatio, boolean thresholded, PEAK_FIT_TYPE fit,
			int numOfGroups)
		{
			List<Peak> ResultPeaks = new ArrayList<>(); //Peak picking��� ����
			peakPs.setPeaksProfileType(profileScan); //Raw data�κ��� profile ����
			peakPs.setOptions(snt, 0, thresholded, fit);//����� �Է� parameter ����(snt,bgr,fittype)
			
//			System.out.println("Start Peak picking : " + profileScan);
			// GroupIndex����
			List<Integer> vGroupIndex = new ArrayList<>();
			int numPeaksinGroup = pl.size() / numOfGroups;
			if( numPeaksinGroup == 0 ) {
				numOfGroups= 1;
			}
			
			// profile ��忡�� peak profile ���ο��� partition�� ���������� �� ��� ���� �ذ��ϰ��� �ڵ� ����

			int index = 0;
			for (int Gidx = 1; Gidx < numOfGroups; Gidx++)
			{
				if (profileScan == false){
					vGroupIndex.add(numPeaksinGroup * Gidx);
				//	System.out.println("G index : " + numPeaksinGroup * Gidx);
				}
				else
				{
					index = numPeaksinGroup * Gidx;
					index++;
					
					vGroupIndex.add(index);
				}
			}
			vGroupIndex.add(pl.size() - 1);
		//	System.out.println("G index : " + (pl.size() - 1));
		//	System.out.println("G index = 0 : " + vGroupIndex.get(0) +" / "+ numPeaksinGroup);
			
			// profile ��忡�� peak profile ���ο��� partition�� ���������� �� ��� ���� �ذ��ϰ��� �ڵ� ����

			List<Peak> vGrouppk = new ArrayList<>();
			int Gidx = 0;
			for (int i = 0; i < pl.size(); i++)
			{
				vGrouppk.add(pl.get(i));
				
				//numOfGroups�� ��ŭ�� �������� spectrum�� ó���Ѵ�. 
				if (i == vGroupIndex.get(Gidx))
				{
//					System.out.println("Group peak counts : " + vGroupmz.size());
					//���� �������� peak intensity�� ����� ���ϰ� �̸� avg�� �Ѵ�.
					//avg*5�� ������ ���� intensity�� ���� peak intensity ����� background intensity�� �Ѵ�.
					double bgIntensity = background_intensity(vGrouppk, snt);

//					System.out.println("GB_intensity " + bgIntensity);
//					System.out.println("G index : " + vGroupIndex.get(Gidx));
					
					//���� ������ threshold�� �����Ѵ�. 
					double threshold = bgIntensity* PeakBackgroundRatio;

					peakPs.setPeakIntensityThreshold(threshold); //back ground intensity����
					peakPs.discoverPeaks(vGrouppk); //������ ������ peak discovery�� �����Ѵ�.
					PeakData pd = peakPs.mobj_peak_data;

//					System.out.println("result count : " + peakPs.mobj_peak_data.getNumPeaks());	
					
					// ���͸� ����� �����Ѵ�.
					for (int j = 0; j < pd.getNumPeaks(); ++j)
					{
						Peak onePeak = new Peak();
						onePeak = pd.mvect_peak_tops.get(j);

						ResultPeaks.add(onePeak);
					}

					//Group vector ����
					vGrouppk = new ArrayList<>();

					//Group index ����
					Gidx++;
					if ( Gidx != numOfGroups ) //(Gidx != 5)[]byna
						i = Math.max(0, i - 2);
					
//					System.out.println();
				}
			}

			peakPs.clear();
			return ResultPeaks;
		}
	
	public boolean isUsageRight(String[] args){
		int argc = args.length;
		
		
		if( argc != 1 ){
			System.out.println("RAPID v1.1: Ratio based Algorithm for Polypeptide Isotopic cluster Determination");
			System.out.println("Usage: <parameter file>");
			System.out.println("The <File::DataList> are required in the parameter file.");
			System.out.println("Example>");
			System.out.println("File::DataList");
			System.out.println("\tC:\\data\\18proteins.RAW");
			System.out.println("See readme.txt for other parameters.");
			return false;
		}
		return true;
	}


	public void process_rawfile(String ifname, String outputFormat , Param param, int start_idx, int end_idx) throws IOException {
		String of, of2;
		BufferedWriter fout, fout2 = null, fout3 = null, fout4 = null;
		String fileName = ifname.substring(0, ifname.lastIndexOf('.'));
		if (outputFormat.compareTo("CSV") == 0) {
			of = fileName + "_scans.csv";
			of2 = fileName + "_isos.csv";

			fout2 = new BufferedWriter(new FileWriter(new File(of2)));
		}
		else {
			of = fileName + ".pek";
		}
		fout = new BufferedWriter(new FileWriter(new File(of)));
		if (Constants.PrintPeak) {
			of2 = fileName + "_peak.csv"; // peak list ��¿�
			fout3 = new BufferedWriter(new FileWriter(new File(of2)));
		}
		if (Constants.MakeMSMS > 0) {
			of2 = fileName + ".msms"; // MSMS ����
			fout4 = new BufferedWriter(new FileWriter(new File(of2)));
		}
		
		RawData prd = ReaderFactory.getRawData(Constants.getFileType(Constants.DataFileType));
		PeakProcessing pp = new PeakProcessing();
		
		if( prd == null ){
			System.out.println("ERROR: File read error");
			if( Constants.DataFileType == 1 ){
				System.out.println("ERROR: Library to process raw data exist?");
			}
			System.exit(-1);
		}
		
		prd.load(ifname);
		
		// default parameter
		double snt = 3.0;
		double bgr = 5.0;
		double bgr2 = 0.0;
		boolean tdata = true;
		Constants.PEAK_FIT_TYPE fit = PEAK_FIT_TYPE.QUADRATIC;
		
		// load parameter
		snt = Double.parseDouble(param.get("PeakPicking::SNRThreshold"));
		bgr = Double.parseDouble(param.get("PeakPicking::BackgroundRatio"));
		bgr2 = Double.parseDouble(param.get("PeakPicking::BackgroundRatio2"));
		String fittype = param.get("PeakPicking::FitType");
		if( fittype != null ){
			if( fittype.compareTo("APEX") == 0 ){
				fit = PEAK_FIT_TYPE.APEX;
			}else if( fittype.compareTo("LORENTZIAN") == 0 ){
				fit = PEAK_FIT_TYPE.LORENTZIAN;
			}else{
				fittype = null;
			}
		}
		
		pp.setOptions(snt, 0, tdata, fit);
		
		// process all scans
		int rawtime;
		long cmt = System.currentTimeMillis();
		stime = cmt;
//		System.out.println("Start time: " + cmt);
		if ( Constants.OutputFormat.compareTo("PEK") == 0) {
			fout.append("PeakPicking::SNRThreshold = "+snt+"\n");
			fout.append("PeakPicking::BackgroundRatio = "+bgr+"\n");
			if (fittype != null)
				fout.append("PeakPicking::FitType = "+fittype+"\n");
			else
				fout.append("PeakPicking::FitType = QUADRATIC\n");
			fout.append("DeconvPep::MaxCharge = "+Constants.maxCharge+"\n");
			fout.append("DeconvPep::ThScore = "+Constants.thScore+"\n");
			fout.append("DeconvPep::Target = "+Constants.Target+"\n");
//			fout.append("DeconvPep::Truncated = %s\n", Truncated);
			fout.append("DeconvPep::ResultOrder = "+Constants.ResultOrder+"\n");
			fout.append("AdvDeconv::ScanNoModifier = "+Constants.scanNoModifier+"\n");
			fout.append("AdvDeconv::MaxAbundancePeak = "+Constants.maxAbundancePeak+"\n");
			fout.append("AdvDeconv::MaxMissPeak = "+Constants.maxMissPeak+"\n");
			fout.append("AdvDeconv::MassErr = "+Constants.massErr+"\n");
//			fout.append("AdvDeconv::ThClustExt = %.2lf\n", CDeconvPep::ThClustExt);
//			fout.append("AdvDeconv::IntsRangeErr = %.2lf\n\n", CIntsRangeFunc::FctRange);
		}
		else {
			fout.append("scan_num,scan_time,type,bpi,bpi_mz,tic,num_peaks,num_deisotoped\n");
			fout2.append("scan_num,charge,abundance,mz,fit,average_mw,monoisotopic_mw,mostabundant_mw,fwhm,signal_noise,mono_abundance,mono_plus2_abundance,RT\n");
		}
		
		if( start_idx == -1 )
			start_idx = 1;
		if (end_idx == -1)
			end_idx = prd.getNumScans();
		
		System.out.println("The number of scans: " + prd.getNumScans());
		for (int i = start_idx; i <= end_idx; ++i) {
			if( !prd.readyScan(i) ) continue;
		//	System.out.println("Read " +i);
			double rt = prd.getScanTime(i);
			if (prd.isMSScan(i)) {
//				System.out.println("Read");
				if (i % 10 == 0)
					System.out.println("scan = "+i + "/"+prd.getNumScans()  + "\t" + (System.currentTimeMillis()-stime)/1000);
				List<Peak> pks = new ArrayList<>();
				if (prd.getRawData(pks, i)) {
					cmt = System.currentTimeMillis();

					// mzXML�ϰ�� vmz���������� ��Ʈ!! 
					// �ð��� ����� �����ɸ�!! - ���Ŀ� ���� �ٽ� �߰� Ȥ�� ���ٰ�
					//if( DataFileType == 6 ) SortMZVec(vin, vmz); 

					double tic_intensity, bpi[] = {0.0}, bp_mz[] = {0.0};
					Constants.backgroundIntensity = background_intensity(pks, snt);

					tic_intensity = get_TIC(pks.get(0).mz(), pks.get(pks.size()-1).mz(), pks,
						(float)(Constants.backgroundIntensity * bgr), bpi, bp_mz);

					List<Peak> pkl = new ArrayList<>();

//new peak picking
//					System.out.println("Start peakPicking\t" + (System.currentTimeMillis()-stime)/1000);
					
					if (bgr >0)
						pkl = peakPicking(prd.isProfileScan(i), pks, snt, bgr, tdata, fit, Constants.Partition);
					else {
						for (int j = 0; j<pks.size(); j++) {
							pkl.add(pks.get(j));
						}
					}
//					System.out.println("End peakPicking\t" + (System.currentTimeMillis()-stime)/1000);
					
					if (Constants.PrintPeak) {
						fout3.append("Scan = "+(i + Constants.scanNoModifier)+"\n");
						for (int j = 0; j<pkl.size(); j++) {
							fout3.append(String.format("%.6f", pkl.get(j).mz()) + "," + String.format("%.6f", pkl.get(j).in())+"\n");
						}
						fout3.append("\n");
					}

//old peak picking
/*					pp.SetPeakIntensityThreshold(CDeconvPep::BackgroundIntensity*bgr);
					pp.SetPeaksProfileType(prd->IsProfileScan(i));
//					if( DataFileType == 6 ) pp.SetPeaksProfileType(false); // mxXML�ϰ�� profile����

					pp.DiscoverPeaks(&vmz, &vin);
					Engine::PeakProcessing::PeakData & pd = *(pp.mobj_peak_data);
//					fprintf(fout2, "Scan = %d\n", i + CDeconvPep::ScanNoModifier);
					for(int j = 0; j < pd.GetNumPeaks(); ++j) {
					CPeak pk(&pd.mvect_peak_tops[j]);
					pkl.push_back(pk);
//						fprintf(fout2, "%d\t%lf\t%lf\n", j, pkl.get(j).mz(), pkl.get(j).in());
					}
*/
//					System.out.println("Start Clustering\t" + (System.currentTimeMillis()-stime)/1000);
					List<IsoCluster> vic = dp.process(pkl, pks);
//					System.out.println("End Clustering\t" + (System.currentTimeMillis()-stime)/1000);
//					System.out.println("Sort Clusters\t" + (System.currentTimeMillis()-stime)/1000);
					Collections.sort(vic);
//					System.out.println("End Sorting\t" + (System.currentTimeMillis()-stime)/1000);

					if (Constants.OutputFormat.compareTo("PEK") == 0){
						//print head
						fout.append("Scan = "+(i + Constants.scanNoModifier)+"\n");
						fout.append("Scan Type = ");
						if (prd.isMSScan(i)) fout.append("MS\n");
						else fout.append("MS/MS\n");
						fout.append("Peak = "+bp_mz[0]+"\n");
						fout.append("RT = "+prd.getScanTime(i)+"\n");
						fout.append("NL = "+bpi[0]+"\n");
						fout.append("TIC = "+tic_intensity+"\n");
						int tm_hour = ((int)(cmt/3600000));
						int tm_min = ((int)(cmt/60000)) - tm_hour*60;
						int tm_sec = ((int)(cmt/1000)) - tm_hour*3600 - tm_min*60;
						fout.append("Processing start time: "+tm_hour+":"+tm_min+":"+tm_sec+"\n");
						fout.append("Isotopic mass transform results:\n");
						fout.append("CS,\tAbundance,\tm/z,\tScore,\tAverage MW,\tMonoisotopic MW,\tMost abundant MW\n");

						//print body
						for (int j = 0; j < (int)vic.size(); ++j) {
							fout.append(vic.get(j).charge + "\t");
							fout.append(vic.get(j).sum + "\t");
							fout.append(vic.get(j).plist.get(vic.get(j).abundant).mz() + "\t");
							fout.append(vic.get(j).score + "\t");
							fout.append(vic.get(j).avg + "\t");
							fout.append(vic.get(j).mass + "\t");
							fout.append(vic.get(j).plist.get(vic.get(j).abundant).mass(vic.get(j).charge) + "\n");
						}

						//print tail
						cmt = System.currentTimeMillis();
						tm_hour = ((int)(cmt/3600000));
						tm_min = ((int)(cmt/60000)) - tm_hour*60;
						tm_sec = ((int)(cmt/1000)) - tm_hour*3600 - tm_min*60;
						fout.append("Processing stop time: "+tm_hour+":"+tm_min+":"+tm_sec+"\n");
						fout.append("Number of peaks in spectrum = "+pkl.size()+"\n");
						fout.append("Number of isotopic distributions detected = "+vic.size()+"\n");
						fout.append("\n");
					}
					else {
						fout.append((i + Constants.scanNoModifier) + ","+String.format("%.6f", prd.getScanTime(i))+","+prd.getMSLevel(i)+","+String.format("%.6f", bpi[0])+","+String.format("%.6f", bp_mz[0])+","+String.format("%.6f", tic_intensity)+","+pkl.size()+","+vic.size()+"\n");
						for (int j = 0; j < vic.size(); ++j) {
							fout2.append((i + Constants.scanNoModifier)+","+vic.get(j).charge+","+String.format("%.6f", vic.get(j).sum)+","+String.format("%.6f", vic.get(j).plist.get(vic.get(j).abundant).mz())+","+String.format("%.6f", vic.get(j).score)+","
									+String.format("%.6f", vic.get(j).avg)+","+String.format("%.6f", vic.get(j).mass)+","+String.format("%.6f", vic.get(j).plist.get(vic.get(j).abundant).mass(vic.get(j).charge))+",0,"
									+String.format("%.6f", (vic.get(j).plist.get(vic.get(j).abundant).in() / Constants.backgroundIntensity))+","
									+String.format("%.6f", ((vic.get(j).peakNum == 0) ? vic.get(j).plist.get(0).in() : 0))+","
									+String.format("%.6f", ((vic.get(j).plist.size() + vic.get(j).peakNum > 2 && vic.get(j).peakNum < 3) ? vic.get(j).plist.get(2 - vic.get(j).peakNum).in() : 0)));
							fout2.append(","+rt+"\n");
						}
					}
					
//					System.out.println("Deconvolution\t" + (System.currentTimeMillis()-stime)/1000);
					int loop_flag = 0;
					while (Constants.MakeMSMS > 0 && i + 1 <= prd.getNumScans()) {
						if( !prd.readyScan(i+1) ) continue;
						
						if (!prd.isMSScan(i + 1)) {
							i++;
//							if (i % 1000 == 0)
//								System.out.println("scan = "+i + "\t" + System.currentTimeMillis());

							
//							if( i == 5152 ) {
//								System.out.println("Scan 5152");
//							}
							
							//determine precursor m/z
//							DeconvolutedPair dp = determineMono(monomz, vic);
							double m = 0, in = 0;
							int cs = 0;
							List<Peak> res_peak = new ArrayList<>();

							double monomz;
							if (Constants.DataFileType == 6) {
								monomz = prd.getParentMz(i);
								//cs = prd->GetParentCharge(i);								
							}
							else {
								monomz = prd.getMonoMZFromHeader(i);
							}
							
							//determine precursor m/z
							for (int j = 0; j < vic.size(); ++j) {
								int check_flag = 0;
								for (int k = 0; k<vic.get(j).plist.size(); k++) {
//									monomz = plist.get(j).plist.get(k).mdbl_mz;
									if (vic.get(j).plist.get(k).mz() > monomz*(1 - Constants.massErr/*1.0E-05*/) && vic.get(j).plist.get(k).mz() < monomz*(1 + Constants.massErr) && in < vic.get(j).sum) {
										if (in == 0) {
											in = vic.get(j).plist.get(vic.get(j).abundant).in();
											m = vic.get(j).mass;
											cs = vic.get(j).charge;
											check_flag = 1;
										}
										else {
											if (in < vic.get(j).plist.get(vic.get(j).abundant).in()) {
												in = vic.get(j).plist.get(vic.get(j).abundant).in();
												m = vic.get(j).mass;
												cs = vic.get(j).charge;
												check_flag = 1;
											}
										}
									}
								}

								if (vic.get(j).peakNum > 0 && check_flag == 0) {
									double misspeak_mz[] = new double[3];
//									int cnt_misspeak;

									misspeak_mz[0] = (vic.get(j).mass + vic.get(j).charge*1.00235) / vic.get(j).charge;
									misspeak_mz[1] = misspeak_mz[0] + 1 / vic.get(j).charge;
									misspeak_mz[2] = misspeak_mz[1] + 1 / vic.get(j).charge;

									if (misspeak_mz[0] > monomz*(1 - Constants.massErr) && misspeak_mz[0] < monomz*(1 + Constants.massErr)) {
										if (in == 0) {
											in = vic.get(j).plist.get(vic.get(j).abundant).in();
											m = vic.get(j).mass;
											cs = vic.get(j).charge;
										}
										else {
											if (in < vic.get(j).plist.get(vic.get(j).abundant).in()) {
												in = vic.get(j).plist.get(vic.get(j).abundant).in();
												m = vic.get(j).mass;
												cs = vic.get(j).charge;
											}
										}
									}
									else if (misspeak_mz[1] > monomz*(1 - Constants.massErr) && misspeak_mz[1] < monomz*(1 + Constants.massErr)) {
										if (vic.get(j).peakNum > 1) {
											if (in == 0) {
												in = vic.get(j).plist.get(vic.get(j).abundant).in();
												m = vic.get(j).mass;
												cs = vic.get(j).charge;
											}
											else {
												if (in < vic.get(j).plist.get(vic.get(j).abundant).in()) {
													in = vic.get(j).plist.get(vic.get(j).abundant).in();
													m = vic.get(j).mass;
													cs = vic.get(j).charge;
												}
											}
										}
									}
									else if (misspeak_mz[2] > monomz*(1 - Constants.massErr) && misspeak_mz[2] < monomz*(1 + Constants.massErr)) {
										if (vic.get(j).peakNum == 3) {
											if (in == 0) {
												in = vic.get(j).plist.get(vic.get(j).abundant).in();
												m = vic.get(j).mass;
												cs = vic.get(j).charge;
											}
											else {
												if (in < vic.get(j).plist.get(vic.get(j).abundant).in()) {
													in = vic.get(j).plist.get(vic.get(j).abundant).in();
													m = vic.get(j).mass;
													cs = vic.get(j).charge;
												}
											}
										}
									}
								}
								Peak rp = new Peak();
								rp.mdbl_mz = vic.get(j).mass;
								rp.mdbl_intensity = vic.get(j).sum;
								rp.mdbl_charge = vic.get(j).charge;
								res_peak.add(rp);
							}
//							if( i == 5152 ) {
//								System.out.println("Scan 5152");
//								System.out.println("loop_flag: " + loop_flag);
//								System.out.println("cs: " + cs);
//								System.out.println("m: " + m);
//							}
							if (loop_flag == 1) {	//mzXML�� charge�� �������� ���� ��� charge�� 2,3�����ؼ� �ι� ���
								cs = 3;
								loop_flag = 2;
								m = monomz*cs - cs*1.00235;
							}

							if (cs == 0) {
								cs = 2;
								loop_flag = 1;
								m = monomz*cs - cs*1.00235;
							}

							List<Peak> vpk2 = new ArrayList<>();
							if (prd.getRawData(vpk2, i)) {
								fout4.append("Scan = "+(i + Constants.scanNoModifier)+"\n");
								fout4.append(String.format("%.6f", m) + "\t"+cs+"\n");
								if (Constants.MakeMSMS == 1) {			//msms scan�� peak�� ���
									for (int j = 0; j<vpk2.size(); j++) {
										fout4.append(String.format("%.6f", vpk2.get(j).mz()) + "\t"+String.format("%.6f", vpk2.get(j).in())+"\n");
									}
								}
								else {						//msms Deisotoped
									Constants.backgroundIntensity = background_intensity(vpk2, snt);
									List<Peak> pkl2 = new ArrayList<>();
									if (bgr2 >= 0)			//if (bgr2 > 0)  ->  if (bgr2 >= 0)  by slee
										pkl2 = peakPicking(prd.isProfileScan(i), vpk2, snt, bgr2, tdata, fit, Constants.Partition);
									else {
										for (int j = 0; j<vpk2.size(); j++) {
											pkl2.add(vpk2.get(j));
										}
									}
									//printf("# of msms pick = %d\n", pkl2.size());
									DeconvPept dp2 = new DeconvPept(pkl2, vpk2);
									List<IsoCluster> vic2 = dp2.process(pkl2, vpk2);

									for (int j = 0; j<vic2.size(); j++) {
										fout4.append(String.format("%.6f", vic2.get(j).mass + 1.0078246) + "\t"+String.format("%.6f", vic2.get(j).sum)+"\t"+vic2.get(j).charge+"\n");	//monoisotopic_mw, sum of intensity				
									}
								}
								fout4.append("\n");
							}
							else {
								System.out.println("ERROR: GetRawData at scan "+i+"\n\n");
							}

						} else {
							break;
						}

						if (loop_flag == 1) {
							i--;
						}
						else if (loop_flag == 2) {
							loop_flag = 0;
						}
					}
//					System.out.println("End deconvolution\t" + (System.currentTimeMillis()-stime)/1000);

				} else {
					System.out.println("ERROR: GetRawData at scan "+i+"\n\n");
				}
			}
		}
		fout.close();
		fout2.close();
		if (Constants.PrintPeak) {
			fout3.close();
		}
		if (Constants.MakeMSMS > 0) {
			fout4.close();
		}
	}
	
}
