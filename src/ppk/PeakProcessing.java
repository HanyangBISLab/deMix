package ppk;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import utils.Constants.*;

public class PeakProcessing {
	
	//! tolerance width used while searching in attention list in the peaks. 
	private double mdbl_mz_tolerance_attention_list;
	
	//! flag to check attention list for peaks that we want to necessarily look for.
	private boolean mbln_chk_attention_list;
	
	//! map of attention list of m/z's. 
	/*!
	\remark This was implemented as a map because it would allow for direct sorted lookup for an m/z value.
	*/
	private Map<Double, Integer> mmap_attention_list;

	//! background intensity. When user sets min signal to noise, and min intensity, 
	// this value is set as min intensity / min signal to noise. 
	private double mdbl_background_intensity;

	//! minimum intensity for a point to be considered a peak. 
	private double mdbl_peak_intensity_threshold;
	
	//! signal to noise threshold for a Peak to be considered as a peak.
	private double mdbl_signal_2_noise_threshold;
	
	//! fit type used in the PeakFit class to detect peaks in the raw data object. 
	private PEAK_FIT_TYPE menm_peak_fit_type;
	private PEAK_PROFILE_TYPE menm_profile_type;

	//! This variable helps find the m/z value of a peak using the specified fit function. 
	private PeakFit mobj_peak_fit;


	//! tests if a specific m/z value is in the attention list.
	/*!
	\param mz_val is the m/z value we are looking for. 
	\return returns true if mz_val is found in attention list, and false otherwise.
	\note The tolerance used in searching is specified in PeakProcessor::mdbl_mz_tolerance_attention_list
	*/
	private boolean isInAttentionList(double mz_val){
		if (mmap_attention_list.isEmpty())
			return false;
//		(
		List<Double> mmap_key_list = new ArrayList<>();
		mmap_key_list.addAll(mmap_attention_list.keySet());
		Collections.sort(mmap_key_list);
		Iterator<Double> closest_iterator = mmap_key_list.iterator();
		double tmp  = closest_iterator.next();
		while( tmp < mz_val - mdbl_mz_tolerance_attention_list ) tmp  = closest_iterator.next();
		if (tmp <= mz_val + mdbl_mz_tolerance_attention_list)
			return true;
		return false;
	}

	//! Is the data centroided ? 
	private boolean mbln_centroided_data;
	
	// if the data is thresholded, the ratio is taken as the ratio to background intensity.
	private boolean mbln_thresholded_data;

	//! pointer to PeakData instance that stores the peaks found by an instance of this PeakProcessor.
	public PeakData mobj_peak_data;
	//! pointer to PeakStatistician instance that is used to calculate signal to noise and full width at half maximum for the peaks in the raw data. 
	public PeakStatistician mobj_peak_statistician;
	
	/**
	 * sets the threshold for signal to noise for a peak to be considered as real.
	 * @param s2n is the signal to noise threshold value.
	 * remarks For a peak to be considered real it has to pass two criterias:
	- Its signal to noise must be greater than the threshold (mdbl_signal_2_noise_threshold)
	- Its intensity needs to be greater than the threshold (mdbl_peak_intensity_threshold)
	 */
	public void setSignal2NoiseThreshold(double s2n){
		this.mdbl_signal_2_noise_threshold = s2n;
		if (mbln_thresholded_data)
		{
			if (mdbl_signal_2_noise_threshold != 0)
			{
				mdbl_background_intensity = mdbl_peak_intensity_threshold / mdbl_signal_2_noise_threshold;
			}
			else
				mdbl_background_intensity = 1;
		}
	}
	
	/**
	 * sets the threshold intensity for a peak to be considered a peak.
	 * @param threshold is the threshold peak intensity.
	 * If threshold is less than zero, then the Helpers::absolute value of the threshold is used as the cutoff intensity.
	However, if threshold is greater than equal to zero, otherwise it is proportional to threshold * background intensity in scan.
	remarks For a peak to be considered real it has to pass two criterias:
	Its signal to noise must be greater than the threshold (mdbl_signal_2_noise_threshold)
	Its intensity needs to be greater than the threshold (mdbl_peak_intensity_threshold)
	calculatePeakIntensityThreshold
	 */
	public void setPeakIntensityThreshold(double threshold){
		this.mdbl_peak_intensity_threshold = threshold;
		if (mbln_thresholded_data)
		{
			if (mdbl_signal_2_noise_threshold != 0)
				mdbl_background_intensity = threshold / mdbl_signal_2_noise_threshold;
			else if (threshold != 0)
				mdbl_background_intensity = threshold;
			else
				mdbl_background_intensity = 1;
		}

		//			printf("mdbl_background_intensity = %f \n",mdbl_background_intensity);
	}

	/**
	 * sets the type of peak fitting used to find m/z values for peaks. 
	 * @param type specifies the type of peak fitting.
	 */
	public void setPeakFitType(PEAK_FIT_TYPE type){
		this.menm_peak_fit_type = type;
		mobj_peak_fit.setOptions(type);
	}
	
	//! sets the type of profile 
	/*!
	\param profile 
	*/
	/**
	 * 
	 * @param profile is  a boolean, tru if profile data, false if centroided
	 */
	public void setPeaksProfileType(boolean profile){
		if (profile == true)
			menm_profile_type = PEAK_PROFILE_TYPE.PROFILE;
		else
			menm_profile_type = PEAK_PROFILE_TYPE.CENTROIDED;
	}

	/**
	 *  sets the options for this instance.
	 * @param s2n sets the threshold signal to noise value.
	 * @param thresh sets the peak intensity threshold.
	 * @param thresholded
	 * @param type sets the type of peak fitting algorithm used.
	 */
	public void setOptions(double s2n, double thresh, boolean thresholded, PEAK_FIT_TYPE type){
		mbln_thresholded_data = thresholded;
		// signal to noise should ideally be set before PeakIntensityThreshold
		setSignal2NoiseThreshold(s2n);
		setPeakIntensityThreshold(thresh);
		setPeakFitType(type);
	}
	
	/** 
	 * clears the PeakData member variable PeakProcessor::mobj_peak_data and  PeakProcessor::mmap_attention_list
	 */
	public void clear(){
		mmap_attention_list = new TreeMap<>();
		mobj_peak_data.clear();
	}
	
	/**
	 *  Function discovers peaks in the m/z and intensity vectors supplied. 
	 * @param vect_mz is the pointer to std::vector of m/z values
	 * @param vect_intensity intensity is the pointer to std::vector of intensity values
	 * @return returns the number of peaks that were found in the vectors.
	 * remarks The function uses PeakStatistician::findFWHM, and PeakStatistician::findSignalToNoise functions
	to discover the full width at half maximum and signal to noise values for a peak. The signal to noise of a
	peak is tested against the threshold value before its accepted as a peak. All peaks are used during the process,
	but once generated only those which are above mdbl_peak_intensity_threshold are tested for peptidicity by Deconvolution::HornMassTransform
	 */
	public int discoverPeaks(List<Peak> vect_pk){
		if (vect_pk.size() == 0)
			return 0;
		double min_mz = vect_pk.get(0).mz();

		int num_elements = vect_pk.size() - 1;
		double max_mz = vect_pk.get(num_elements).mz();

		return discoverPeaks(vect_pk, min_mz, max_mz);
	}
	
	//! Function discovers peaks in the m/z and intensity vectors supplied within the supplied m/z window. 
	/*!
	\param vect_mz is the pointer to std::vector of m/z values
	\param vect_intensity is the pointer to std::vector of intensity values
	\param start_mz minimum m/z of the peak.
	\param stop_mz maximum m/z of the peak.
	\return returns the number of peaks that were found in the vectors.
	\remarks The function uses PeakStatistician::FindFWHM, and PeakStatistician::FindSignalToNoise functions
	to discover the full width at half maximum and signal to noise values for a peak. The signal to noise of a
	peak is tested against the threshold value before its accepted as a peak. All peaks are used during the process,
	but once generated only those which are above mdbl_peak_intensity_threshold are tested for peptidicity by Deconvolution::HornMassTransform
	*/
	public int discoverPeaks(List<Peak> vect_pk, double start_mz, double stop_mz){
		if (vect_pk.size() < 1)
			return 0;

		mobj_peak_data.clear();
		int num_data_pts = vect_pk.size();
//		System.out.println("num_data_pts before: " + num_data_pts);
		/* charge 결정하는 부분 추가 형서*/
		//			vector<double> vmz, vin, vsc; // vsc는 score저장
		//			vector<int> vid;
		//			vector<short> vch; // charge 저장
		//			vmz.clear();	vin.clear();	vid.clear();	vch.clear();	vsc.clear();
		//			vmz.reserve(num_data_pts);	vin.reserve(num_data_pts);	vid.reserve(num_data_pts);	vch.reserve(num_data_pts);	
		//			vsc.reserve(num_data_pts);

		/* charge 결정하는 부분 추가 형서*/

		int ilow;
		int ihigh;

		
		
		int start_index = PeakIndex.getNearestBinary(vect_pk, start_mz, 0, num_data_pts - 1);
//		System.out.println("???");
		int stop_index = PeakIndex.getNearestBinary(vect_pk, stop_mz, start_index, num_data_pts - 1);
//		System.out.println("num_data_pts after: " + num_data_pts);
//		System.out.println("Discover peaks : " + start_mz + "\t"+ stop_mz);
		
		if (start_index <= 0)
			start_index = 1;
		if (stop_index >= vect_pk.size() - 2)
			stop_index = vect_pk.size() - 2;

		/* charge 결정하는 부분 추가 형서*/
		/*
		for (int i=0 ; i < num_data_pts ; i++)
		{
		vmz.push_back((*vect_mz)[i]);
		vin.push_back((*vect_intensity)[i]);
		vid.push_back(i);
		}

		int mIdx = 0;
		double t1,t2 = 0.0;
		int t3 = 0;



		for(int i=0;i<vmz.size();i++)   // intensity 값을 기준으로 sort 내림차순, vid는 빠른 인덱싱을 위해
		{
		t1 = vmz[i]; t2 = vin[i];  t3 = vid[i];
		mIdx = i;
		for(int j=i;j<vmz.size();j++)
		{
		if (vin[mIdx]<=vin[j]) mIdx = j;
		}

		vmz[i] = vmz[mIdx]; vin[i] = vin[mIdx];	vid[i] = vid[mIdx];
		vmz[mIdx] = t1;		vin[mIdx] = t2;		vid[mIdx] = t3;
		}

		double score;
		double current_mz, last_mz, next_mz, compare_A, compare_B;
		double current_mass;
		double expect_mz[11]={0};
		short max_p_c;


		for(int i=0; i<num_data_pts; i++)
		{
		current_mz = (*vect_mz)[vid[i]];

		if(vid[i]== 0)
		{
		last_mz = 0;
		next_mz = (*vect_mz)[vid[i]+1];
		}

		else if(vid[i]== num_data_pts-1)
		{
		last_mz = (*vect_mz)[vid[i]-1];
		next_mz = 0;
		}

		else
		{
		last_mz = (*vect_mz)[vid[i]-1];
		next_mz = (*vect_mz)[vid[i]+1];
		}

		compare_A=current_mz-last_mz;
		compare_B=next_mz-current_mz;

		max_p_c = zmax_determination((compare_A <= compare_B) ? compare_A : compare_B);

		if(max_p_c==0 || max_p_c==1)
		vch.push_back(max_p_c);
		else
		{
		for(int k=max_p_c; k>=1; k--)
		{
		current_mass=current_mz*k;
		score=0.0;
		int cnt=-5;
		for(int j=0; j<=10; j++)
		{
		expect_mz[j]=(current_mass+cnt)/k;
		int test_index=mobj_pk_index.GetNearestBinary(*vect_mz,expect_mz[j], 0, num_data_pts-1) ;
		double distance=Helpers::absolute(expect_mz[j]-(*vect_mz)[test_index]);
		if (distance==0.0)
		score+=1000.0;
		else
		{
		score+=0.6/distance;
		}
		cnt++;
		}
		vsc.push_back(score);
		}

		double max_score=0.0;
		int max_index=0;

		for(int m=(vsc.size()-1); m>=0; m--)
		{
		if(max_score<=vsc[m])
		{
		max_score=vsc[m];
		max_index=m;
		}
		vsc.pop_back();
		}
		vch.push_back(max_p_c-max_index);

		}

		}

		short t4=0;

		for(int i=0;i<vmz.size();i++)   // vmz 값을 기준으로 sort 내림차순, vid는 빠른 인덱싱을 위해
		{
		t1 = vmz[i]; t2 = vin[i];  t3 = vid[i];  t4 = vch[i];
		mIdx = i;
		for(int j=i;j<vmz.size();j++)
		{
		if (vmz[mIdx]>=vmz[j]) mIdx = j;
		}

		vmz[i] = vmz[mIdx]; vin[i] = vin[mIdx];	vid[i] = vid[mIdx];	vch[i] = vch[mIdx];
		vmz[mIdx] = t1;		vin[mIdx] = t2;		vid[mIdx] = t3;		vch[mIdx] = t4;
		}

	

		*/
		
		for (int index = start_index; index <= stop_index; index++)
		{
			Peak peak = new Peak();
			double FWHM = -1;
			double current_intensity = vect_pk.get(index).in();
			double last_intensity = vect_pk.get(index-1).in();
			double next_intensity = vect_pk.get(index+1).in();
			double current_mz = vect_pk.get(index).mz();
//			System.out.println(current_intensity + " vs. " + mdbl_peak_intensity_threshold);
			if (menm_profile_type == PEAK_PROFILE_TYPE.CENTROIDED)
			{
				
				if (current_intensity >= mdbl_peak_intensity_threshold)
				{
//					System.out.println("Do this");
					double mass_ = vect_pk.get(index).mz();
					double SN = current_intensity / mdbl_peak_intensity_threshold;
					FWHM = 0.6;

					peak.set(mass_, current_intensity, SN, mobj_peak_data.getNumPeaks(), index, FWHM);
					mobj_peak_data.addPeak(peak);
				}
			}

			else if (menm_profile_type == PEAK_PROFILE_TYPE.PROFILE)
			{
				if (current_intensity >= last_intensity && current_intensity >= next_intensity
					&& current_intensity >= this.mdbl_peak_intensity_threshold)
				{
					//See if the peak meets the conditions.
					//The peak data will be found at _transformData->begin()+i+1.
					double SN = 0;

					if (!mbln_thresholded_data)
						SN = this.mobj_peak_statistician.findSignalToNoise(current_intensity, (vect_pk), index);
					else
						SN = current_intensity / mdbl_background_intensity;

					// Run Full-Width Half-Max algorithm to try and squeak out a higher SN
					if (SN < this.mdbl_signal_2_noise_threshold)
					{
						// 수정
						System.out.printf("SN = %lf mdbl_signal_2_noise_threshold = %lf ", SN, this.mdbl_signal_2_noise_threshold);
//						double mass_ = vect_pk.get(index).mz();
						FWHM = this.mobj_peak_statistician.findFWHM(vect_pk, index, SN);
						if (FWHM > 0 && FWHM < 0.5)
						{
							ilow = PeakIndex.getNearestBinary(vect_pk, current_mz - FWHM, 0, index);
							ihigh = PeakIndex.getNearestBinary(vect_pk, current_mz + FWHM, index, stop_index);

							double low_intensity = vect_pk.get(ilow).in();
							double high_intensity = vect_pk.get(ihigh).in();

							double sum_intensity = low_intensity + high_intensity;
							if (sum_intensity > 0)
								SN = (2.0 * current_intensity) / sum_intensity;
							else
								SN = 10;
						}
					}
					// Found a peak, make sure it is in the attention list, if there is one.
					if (SN >= this.mdbl_signal_2_noise_threshold && (!this.mbln_chk_attention_list || this.isInAttentionList(current_mz)))
					{
						// Find a more accurate m/z location of the peak.
						double fittedPeak = mobj_peak_fit.fit(index, vect_pk);
						if (FWHM == -1)
						{
							FWHM = this.mobj_peak_statistician.findFWHM(vect_pk, index, SN);
						}

						if (FWHM > 0)
						{
							peak.set(fittedPeak, current_intensity, SN, mobj_peak_data.getNumPeaks(), index, FWHM);
							mobj_peak_data.addPeak(peak);
						}
						// move beyond peaks have the same intensity.
						while (index < num_data_pts && vect_pk.get(index).in() == current_intensity)
							index++;
					}
				}
			}
		}
		mobj_peak_data.mptr_vect_peaks = vect_pk;

		/*  charge 출력 부분
		vector<short>::iterator iter1;

		iter1=vc.begin();
		int a=0;

		for(iter1=vc.begin(); iter1!=vc.end(); iter1++)
		{
		printf("charge = %d\t", *iter1);
		a++;
		}
		printf(" 마지막 값은 = %d\n", a);
		*/

		return mobj_peak_data.getNumPeaks();
	}

	// charge 결정에 필요한 함수 형서
	public short zmax_determination(double charge){
		short CS;
		if (charge >= 1.0)
		{
			CS = 1;
			return CS;
		}
		else if (charge<1.0 && charge >= 0.5)
		{
			CS = 2;
			return CS;
		}
		else if (charge<0.5 && charge >= 0.33)
		{
			CS = 3;
			return CS;
		}
		else if (charge<0.33 && charge >= 0.25)
		{
			CS = 4;
			return CS;
		}
		else if (charge<0.25 && charge >= 0.2)
		{
			CS = 5;
			return CS;
		}
		else if (charge<0.2 && charge >= 0.16)
		{
			CS = 6;
			return CS;
		}
		else if (charge<0.16 && charge >= 0.14)
		{
			CS = 7;
			return CS;
		}
		else if (charge<0.14 && charge >= 0.125)
		{
			CS = 8;
			return CS;
		}
		else if (charge<0.125 && charge >= 0.11)
		{
			CS = 9;
			return CS;
		}
		else if (charge<0.11 && charge >= 0.1)
		{
			CS = 10;
			return CS;
		}
		else if (charge<0.1 && charge >= 0.09)
		{
			CS = 11;
			return CS;
		}
		else if (charge<0.09 && charge >= 0.083)
		{
			CS = 12;
			return CS;
		}
		else if (charge<0.083 && charge >= 0.076)
		{
			CS = 13;
			return CS;
		}
		else if (charge<0.076 && charge >= 0.071)
		{
			CS = 14;
			return CS;
		}
		else if (charge<0.071 && charge >= 0.066)
		{
			CS = 15;
			return CS;
		}
		else if (charge<0.066 && charge >= 0.0625)
		{
			CS = 16;
			return CS;
		}
		else if (charge<0.0625 && charge >= 0.059)
		{
			CS = 17;
			return CS;
		}
		else if (charge<0.059 && charge >= 0.055)
		{
			CS = 18;
			return CS;
		}
		else if (charge<0.055 && charge >= 0.052)
		{
			CS = 19;
			return CS;
		}
		else if (charge<0.052 && charge >= 0.05)
		{
			CS = 20;
			return CS;
		}
		else
			return -1;
	}
	
	public PeakProcessing() {
		mobj_peak_statistician = new PeakStatistician();
		mmap_attention_list = new TreeMap<>();
		mdbl_mz_tolerance_attention_list = 5.0;
		mbln_chk_attention_list = false;
		mobj_peak_data = new PeakData();
		mbln_centroided_data = false;
		menm_profile_type = PEAK_PROFILE_TYPE.PROFILE;
		mobj_peak_fit = new PeakFit();
		setPeakFitType(PEAK_FIT_TYPE.QUADRATIC);
	}
}
