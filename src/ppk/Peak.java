package ppk;

import utils.Constants;

public class Peak {
	
	public Peak(){
		mdbl_mz = 0;
		mdbl_intensity = 0;

		/* charge 결정하는 부분 추가 형서*/
		//				mdbl_charge = 0;
		/* charge 결정하는 부분 추가 형서*/

		mdbl_SN = 0;
		mint_peak_index = -1;
		mint_data_index = -1;
		mdbl_FWHM = 0;
	}
	
	public Peak(Peak a){
		mdbl_mz = a.mdbl_mz;
		mdbl_intensity = a.mdbl_intensity;

		/* charge 결정하는 부분 추가 형서*/
		//				mdbl_charge = a.mdbl_charge;
		/* charge 결정하는 부분 추가 형서*/

		mdbl_SN = a.mdbl_SN;
		mint_peak_index = a.mint_peak_index;
		mint_data_index = a.mint_data_index;
		mdbl_FWHM = a.mdbl_FWHM;
	}
	
	//! mz of the peak.
	public	double mdbl_mz;
	//!  intensity of peak.
	public	double mdbl_intensity;


	/* charge 결정하는 부분 추가 형서*/
	public double mdbl_charge;
	/* charge 결정하는 부분 추가 형서*/

	//! Signal to noise ratio
	public double mdbl_SN;
	//! index in PeakData::mvect_peak_tops std::vector. 
	public int mint_peak_index;
	//! index in mzs, intensity vectors that were used to create the peaks in PeakProcessor::DiscoverPeaks.
	public int mint_data_index;
	//! Full width at half maximum for peak.
	public double mdbl_FWHM;

	/**
	 *  Assignment operator
	 * @param a Peak we want to assign to this Peak.
	 * @return
	 */
	public Peak operator(Peak a){
		mdbl_mz = a.mdbl_mz;
		mdbl_intensity = a.mdbl_intensity;

		/* charge 결정하는 부분 추가 형서*/
		//				mdbl_charge = a.mdbl_charge;
		/* charge 결정하는 부분 추가 형서*/

		mdbl_SN = a.mdbl_SN;
		mint_peak_index = a.mint_peak_index;
		mint_data_index = a.mint_data_index;
		mdbl_FWHM = a.mdbl_FWHM;
		return this;
	}
	
	/**
	 *  Sets the members of the Peak.
	 * @param mz m/z of the peak.
	 * @param intensity intensity of the peak
	 * @param signal2noise  signal2noise of the peak look at PeakProcessor::PeakStatistician::FindSignalToNoise
	 * @param peak_idx index of the peak in PeakData::mvect_peak_tops std::vector of the PeakData instance that was used to generate these peaks.
	 * @param data_idx index of the peak top in the mz, intensity lists that are the raw data input into PeakData::DiscoverPeaks
	 * @param fwhm full width half max of the peak. For details about how this is calculated look at PeakProcessor::PeakStatistician::FindFWHM.
	 */
	public void set(double mz, double intensity, double signal2noise, int peak_idx, int data_idx, double fwhm){
		mdbl_mz = mz;
		mdbl_intensity = intensity;

		/* charge 결정하는 부분 추가 형서*/
		//				mdbl_charge = charge;
		/* charge 결정하는 부분 추가 형서*/

		mdbl_SN = signal2noise;
		mint_peak_index = peak_idx;
		mint_data_index = data_idx;
		mdbl_FWHM = fwhm;
	}

	/* charge 결정하는 부분 추가 형서*/
	//			void Set(double mz, double intensity,double charge ,double signal2noise, int peak_idx, int data_idx, double fwhm);
	/* charge 결정하는 부분 추가 형서*/

	/**
	 * @return	mz of the peak
	 */
	public double mz() {
		return mdbl_mz;
	}

	/**
	 * @return  intensity of the peak
	 */
	public double in() {
		return mdbl_intensity;
	}

	/**
	 * @return  Signal to Noise ratio
	 */
	public double sn() {
		return mdbl_SN;
	}

	/**
	 * 
	 * @return Full width at half maximum of the peak
	 */
	public double fwhm() {
		return mdbl_FWHM;
	}

	public int idx() {
		return mint_peak_index;
	}

	/**
	 * @param chrg is charge of peak
	 * @return  mass of given charge
	 */
	public double mass(int chrg) {		//이상해서 고쳐봄 by slee
		if (chrg == 0) {
			return mz();
		}
		else {
			return (mz() - (Constants.MASS_H - Constants.MASS_ELECTRON))*chrg;
		}
	}
	
}
