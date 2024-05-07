package ppk;

import java.util.List;

import utils.Constants;
import utils.Constants.PEAK_FIT_TYPE;

/**
 * Used for detecting peaks in the data. 
 * @author user_Lee
 *
 */
public class PeakFit {
	/**
	 *  Type of fit function used to find the peaks
	 */
	private PEAK_FIT_TYPE menm_peak_fit_type;
	
	/**
	 *  member variable to find out information about peaks such as signal to noise and full width at half maximum.
	 */
	private PeakStatistician mobj_peak_statistician;

	/**
	 *  Gets the peak that fits the point at a given index with a quadratic fit.
	 * @param mzs List of raw data of m\zs.
	 * @param intensities List of raw data of intensities.
	 * @param index index of the point in the m/z lists which is the apex of the peak.
	 * @return returns the m/z of the peak.
	 */
	private double quadraticFit(List<Peak> pks, int index){
		double x1, x2, x3;
		double y1, y2, y3;
		double d;
		
		if (index <1)
			return pks.get(0).mz();
		if (index >= (int)pks.size() - 1)
			return pks.get(pks.size()-1).mz();

		x1 = pks.get(index - 1).mz();
		x2 = pks.get(index).mz();
		x3 = pks.get(index + 1).mz();
		y1 = pks.get(index - 1).in();
		y2 = pks.get(index).in();
		y3 = pks.get(index + 1).in();

		d = (y2 - y1)*(x3 - x2) - (y3 - y2)*(x2 - x1);
		if (d == 0)
			return x2;  // no good.  Just return the known peak
		d = ((x1 + x2) - ((y2 - y1) * (x3 - x2) * (x1 - x3)) / d) / 2.0;
		return d;	// Calculated new peak.  Return it.
	}
	
	/**
	 *  Gets the peak that fits the point at a given index with a Lorentzian fit.
	 * @param mzs List of raw data of m\zs.
	 * @param intensities List of raw data of intensities.
	 * @param index index of the point in the m/z lists which is the apex of the peak.
	 * @param FWHM
	 * @return
	 */
	private double lorentzianFit(List<Peak> pks, int index, double FWHM){
		double A = pks.get(index).in();
		double Vo = pks.get(index).mz();
		int lstart, lstop;

		double E, CE, le;

		E = (Vo - pks.get(index + 1).mz()) / 100;
		E = Constants.absolute(E);

		if (index <1)
			return pks.get(index).mz();
		if (index == pks.size())
			return pks.get(index).mz();

		lstart = PeakIndex.getNearest(pks, Vo + FWHM, index) + 1;
		lstop = PeakIndex.getNearest(pks, Vo - FWHM, index) - 1;

		CE = lorentzianLS(pks, A, FWHM, Vo, lstart, lstop);
		for (int i = 0; i < 50; i++)
		{
			le = CE;
			Vo = Vo + E;
			CE = lorentzianLS(pks, A, FWHM, Vo, lstart, lstop);
			if (CE > le)
				break;
		}

		Vo = Vo - E;
		CE = lorentzianLS(pks, A, FWHM, Vo, lstart, lstop);
		for (int i = 0; i < 50; i++)
		{
			le = CE;
			Vo = Vo - E;
			CE = lorentzianLS(pks, A, FWHM, Vo, lstart, lstop);
			if (CE > le)
				break;
		}
		Vo = Vo + E;
		return Vo;
	}
	
	/**
	 *  Gets the peak that fits the point at a given index with a Lorentzian least square fit.
	 * @param mzs List of raw data of m\zs.
	 * @param intensities List of raw data of intensities.
	 * @param A
	 * @param FWHM
	 * @param Vo
	 * @param lstart
	 * @param lstop
	 * @return returns the m/z of the fit peak.
	 */
	private double lorentzianLS(List<Peak> pks, double A, double FWHM, double Vo, int lstart, int lstop){
		double u;
		double Y1, Y2;
		double RMSerror = 0;

		for (int index = lstart; index <= lstop; index++)
		{
			u = 2 / FWHM * (pks.get(index).mz() - Vo);
			Y1 = (int)(A / (1 + u * u));
			Y2 = pks.get(index).in();
			RMSerror = RMSerror + (Y1 - Y2)*(Y1 - Y2);
		}
		return RMSerror;
	}

	/**
	 *  Sets the type of fit.
	 * @param type sets the type of fit function that this instance uses.
	 */
	public void setOptions(PEAK_FIT_TYPE type){
		menm_peak_fit_type = type;
	}
	//! Gets the peak that fits the point at a given index by the specified peak fit function.
	/*!
	\param index index of the point in the m/z vectors which is the apex of the peak.
	\param mzs std::vector of raw data of m\zs.
	\param intensities std::vector of raw data of intensities.
	\return returns the m/z of the peak.
	*/
	public double fit(int index, List<Peak> pks){
		if (menm_peak_fit_type == PEAK_FIT_TYPE.APEX)
			return pks.get(index).mz();
		else if (menm_peak_fit_type == PEAK_FIT_TYPE.QUADRATIC)
			return this.quadraticFit(pks, index);
		else if (menm_peak_fit_type == PEAK_FIT_TYPE.LORENTZIAN)
		{
			double FWHM = this.mobj_peak_statistician.findFWHM(pks, index, 0.0);
			if (FWHM != 0)
				return this.lorentzianFit(pks, index, FWHM);
			return pks.get(index).mz();
		}
		return 0.0;
	}
}
