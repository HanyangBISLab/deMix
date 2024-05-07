package ppk;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;
import utils.Constants;

public class PeakStatistician {
	public PeakStatistician() {
		mvect_mzs = new ArrayList<>();
		mvect_intensities = new ArrayList<>();
	}
	
	/**
	 * internal variable to store temporary m/z values. It also guarantees a workspace that doesn't need to be reallocated all the time.
	 */
	private List<Double> mvect_mzs;
	
	/**
	 *  internal variable to store temporary intensity values. It also guarantees a workspace that doesn't need to be reallocated all the time.
	 */
	private List<Double> mvect_intensities;
		
	/** Find signal to noise value at position specified. 
	 * remarks Looks for local minima on the left and the right hand sides and calculates signal to noise of peak relative to the minimum of these shoulders.
	 * @param yValue is intenstiy at specified index.
	 * @param vect_intensities is list of intensities.
	 * @param index is position of point at which we want to calculate signal to noise.
	 * @return returns computed signal to nosie value.
	 */
	public double findSignalToNoise(double yValue, List<Peak> vect_pks, int index) { // The place in arrDerivative the derivative crossed 0
		double minIntensityLeft = 0, minIntensityRight = 0;
		int num_data_pts = vect_pks.size();
		if (yValue == 0)
			return 0;

		if (index <= 0 || index >= num_data_pts - 1)
			return 0;

		// Find the first local minimum as we go down the m/z range.
		boolean found = false;
		for (int i = index; i > 0; i--) {
			if (vect_pks.get(i + 1).in() >= vect_pks.get(i).in()
					&& vect_pks.get(i - 1).in() > vect_pks.get(i).in()) // Local minima here 
			{
				minIntensityLeft = vect_pks.get(i).in();
				found = true;
				break;
			}
		}
		if (!found)
			minIntensityLeft = vect_pks.get(0).in();

		found = false;
		//// Find the first local minimum as we go up the m/z range.
		for (int i = index; i < num_data_pts - 1; i++) {
			if (vect_pks.get(i + 1).in() >= vect_pks.get(i).in()
					&& vect_pks.get(i - 1).in() > vect_pks.get(i).in()) // Local minima here 
			{
				minIntensityRight = vect_pks.get(i).in();
				found = true;
				break;
			}
		}
		if (!found)
			minIntensityRight = vect_pks.get(num_data_pts - 1).in();
		if (minIntensityLeft == 0) {
			if (minIntensityRight == 0)
				return 100;
			return (1.0 * yValue) / minIntensityRight;
		}
		if (minIntensityRight < minIntensityLeft && minIntensityRight != 0)
			return (1.0 * yValue) / minIntensityRight;

		return (1.0 * yValue) / minIntensityLeft;

	}

//	private double mse;
	
	/** 
	 * coe is coefficients found by curve regression.
	 */
	private double coe[] = new double[2];
	
	/** Find full width at half maximum value at position specified. 
	 * remarks Looks for half height locations at left and right side, and uses twice of that value as the FWHM value. If half height
			locations cannot be found (because of say an overlapping neighbouring peak), we perform interpolations.
	 * @param vect_mzs list of mzs.
	 * @param vect_intensities list of intensities.
	 * @param data_index is position of point at which we want to calculate FWHM.
	 * @return returns computed FWHM.
	 */
	double findFWHM(List<Peak> vect_pks, int data_index, double signalToNoise) {
		double upper, lower, mass;
		double X1, X2;
		double Y1, Y2, peak, peakHalf;
		int iStat;
		int points;
		

		points = 0;
		peak = vect_pks.get(data_index).in();

		peakHalf = ((double) peak / 2.0);

		mass = vect_pks.get(data_index).mz();

		if (peak == 0.0)
			return 0.0;

		int num_input_pts = (int) vect_pks.size();

		if (data_index <= 0 || data_index >= num_input_pts - 1)
			return 0;

		upper = vect_pks.get(0).mz();

		for (int index = data_index; index >= 0; index--) {
			double current_mass = vect_pks.get(index).mz();
			Y1 = vect_pks.get(index).in();
			if ((Y1 < peakHalf) || (Constants.absolute(mass - current_mass) > 5.0)
					|| ((index < 1 || vect_pks.get(index - 1).in() > Y1)
							&& (index < 2 || vect_pks.get(index - 2).in() > Y1) && (signalToNoise < 4.0))) {
				Y2 = vect_pks.get(index + 1).in();
				X1 = vect_pks.get(index).mz();
				X2 = vect_pks.get(index + 1).mz();
				if ((Y2 - Y1 != 0) && (Y1 < peakHalf)) {
					upper = X1 - (X1 - X2) * (peakHalf - Y1) / (Y2 - Y1);
				} else {
					upper = X1;
					points = data_index - index + 1;
					if (points >= 3) {
						mvect_mzs.clear();
						mvect_intensities.clear();
						// if (mvect_mzs.capacity() < points)
						// {
						// mvect_mzs.reserve(points);
						// mvect_intensities.reserve(points);
						// }

						int j = points - 1;
						for (; j >= 0; j--) {
							mvect_mzs.add(vect_pks.get(data_index - j).mz());
							mvect_intensities.add(vect_pks.get(data_index - j).in());
						}
						for (j = 0; j < (int) points && mvect_intensities.get(0) == mvect_intensities.get(j); j++)
							;
						if (j == points)
							return 0.0;
						iStat = curvReg(mvect_intensities, mvect_mzs, points, 1);
						// only if successful calculation of peak was done,
						// should we change upper.
						if (iStat != -1)
							upper = coe[1] * peakHalf + coe[0];
					}
				}
				break;

			}
		}

		lower = vect_pks.get(num_input_pts - 1).mz();
		for (int index = data_index; index < num_input_pts; index++) {
			double current_mass = vect_pks.get(index).mz();
			Y1 = vect_pks.get(index).in();
			if ((Y1 < peakHalf) || (Constants.absolute(mass - current_mass) > 5.0)
					|| ((index > num_input_pts - 2 || vect_pks.get(index + 1).in() > Y1)
							&& (index > num_input_pts - 3 || vect_pks.get(index + 2).in() > Y1)
							&& signalToNoise < 4.0)) {
				Y2 = vect_pks.get(index - 1).in();
				X1 = vect_pks.get(index).mz();
				X2 = vect_pks.get(index - 1).mz();

				if ((Y2 - Y1 != 0) && (Y1 < peakHalf)) {
					lower = X1 - (X1 - X2) * (peakHalf - Y1) / (Y2 - Y1);

				} else {
					lower = X1;
					points = index - data_index + 1;
					if (points >= 3) {
						mvect_mzs.clear();
						mvect_intensities.clear();
						// if (mvect_mzs.capacity() < points)
						// {
						// mvect_mzs.reserve(points);
						// mvect_intensities.reserve(points);
						// }
						for (int k = points - 1; k >= 0; k--) {
							mvect_mzs.add(vect_pks.get(index - k).mz());
							mvect_intensities.add(vect_pks.get(index - k).in());
						}
						int j = 0;
						for (j = 0; j < (int) points && mvect_intensities.get(0) == mvect_intensities.get(j); j++)
							;
						if (j == points)
							return 0.0;

						iStat = curvReg(mvect_intensities, mvect_mzs, points, 1);
						// only if successful calculation of peak was done,
						// should we change lower.
						if (iStat != -1)
							lower = coe[1] * peakHalf + coe[0];
					}
				}
				break;
			}
		}

		if (upper == 0.0)
			return 2 * Constants.absolute(mass - lower);
		if (lower == 0.0)
			return 2 * Constants.absolute(mass - upper);
		return Constants.absolute(upper - lower);
	}

	/**
	 *  Calculate Least Square error mapping y = f(x). 
	 * @param x list of x values.
	 * @param y list of y values.
	 * @param n number of points in list.
	 * @param nterms order of the function y = f(x).
	 * @return returns 0 if successful an -1 if not.
	 */
	private int curvReg(List<Double> x, List<Double> y, int n, int nterms){
		Matrix I_At_At_T, At_At_T, At_T;
		int i, j;
		double w[];
		Matrix at, b, z, out;
		double yfit, xpow;

		// weights
		w = new double[n];
//		if (!w) return(-1);
		for (i = 0; i < n; i++) w[i] = 1.0;

		// weighted powers of x matrix transpose = At 
		at = new Matrix(nterms + 1, n);
		for (i = 0; i < n; i++) {
			at.set(0, i, w[i]);
			for (j = 1; j < (nterms + 1); j++)
				at.set(j, i, at.get(j-1, i)*x.get(i));
		}
		
		// Z = weighted y std::vector 
		z = new Matrix(n, 1);
		
		for (i = 0; i < n; i++)
			z.set(i, 0, w[i] * y.get(i));
		
		At_T = at.transpose();
		At_At_T = at.times(At_T);
		I_At_At_T = At_At_T.inverse();
		if (I_At_At_T == null)
		{
			return -1;
		}

		Matrix At_Ai_At = I_At_At_T.times(at);

		b =  At_Ai_At.times(z);

		// make a matrix with the fit y and the difference 
		out = new Matrix(2, n);

		// calculate the least squares y values 
//		double mse = 0.0;
		for (i = 0; i < n; i++) {
			coe[0] = b.get(0, 0);
			yfit = b.get(0, 0);
			xpow = x.get(i);
			for (j = 1; j <= nterms; j++) {
				coe[j] = b.get(j, 0);
				yfit += b.get(j, 0) * xpow;
				xpow = xpow * x.get(i);
			}
			out.set(0, i, yfit);
			out.set(1, i, y.get(i) - yfit);
//			mse += y.get(i) - yfit;
		}
		
		return 0;
	}

}
