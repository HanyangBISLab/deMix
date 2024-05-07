package cal;

import utils.Constants.CALIBRATION_TYPE;;

public class CalibratorImp {
	// Type 0: m/z = A/f + B/f^2 + C/f^3
	// Type 1: m/z = A/f + |Vt|B/f^2
	// Type 2: m/z = A/f + |Vt|B/f^2 + I|Vt|C/f^2
	// Type 3: m/z = A/f + |Vt|B/f^2 + C
	// Type 4: m/z = f 
	// Type 5: m/z = A/(f+B)
	// Type 6: m/z = A/(f+B+CI)
	// Type 7: t = A*t^2 + B*t + C
	// Type 9: This is to support bruker calmet 1
	//         m/z = (-A - SQRT(A^2 - 4(B-f)C))/2(B-f)
	//         f   = A/mz + C/mz^2 +B
	
	
	final double MAX_MZ = 100000;
	final double MIN_MZ = 392.2944;
	final double MAX_MASS = 10000000;
	
	CALIBRATION_TYPE menm_calibration_type;
	// Normalizer for calibrator type 2.
	double mdbl_intensity_normalizer;
	
	public  int mint_num_points_in_scan;
	public  double mdbl_sample_rate;

	public  double mdbl_low_mass_frequency;
	public  double mdbl_frequency_shift; // don't know what it is yet exatcly. initialize to 0.

	public  double mdbl_calib_const_a;
	public  double mdbl_calib_const_b;
	public  double mdbl_calib_const_c;

	public int mint_byte_order;
	
	public CalibratorImp(CALIBRATION_TYPE type) {
		mdbl_frequency_shift = 0;
		mdbl_low_mass_frequency = 0;
		mdbl_calib_const_a = mdbl_calib_const_b = mdbl_calib_const_c = 0;
		mdbl_intensity_normalizer = 1;
		menm_calibration_type = type;
	}
	
	public CALIBRATION_TYPE getCalibrationType()
	{
		return menm_calibration_type;
	}
}
