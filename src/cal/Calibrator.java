package cal;

import utils.Constants.CALIBRATION_TYPE;

public class Calibrator {
	public CalibratorImp mobj_calibrator_imp;
	public CALIBRATION_TYPE menm_calibration_type;
	
	public Calibrator()
	{
		mobj_calibrator_imp = null;
		menm_calibration_type = CALIBRATION_TYPE.UNDEFINED;
	}
	public Calibrator(CALIBRATION_TYPE cal_type)
	{
		initialize(cal_type);
	}

	public void initialize(CALIBRATION_TYPE cal_type)
	{
		menm_calibration_type = cal_type;
		mobj_calibrator_imp = new CalibratorImp(cal_type);
	}
}
