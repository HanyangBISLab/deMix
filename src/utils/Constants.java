package utils;

import java.util.List;
import java.util.Map;

import ppk.PeakIndex;

public class Constants {
	public Constants() {
		new PeakIndex();
	}
	
	public static int NEGINF = -999999999;
	public static int POSINF = 999999999;
	
	public static int maxMissPeak = 3;				// Distribution 앞부분에 사라지는 peak 개수
	public static int maxCharge = 10;				// 찾아볼 charge state 최대값
	public static int scanNoModifier = 0;			// 출력되는 scan number에 더해줄 값
	public static int maxAbundancePeak = 3;			// 출력되는 abundance에 더해질 peak의 개수
	public static double massErr = 1.0E-05;			// 질량 오차 허용범위
	public static double thMassAbn1 = 300;			// 첫번째 peak가 가장 큰 질량 범위
	public static double thMassAbn2 = 1800;			// 두번째 peak가 가장 큰 질량 범위
	public static double thClustExt = -2.0;			// extendCluster에서 사용하는 threshold
	public static double thScore = -0.1;			// 답으로 리포트하는 threshold			
													// 0.0 -> -0.1 by slee
	public static double wgtR1 = 3.0;				// I2/I1에 가중치
	public static boolean truncated = false;
	public static double backgroundIntensity = 0;
	public static double thMonoMassRange = 0.03;	//monoisotoped mass 가 일정범위 이하일때 같은걸로 취급	
	
	public static double MASS_ONE		=	1.00235;
	public static double MASS_H			=	1.0078246;
	public static double MASS_ELECTRON	=	0.00054858;
	
	public static double T_MASS1 = 400;
	public static double T_MASS2 = 1800;
	public static double T_MASS3 = 4000;

	public static int BACKGROUND_RATIO_FOR_TIC = 3;
	public static double MIN_MZ = 400;
	public static double MAX_MZ = 2000;
	
	public static Map<String, String> _tbl;
	public static Map<String, List<String>> _tbl2;
	
	public static String ResultOrder;
	public static String OutputFormat;
	public static String Target;
	public static String MSMS;

	public static int DataFileType = 1; // default : Finnigan Raw
	public static int Partition = 5;
	public static int MakeMSMS = 0;
	public static boolean PrintPeak = true;
	
	public static FILE_TYPE getFileType(int file_Type){
		switch (file_Type) {
		case 0:
			return FILE_TYPE.BRUKER;
		case 1:
			return FILE_TYPE.FINNIGAN;
		case 2:
			return FILE_TYPE.MICROMASSRAWDATA;
		case 3:
			return FILE_TYPE.AGILENT_TOF;
		case 4:
			return FILE_TYPE.SUNEXTREL;
		case 5:
			return FILE_TYPE.ICR2LSRAWDATA;
		case 6:
			return FILE_TYPE.MZXMLRAWDATA;
		case 7:
			return FILE_TYPE.PNNL_IMS;
		case 8:
			return FILE_TYPE.BRUKER_ASCII;
		case 9:
			return FILE_TYPE.ASCII;
		default:
			break;
		}
		return FILE_TYPE.FINNIGAN;
	}
	
//	public static String getRound(double d, int pos) {
//		if( d == 0 ) return "0.000000"; 
//		String out = ""+((double)Math.round(d*(Math.pow(10, pos))))/Math.pow(10, pos);
//		int o_len = out.length()-out.lastIndexOf('.'); 
//		if( o_len < 7 ) {
//			for( int i = 0; i < 7-o_len; i++ ) {
//				out += "0";
//			}
//		}
//		return out;
//	}
	
	public static int getMaxMissPeak() {
		return maxMissPeak;
	}

	public static void setMaxMissPeak(int maxMissPeak) {
		Constants.maxMissPeak = maxMissPeak;
	}

	public static int getScanNoModifier() {
		return scanNoModifier;
	}

	public static void setScanNoModifier(int scanNoModifier) {
		Constants.scanNoModifier = scanNoModifier;
	}

	public static int getMaxAbundancePeak() {
		return maxAbundancePeak;
	}

	public static void setMaxAbundancePeak(int maxAbundancePeak) {
		Constants.maxAbundancePeak = maxAbundancePeak;
	}

	public static double getMassErr() {
		return massErr;
	}

	public static void setMassErr(double massErr) {
		Constants.massErr = massErr;
	}

	public static double getThClustExt() {
		return thClustExt;
	}

	public static void setThClustExt(double thClustExt) {
		Constants.thClustExt = thClustExt;
	}

	public static double getThScore() {
		return thScore;
	}

	public static void setThScore(double thScore) {
		Constants.thScore = thScore;
	}

	public static int getMaxCharge() {
		return maxCharge;
	}

	public static void setMaxCharge(int maxCharge) {
		Constants.maxCharge = maxCharge;
	}

	
		//! enumeration for type of fit. 
	public static enum PEAK_FIT_TYPE
	{
		APEX,  /*!< The peak is the m/z value higher than the points before and after it */
		QUADRATIC, /*!< The peak is the m/z value which is a quadratic fit of the three points around the apex */
		LORENTZIAN /*!< The peak is the m/z value which is a lorentzian fit of the three points around the apex */
	};
	public static enum PEAK_PROFILE_TYPE
	{
		PROFILE, CENTROIDED
	};
	
	public static enum FILE_TYPE
	{ 
		BRUKER,	FINNIGAN, MICROMASSRAWDATA, AGILENT_TOF, SUNEXTREL, ICR2LSRAWDATA, MZXMLRAWDATA, PNNL_IMS, BRUKER_ASCII, ASCII 
	};

	public static enum CALIBRATION_TYPE
	{
		A_OVER_F_PLUS_B_OVER_FSQ_PLUS_C_OVERFCUBE, A_OVER_F_PLUS_B_OVER_FSQ,
		A_OVER_F_PLUS_B_OVER_FSQ_PLUS_CI_OVERFSQ, A_OVER_F_PLUS_B_OVER_FSQ_PLUS_C,
		AF_PLUS_B, F, A_OVER_F_PLUS_B, A_OVER_F_PLUS_B_PLUS_CI,
		TIME_A_TSQ_PLUS_B_T_PLUS_C, BRUKER_CALMET, UNDEFINED
	};
	
	/**
	 * 
	 * @param x
	 * @return returns positive value 
	 */
	public static double absolute(double x)
	{
		if (x >= 0)
			return x;
		return -1*x;
	}
}
