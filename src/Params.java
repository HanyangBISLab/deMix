package kr.ac.hanyang.bislab.demix;

public class Params {

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public static final double	Electron = 0.000549;
	public static final double	Hydrogen = 1.007825035;
	public static final double	Oxygen = 15.99491463;
	public static final double	Nitrogen = 14.003074;
	public static final double	Proton = Hydrogen-Electron;
	public static final double	HO = Hydrogen + Oxygen;	
	public static final double	H2O = Hydrogen*2 + Oxygen;	
	public static final double	NH3 = Hydrogen*3 + Nitrogen;		
//	public static final double	yionOffset = H2O+Proton;
	public static final double	A_ION_OFFSET = Oxygen + 12.;
	
	public static final double	IsotopeSpace = 1.00235;
	public static final double  DeutriumIsotopeSpace= 1.005;	//1.00782465934753;

	
	public static final double  IMM_OFFSET = -A_ION_OFFSET + Proton;////////////////////////////////////////////////
	
	
	public static int minObservedScans= 3;
	public static int maxElutionHoleScans= 2;

	public static int minObservedHDXPeaks= 5;

	public static int extendibleScanRange= 40;
	
	public static double minPeaksIntensiy= 25;
	
	
	
	
	
}
