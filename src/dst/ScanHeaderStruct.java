package dst;

public class ScanHeaderStruct {
	public int  msLevel;
	public int  peaksCount;

	/**********************hanynag ISA************************/
	public int centroided;
	/**********************hanynag ISA************************/

	public double retentionTime;        /* in seconds */
	public double lowMZ;
	public double highMZ;
	public double precursorMZ;  /* only if MS level > 1 */
	/********************* slee ****************************/
	public int precursorCharge;
	/********************* slee ****************************/
	public double basePeakMz;
	public double basePeakIntensity;
	public double totIonCurrent;
}
