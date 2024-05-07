package rds;

import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.List;

import org.systemsbiology.jrap.stax.MSXMLParser;
import org.systemsbiology.jrap.stax.Scan;
import org.systemsbiology.jrap.stax.ScanHeader;

import dst.InstrumentStruct;
import dst.ScanHeaderStruct;
import ppk.Peak;
import utils.Constants.FILE_TYPE;

public class MZXmlRawData extends RawData{
	public MZXmlRawData() {
		mint_num_scans = 0;
		mint_current_scan = -1;
		parser = null;
		currrnt_scan = null;
	}
	
	String mstr_file_name;

	private int mint_num_scans;
	int mint_current_scan;
	Scan currrnt_scan;
	ScanHeader currrnt_scan_header;
	double mdbl_signal_level;
	MSXMLParser parser;
	
	@Override
	public String getFileName() {
		return mstr_file_name;
	}
	@Override
	public FILE_TYPE getFileType() {
		return menm_file_type;
	}
	@Override
	public boolean getRawData(List<Peak> pks, int scan_num) {
		return getRawData(pks, scan_num, -1);
	}
	@Override
	public boolean getRawData(List<Peak> pks, int scan_num, int num_pts) {
	/*	if (scan_num > mint_num_scans)
		{
			String mesg = "";
			mesg +=  "File only has " + mint_num_scans + " scans. Cannot read to scan number: " + scan_num;
			System.exit(-1);
		}

		mint_current_scan = scan_num;
		currrnt_scan = parser.rap(mint_current_scan);
		currrnt_scan_header = parser.rapHeader(scan_num);//*/

		pks.clear();

		int num_points = currrnt_scan_header.getPeaksCount();
		if (num_pts > 0 && num_pts < num_points)
		{
			num_points = num_pts;
		}
		
		List<Peak> peaks = readPeaks();
		Peak fpeak;
		double max_intensity = -1 * Double.MAX_VALUE, min_intensity = Double.MAX_VALUE;
		for (int pt_num = 0; pt_num < num_points; pt_num ++)
		{
			fpeak = peaks.get(pt_num);
			if (fpeak.mdbl_intensity > max_intensity)
				max_intensity = fpeak.mdbl_intensity;
			if (fpeak.mdbl_intensity > min_intensity)
				min_intensity = fpeak.mdbl_intensity;
			pks.add(fpeak);
		}

		mdbl_signal_level = max_intensity - min_intensity;

		return true;
	}
	
	private List<Peak> readPeaks() {
		// TODO Auto-generated method stub
		double[][] peaklist = currrnt_scan.getMassIntensityList();
		List<Peak> vec_pk = new ArrayList<>();
		for( int i = 0; i < peaklist[0].length; i++ ){
			Peak pk = new Peak();
			pk.mdbl_mz = peaklist[0][i];
			pk.mdbl_intensity = peaklist[1][i];
			vec_pk.add(pk);
		}
		return vec_pk;
	}
	@Override
	public boolean load(String file_name) {
		mstr_file_name = file_name;

		/**********************For JAVA************************/
		parser = new MSXMLParser(mstr_file_name);
		mint_num_scans = parser.getMaxScanNumber();//parser.getScanCount();[by na]
		/**********************For JAVA************************/
		
	/*	mint_current_scan = 1;
		currrnt_scan = parser.rap(mint_current_scan);
		currrnt_scan_header = currrnt_scan.getHeader();//*/
		
		for(int i=1; i<=mint_num_scans; i++) {//[by na]
			currrnt_scan = parser.rap(i);
			if( currrnt_scan == null ) continue;
			
			currrnt_scan_header = currrnt_scan.getHeader();
			mint_current_scan = i;
			break;
		}
		
//		System.out.println("CollisionEnergy() : " + currrnt_scan_header.getCollisionEnergy());
//		System.out.println("DoubleRetentionTime() : " + currrnt_scan_header.getDoubleRetentionTime());
//		System.out.println("MsLevel() : " + currrnt_scan_header.getMsLevel());
//		System.out.println("RetentionTime() : " + currrnt_scan_header.getRetentionTime());
//		System.out.println("RT() : " + currrnt_scan_header.getRT());
//		System.out.println("TotIonCurrent() : " + currrnt_scan_header.getTotIonCurrent());
		
		/**********************hanynag ISA************************/
		return parser.isMzXML(mstr_file_name);
		/**********************hanynag ISA************************/
	}
	
	@Override
	public int getNumScans() {
		return mint_num_scans;
	}
	
	@Override
	public double getScanTime(int scan_num) {
		if (mint_current_scan == scan_num)
		{
			return currrnt_scan_header.getDoubleRetentionTime();
		}
		else
		{
			return parser.rapHeader(scan_num).getDoubleRetentionTime();
		}
	}
	@Override
	public boolean isMSScan(int scan_num) {
		int ms_level = getMSLevel(scan_num);
		if (ms_level == 1)
			return true;
		else
			return false;
	}
	@Override
	public int getMSLevel(int scan_num) {
	//	System.out.println(mint_current_scan +"      "+ scan_num) ;
		if (mint_current_scan == scan_num)
		{
			return currrnt_scan_header.getMsLevel();
		}
		else
		{
			return parser.rapHeader(scan_num).getMsLevel();
		}
	}
	@Override
	public boolean isProfileScan(int scan_num) {
		if (mint_current_scan == scan_num)
		{
			return !(currrnt_scan_header.getCentroided() > 0);
		}
		else
		{
			return !(parser.rapHeader(scan_num).getCentroided() > 0);
		}
		
	}
	
	public void setParentMZ(int scan_num, double precursorMz) {
		if (mint_current_scan == scan_num)
		{
			currrnt_scan_header.setPrecursorMz((float) precursorMz);
		}
		else
		{
			ScanHeader temp_header = parser.rapHeader(scan_num);
			temp_header.setPrecursorMz((float) precursorMz);
		}
	}
	
	@Override
	public double getParentMz(int scan_num) {
		if (mint_current_scan == scan_num)
		{
			return currrnt_scan_header.getPrecursorMz();
		}
		else
		{
			ScanHeader temp_header = parser.rapHeader(scan_num); 
			return temp_header.getPrecursorMz();
		}
	}
	@Override
	public int getParentCharge(int scan_num) {
		if (mint_current_scan == scan_num)
		{
			return currrnt_scan_header.getPrecursorCharge();
		}
		else
		{
			return parser.rapHeader(scan_num).getPrecursorCharge();
		}
	}
	@Override
	public double getMonoMZFromHeader(int scan_num) {
		if (mint_current_scan == scan_num)
		{
			return currrnt_scan_header.getPrecursorMz();
		}
		else
		{
			return parser.rapHeader(scan_num).getPrecursorMz();
		}
	}
	
	@Override
	public boolean readyScan(int scan_num) {
		if (scan_num == mint_current_scan) return true;
		
		currrnt_scan = parser.rap(scan_num);
		if( currrnt_scan == null ) {
			if (scan_num > mint_num_scans) {
				String mesg = "";
				mesg +=  "File only has " + mint_num_scans + " scans. Cannot read to scan number: " + scan_num;
				System.exit(-1);
			}
			currrnt_scan = null;
			mint_current_scan = -1;
			return false;
		}
		else {
			currrnt_scan_header = currrnt_scan.getHeader();
			mint_current_scan = scan_num;
			return true;
		}
	}
	
}
