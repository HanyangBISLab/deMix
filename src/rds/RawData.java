package rds;

import java.util.ArrayList;
import java.util.List;

import cal.Calibrator;
import ppk.Peak;
import ppk.PeakIndex;
import utils.Constants.FILE_TYPE;

public abstract class RawData {
	public RawData() {
		mobj_calibrator = null;
	}
	public FILE_TYPE menm_file_type;
	public Calibrator mobj_calibrator;
	//		public void GetSummedSpectra(std::vector <double> *bins, std::vector <double> *intensities, int scan, int scan_range, 
	//			double min_mz, double max_mz, double mz_bin) ;
	public void getRawData(List<Peak> vectpks, int scan, double min_mz,	double max_mz) {
		if (max_mz <= min_mz){
			System.out.println("max_mz should be greater than min_mz");
			return;
		}
		vectpks.clear();;
		List<Peak> allpks = new ArrayList<>();
		getRawData(allpks, scan);
		int numPts = allpks.size();
		if (numPts == 0)
			return;
		int startIndex = PeakIndex.getNearestBinary(allpks, min_mz, 0, numPts - 1);
		int stopIndex = PeakIndex.getNearestBinary(allpks, max_mz, 0, numPts - 1);
		allpks.addAll(vectpks.subList(startIndex, stopIndex));
	}
	
	abstract public String getFileName();
	abstract public FILE_TYPE getFileType();
	abstract public boolean getRawData(List<Peak> pks, int scan_num);
	abstract public boolean getRawData(List<Peak> pks, int scan_num, int num_pts);

	/*******************hanyang ISA************************/
	abstract public boolean load(String file_n);
	abstract public boolean readyScan(int scan_num);
	/*******************hanyang ISA************************/

	abstract public int getNumScans();
	abstract public double getScanTime(int scan_num);
//	abstract void getScanDescription(int scan, String description)
//	{
//		strcpy(description, "Scan #");
//		_itoa(scan, &description[strlen(description)], 10);
//	}
	abstract public boolean isMSScan(int scan_num);
	abstract public int getMSLevel(int scan_num);
	abstract public boolean isProfileScan(int scan_num);
	abstract public double getParentMz(int scan_num);
	/*************************slee**************************/
	abstract int getParentCharge(int scan_num);
	/*************************slee**************************/
	abstract public double  getMonoMZFromHeader(int scan_num);
}
