package kr.ac.hanyang.bislab.demix.mzfile;

import java.io.IOException;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;


public abstract class MZFormatParser {
	
	public static int MIN_ASSUMED_CHARGE = 2, MAX_ASSUMED_CHARGE = 4;
	
	public int maxScanNumber;
	public TreeMap<Integer, Double> retentionTimeMap;

	public String fileName;
	public String baseName;

	
	public MZFormatParser(String fName) throws IOException {
		fileName = fName;
		baseName = fileName.replace('\\', '/');
		baseName = baseName.substring(baseName.lastIndexOf('/')+1, baseName.lastIndexOf('.'));
		
		maxScanNumber = 0;
		retentionTimeMap = new TreeMap<Integer, Double>();
	}
	
	
	public String getFileName(){ return fileName; }
	public int getMaxScanNumber(){ return maxScanNumber; }
	
	public Set<Integer> getScanSet(){ return retentionTimeMap.keySet(); }
	
	public TreeMap<Double, Integer> getRt2ScanMap(){ 
		
		TreeMap<Double, Integer> rt2scan= new TreeMap<Double, Integer>();
		for(Map.Entry<Integer, Double> entry : retentionTimeMap.entrySet()) {
			rt2scan.put(entry.getValue(), entry.getKey());
		}
		return rt2scan; 
		
		
	}
	
	
	
	
	public abstract double[][] getPeakList(int sn) throws Exception; 	
	public abstract double getRetentionTime(int sn) throws Exception;
	
	
	public static MZFormatParser get(String specFile) throws Exception {
	//	System.out.print( "Reading MS/MS spectra.....  " );
		String ftype= specFile.toLowerCase();
		
		MZFormatParser parser = null;
		if( ftype.endsWith(".mzxml") ){
			parser = new mzXMLFormatParser( specFile );
		}
		else if( ftype.endsWith(".mzml") || ftype.endsWith(".gz") ){
		//	parser = new mzMLFormatParser( specFile );
			parser = new DPmzML( specFile );
		}	
		
		return parser;
	}
}
