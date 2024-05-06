package kr.ac.hanyang.bislab.demix.mzfile;

import java.io.IOException;

import org.systemsbiology.jrap.stax.MSXMLParser;
import org.systemsbiology.jrap.stax.Scan;
import org.systemsbiology.jrap.stax.ScanHeader;

public class mzXMLFormatParser extends MZFormatParser {

	private MSXMLParser parser;
		
	public mzXMLFormatParser( String fileName ) throws Exception {
		
		super(fileName);	
		try{
		parser = new MSXMLParser(fileName);	
		maxScanNumber = parser.getMaxScanNumber();
			
		for (int i=1; i<=maxScanNumber; i++){
			Scan scan = parser.rap(i);	
			if( scan == null ) continue;			
 			ScanHeader shead = scan.getHeader();
			if( shead.getMsLevel() > 1 ) continue; //(select MS spectrum only)						
			retentionTimeMap.put(i, shead.getDoubleRetentionTime());
		}

		
		} catch (Exception e) {
			System.out.println("Abnormal Termination");
			System.exit(1);
		}
	}
	
	public double[][] getPeakList(int sn) throws Exception {		
		if( retentionTimeMap.get(sn) == null ) return null;
		return parser.rap(sn).getMassIntensityList();
	}


	public double getRetentionTime(int sn) throws IOException {
		return retentionTimeMap.getOrDefault(sn, 0.);
	}


}


