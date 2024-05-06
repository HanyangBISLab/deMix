package kr.ac.hanyang.bislab.demix.mzfile;


import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.digitalproteomics.oss.parsers.mzml.MzMLStAXParser;
import com.digitalproteomics.oss.parsers.mzml.builders.XMLSpectrumBuilder;
import com.digitalproteomics.oss.parsers.mzml.builders.XMLSpectrumHeaderBuilder;
import com.digitalproteomics.oss.parsers.mzml.model.Spectrum;
import com.digitalproteomics.oss.parsers.mzml.model.SpectrumHeader;

public class DPmzML extends MZFormatParser {
	
	private  MzMLStAXParser<Spectrum> parser;
	
	public DPmzML( String fileName ) throws IOException {
		
		super(fileName);		
		try{
		
		MzMLStAXParser<SpectrumHeader> headparser= new MzMLStAXParser<SpectrumHeader>(Paths.get(fileName), XMLSpectrumHeaderBuilder::new);	
		for( SpectrumHeader shead : headparser ) {
			if( shead.getMsLevel() > 1 ) continue;
			
			int scanNumber = convertScanIdToScanNumber(shead.getId());
			
			double rt= shead.getScanStartTime();
			retentionTimeMap.put(scanNumber, rt);
			
			if( maxScanNumber < scanNumber )
				maxScanNumber= scanNumber;
 		}
		
		parser= new MzMLStAXParser<Spectrum>(Paths.get(fileName), XMLSpectrumBuilder::new);	
		
		} catch (Exception e) {
			System.out.println("Abnormal Termination");
			System.exit(1);
		}
		
	}

	@Override
	public double[][] getPeakList(int sn) throws Exception {

		if( retentionTimeMap.get(sn) == null ) return null;
		Spectrum spectrum = parser.getSpectrumById("scan="+sn);
		
		List<Double> mzs= spectrum.getMz();
		List<Double> intensities= spectrum.getIntensities();	
	    if( mzs == null || intensities == null ) return null;

	    double[][] peakList= new double[2][mzs.size()];
	    for (int i = 0; i<mzs.size(); i++) {
	    	peakList[0][i]= mzs.get(i);
	    	peakList[1][i]= intensities.get(i);
	    }
	    
	    return peakList;
	}

	
	public double getRetentionTime(int sn) throws Exception {
		return retentionTimeMap.getOrDefault(sn, 0.);
	}
	
	private int convertScanIdToScanNumber(String scanId) {

		final Pattern pattern = Pattern.compile("scan=([0-9]+)");
		final Matcher matcher = pattern.matcher(scanId);
		boolean scanNumberFound = matcher.find();
	
		// Some vendors include scan=XX in the ID, some don't, such as
		// mzML converted from WIFF files. See the definition of nativeID in
		// http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo
	    if( scanNumberFound )
	    	return Integer.parseInt(matcher.group(1));

	    return 0;
	}
	
	
	
}














