package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.systemsbiology.jrap.stax.ScanAndHeaderParser;

import dst.InstrumentStruct;
import dst.RunHeaderStruct;
import dst.ScanHeaderStruct;

public class Ramp {
	static final boolean PROFILE_MODE = false;
	static final boolean CENTROIDED_MODE  = true;
	
	static final int INSTRUMENT_LENGTH = 2000;
	
	static final int SIZE_BUF = 512;
	
	/****************************************************************
	************** 		Find the Offset of the index	*************
	***************************************************************/
	public int getIndexOffset(File pFI) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(pFI));
		String line, target = "<indexOffset>";
		String sub = "-1";
		while( (line=br.readLine()) != null ){
			if( line.contains(target) ){
				sub = line.substring(line.indexOf('>')+1, line.lastIndexOf('<'));
				break;
			}
		}
		br.close();
		
		return Integer.parseInt(sub);
	}

	/****************************************************************
	* Reads the Scan index in a list				*
	* Returns pScanIndex which becomes property of the caller	*
	* pScanIndex is -1 terminated					
	 * @throws IOException *
	***************************************************************/	
	int[] readIndex(File pFI, int indexOffset,	int iLastScan, int pScanIndex) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(pFI));
		br.skip(indexOffset);
		br.readLine();	//	skip <index name="scan">
		String line;
		int cnt_scan = 0;
		while( (line=br.readLine()) != null ){
			if( !line.contains("<offset ") ){
				break;
			}
			cnt_scan++;
		}
		br.close();
		
		int[] off_scans = new int[cnt_scan];
		
		br = new BufferedReader(new FileReader(pFI));
		br.skip(indexOffset);
		br.readLine();	//	skip <index name="scan">
		
		int i = 0;
		while( (line=br.readLine()) != null ){
			if( !line.contains("<offset ") ){
				break;
			}
			String sub = line.substring(line.indexOf('>')+1, line.lastIndexOf('<'));
			int off = Integer.parseInt(sub);
			off_scans[i++] = off;
		}
		br.close();
		
		return off_scans;
	}
	
	/*
	* Reads scan header information.
	* !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE
	*    RETURNING !
	*/
	public static void readHeader(File pFI, int lScanIndex, ScanHeaderStruct scanHeader, ScanAndHeaderParser shp) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(pFI));
		br.skip(lScanIndex);
		String line;
		
		/*
		* initialize defaults
		*/
		scanHeader.msLevel = 0;
		scanHeader.peaksCount = 0;

		/**********************hanynag ISA************************/
		scanHeader.centroided = 0;
		/**********************hanynag ISA************************/

		scanHeader.retentionTime = 0;
		scanHeader.lowMZ = 1.E6;
		scanHeader.highMZ = 0.0;
		scanHeader.precursorMZ = 0.0;
		scanHeader.basePeakIntensity = 0.0;
		scanHeader.basePeakMz = 0.0;
		scanHeader.totIonCurrent = 0.0;

		while ( (line=br.readLine()) != null )
		{
			if (line.contains("<peaks "))
				return;
			if (line.contains("msLevel=\"")){
				String sub = line.substring(line.indexOf("msLevel=")+9);	/* +9 moves the length of msLevel=" */
				sub = sub.substring(0, sub.indexOf('\"'));	
				scanHeader.msLevel = Integer.parseInt(sub);      
			}
			if (line.contains("peaksCount=\"")){
				String sub = line.substring(line.indexOf("peaksCount=")+12);	/* +12 moves the length of peaksCount=" */
				sub = sub.substring(0, sub.indexOf('\"'));	
				scanHeader.peaksCount = Integer.parseInt(sub); 
			}
			/**********************hanynag ISA************************/
			if (line.contains("centroided=\"")){
				String sub = line.substring(line.indexOf("centroided=")+12);	/* +12 moves the length of centroided=" */
				sub = sub.substring(0, sub.indexOf('\"'));	
				scanHeader.centroided = Integer.parseInt(sub); 
			}
			/**********************hanynag ISA************************/

			if (line.contains("retentionTime=\"")){
				String sub = line.substring(line.indexOf("retentionTime=")+17);	/* +17 moves the length of retentionTime=PT[0-9]+S" */
				sub = sub.substring(0, sub.indexOf('\"')-1);	
				scanHeader.retentionTime = Double.parseDouble(sub); 
			}
			if (line.contains("lowMz=\"")){
				String sub = line.substring(line.indexOf("lowMz=")+7);	/* +7 moves the length of lowMz=" */
				sub = sub.substring(0, sub.indexOf('\"')-1);	
				scanHeader.lowMZ = Double.parseDouble(sub); 
			}
			if (line.contains("highMz=\"")){
				String sub = line.substring(line.indexOf("highMz=")+8);	/* +8 moves the length of highMz=" */
				sub = sub.substring(0, sub.indexOf('\"')-1);	
				scanHeader.highMZ = Double.parseDouble(sub); 
			}
			if (line.contains("basePeakMz=\"")){
				String sub = line.substring(line.indexOf("basePeakMz=")+12);	/* +12 moves the length of basePeakMz=" */
				sub = sub.substring(0, sub.indexOf('\"')-1);	
				scanHeader.basePeakMz = Double.parseDouble(sub); 
			}
			if (line.contains("basePeakIntensity=\"")){
				String sub = line.substring(line.indexOf("basePeakIntensity=")+19);	/* +19 moves the length of basePeakIntensity=" */
				sub = sub.substring(0, sub.indexOf('\"')-1);	
				scanHeader.basePeakIntensity = Double.parseDouble(sub); 
			}
			if (line.contains("totIonCurrent=\"")){
				String sub = line.substring(line.indexOf("totIonCurrent=")+15);	/* +15 moves the length of totIonCurrent=" */
				sub = sub.substring(0, sub.indexOf('\"')-1);	
				scanHeader.totIonCurrent = Double.parseDouble(sub); 
			}
			/*
			* read precursor mass
			*/
			/**************************** slee ****************************/
			if (line.contains("precursorCharge=\"")){
				String sub = line.substring(line.indexOf("precursorCharge=")+17);	/* +17 moves the length of precursorCharge=" */
				sub = sub.substring(0, sub.indexOf('\"'));	
				scanHeader.precursorCharge = Integer.parseInt(sub); 
			}
			/**************************** slee ****************************/

			if (line.contains("activationMethod"))
			{

				/*
				* Find end of tag.
				*/
				while (line.indexOf('>') > 0)
				{
					line = br.readLine();
				}
				/*
				 * Skip past white space.
				 */
				line = line.trim();
				String sub = line.substring(line.indexOf('>')+1);
				sub = sub.substring(0, sub.indexOf('<'));	
				scanHeader.precursorMZ = Double.parseDouble(sub); 
			}
		}
	}
	
	/****************************************************************
	* Reads the MS level of the scan.				*
	* !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
	*    RETURNING !!						
	 * @throws IOException *
	***************************************************************/
	int  readMsLevel(File pFI, int lScanIndex) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(pFI));
		br.skip(lScanIndex);
		String line;
		int mslev = -1;
		while ( (line=br.readLine()) != null )
		{
			if (line.contains("<peaks "))
				break;
			if (line.contains("msLevel=\"")){
				String sub = line.substring(line.indexOf("msLevel=")+9);	/* +9 moves the length of msLevel=" */
				sub = sub.substring(0, sub.indexOf('\"'));	
				mslev = Integer.parseInt(sub);      
			}
		}
		return mslev;
	}
	
	/****************************************************************
	* Reads startMz and endMz of the scan.				*
	* Returns null if startMz was not set. Don't free the memory!	*
	* !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
	*    RETURNING !!						
	 * @throws IOException 
	 * @throws NumberFormatException *
	***************************************************************/
	double readStartMz(File pFI, int lScanIndex) throws NumberFormatException, IOException{
		BufferedReader br = new BufferedReader(new FileReader(pFI));
		br.skip(lScanIndex);
		String line;
		double startMz = 1.E6;
		while ( (line=br.readLine()) != null )
		{
			if (line.contains("<peaks "))
				break;
			if (line.contains("startMz=\"")){
				String sub = line.substring(line.indexOf("startMz=")+9);	/* +9 moves the length of startMz=" */
				sub = sub.substring(0, sub.indexOf('\"'));	
				startMz = Double.parseDouble(sub);      
			}
		}
		return startMz;
	}
	
	/****************************************************************
	* Reads startMz and endMz of the scan.				*
	* Returns null if startMz was not set. Don't free the memory!	*
	* !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
	*    RETURNING !!						
	 * @throws IOException 
	 * @throws NumberFormatException *
	***************************************************************/
	double readEndMz(File pFI, int lScanIndex) throws NumberFormatException, IOException{
		BufferedReader br = new BufferedReader(new FileReader(pFI));
		br.skip(lScanIndex);
		String line;
		double endMz = 0.0;
		while ( (line=br.readLine()) != null )
		{
			if (line.contains("<peaks "))
				break;
			if (line.contains("endMz=\"")){
				String sub = line.substring(line.indexOf("endMz=")+7);	/* +7 moves the length of endMz=" */
				sub = sub.substring(0, sub.indexOf('\"'));	
				endMz = Double.parseDouble(sub);      
			}
		}
		return endMz;
	}
	
	/****************************************************************
	* Reads RT of the scan.					*
	* Return a String that becomes property of the caller!		*
	* Returns null if the RT was not set. Don't free the memory!	*
	* !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
	*    RETURNING !!						
	 * @throws IOException *
	***************************************************************/
	String readRT(File pFI, int lScanIndex) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(pFI));
		br.skip(lScanIndex);
		String line;
		String rt = null;
		while ( (line=br.readLine()) != null )
		{
			if (line.contains("<peaks "))
				break;
			if (line.contains("retentionTime=\"")){
				String sub = line.substring(line.indexOf("retentionTime=")+15);	/* +15 moves the length of retentionTime=PT[0-9]+S" */
				sub = sub.substring(0, sub.indexOf('\"')-1);	
				rt = sub;
			}
		}
		return rt;
	}
	
	/****************************************************************
	* READS the base64 encoded list of peaks.			*
	* Return a float* that becomes property of the caller!		*
	* The list is terminated by -1					*
	* !! THE STREAM IS NOT RESET AT THE INITIAL POSITION BEFORE	*
	*    RETURNING !!						*
	***************************************************************/
//	double[] readPeaks(File pFI, int lScanIndex){
//		int  n;
//		int  peaksCount;
//		int  peaksCountLength;
//		int  peaksLen;       // The length of the base64 section
//		int  trail;
//		String szPeaksCount;
//		String line;
//		double pPeaks[];
//
//		String pData;
//		String pBeginData;
//		String pDecoded;
//		
//		BufferedReader br = new BufferedReader(new FileReader(pFI));
//		br.skip(lScanIndex);
//
//		// Get the num of peaks in the scan and allocate the space we need
//		line = br.readLine();
//		while (!(line.contains("peaksCount=\"")))
//		{
//			line = br.readLine();
//		}
//		
//		// We need to move forward the length of peaksCount="
//		String sub = line.substring(line.indexOf("peaksCount=")+12);
//		sub = sub.substring(0, sub.indexOf('\"'));
//		peaksCount = Integer.parseInt(sub);
//		
//		// We add 100 to compensate for initial white spaces and the opening of the element
//		peaksLen = peaksCount * (64 / 3) + 100;
//		if (peaksCount == 0)
//		{ // No peaks in this scan!!
//			return null;
//		}
//		
//		while (!(line.contains("peaks")))
//		{
//			line = br.readLine();
//		}
//
//		// skip <peaks precision="xx">
//		//pBeginData += 21;
//		while (!(line.contains(">")))
//		{
//			line = br.readLine();
//		}
//
//		// Base64 decoding
//		Base64.b64_decode_mio(pDecoded, pBeginData);
//
//
//		// And byte order correction
//
//		for (n = 0; n < (2 * peaksCount); n++)
//		{
//			((unsigned long *)pPeaks)[n] = ntohl((unsigned long)((unsigned long *)pDecoded)[n]);
////			printf("pPeaks [%d] = %ld \n", n, pPeaks[n]);
////			printf("pDecoded [%d] = %ld \n", n, pDecoded[n]);
//		}
//
////		pDecoded[peaksCount * 8] = '\0';
//
//		free(pData);
////		free(pDecoded);
//
//		return (pPeaks);
//	}
	
	/*
	* walk through each scan to find overall lowMZ, highMZ
	* sets overall start and end times also
	*/
//	void readRunHeader(File pFI, int[] pScanIndex, RunHeaderStruct runHeader, int iLastScan){
//		int  i;
//		ScanHeaderStruct scanHeader = new ScanHeaderStruct();
//		double startMz = 0.0;
//		double endMz = 0.0;
//
//		readHeader(pFI, pScanIndex[0], scanHeader);
//
//
//
//		/*
//		* initialize values to first scan
//		*/
//		// start at scan 1 or 0 ? 
//		i = 1;
//		runHeader.lowMZ = scanHeader.lowMZ;
//		runHeader.highMZ = scanHeader.highMZ;
//		runHeader.dStartTime = scanHeader.retentionTime;
//		runHeader.startMZ = readStartMz(pFI, pScanIndex[i]);
//		runHeader.endMZ = readEndMz(pFI, pScanIndex[i]);
//
//		for (i = 1; i < iLastScan; i++)
//		{
//			readHeader(pFI, pScanIndex[i], scanHeader);
//
//			if (scanHeader.lowMZ < runHeader.lowMZ)
//				runHeader.lowMZ = scanHeader.lowMZ;
//			if (scanHeader.highMZ > runHeader.highMZ)
//				runHeader.highMZ = scanHeader.highMZ;
//			if ((startMz = readStartMz(pFI, pScanIndex[i])) < runHeader.startMZ)
//				runHeader.startMZ = startMz;
//			if ((endMz = readEndMz(pFI, pScanIndex[i])) > runHeader.endMZ)
//				runHeader.endMZ = endMz;
//		}
//
//		runHeader.dEndTime = scanHeader.retentionTime;
//	}

	String setTagValue(String text, int maxlen, String lead, String tail){
		String result = null;
		String term = null;
		String storage = null;
		int len = maxlen - 1;
		int id_lead = text.indexOf(lead);
		if (id_lead > 0 )
		{
			result = text.substring(id_lead);
			int id_tail = result.substring(lead.length()).indexOf(tail);
			if (id_tail > 0)
			{
				term =  result.substring(lead.length()).substring(id_tail);
				if (result.length() - term.length() - lead.length() < len)
					len = result.length() - term.length() - lead.length();
				
				storage = result.substring(lead.length()).substring(0,len);
			} // if term
		}
		return storage;
	}
		
	InstrumentStruct getInstrumentStruct(File pFI) throws IOException {
		InstrumentStruct output = null;
		String result = null;
		String term = null;
		int found[] = { 0, 0, 0, 0, 0, 0 };
		String line;
		/**********************hanynag ISA************************/
		String storeData;
		/**********************hanynag ISA************************/

		output.manufacturer = null;
		output.model = null;
		output.ionisation = null;
		output.analyzer = null;
		output.detector = null;
		/**********************hanynag ISA************************/
		output.dataprocessing = null;            // default ดย profile
		/**********************hanynag ISA************************/

		
		BufferedReader br = new BufferedReader(new FileReader(pFI));
		line = br.readLine();

		while (line != null)  /* this should not be needed if index offset points to correct location */
		{
			if( !line.contains("<msInstrument") && !line.contains("<dataProcessing") )
				break;
			line = br.readLine();
		}
		/**********************hanynag ISA************************/
		while ( line != null)
		{
			if( line.contains("</msInstrument") )
				break;
			
			if (found[0] == 0)
			{
				if (line.contains("<msManufacturer")){
					result =  line.substring(line.indexOf("<msManufacturer"));
					if( (output.manufacturer =	setTagValue(result, INSTRUMENT_LENGTH, "value=\"", "\"")) != null )
						found[0] = 1;
				}
			}
			if (found[1] == 0)
			{
				if (line.contains("<msModel")){
					result =  line.substring(line.indexOf("<msModel"));
					if( (output.model =			setTagValue(result, INSTRUMENT_LENGTH, "value=\"", "\"")) != null )
						found[1] = 1;
				}
					
			}
			if (found[2] == 0)
			{
				if (line.contains("<msIonisation")){
					result =  line.substring(line.indexOf("<msIonisation"));
					if( (output.ionisation =	setTagValue(result, INSTRUMENT_LENGTH, "value=\"", "\"")) != null )
						found[2] = 1;
				}
			}
			if (found[3] == 0)
			{
				if (line.contains("<msMassAnalyzer")){
					result =  line.substring(line.indexOf("<msMassAnalyzer"));
					if( (output.analyzer =		setTagValue(result, INSTRUMENT_LENGTH, "value=\"", "\"")) != null )
						found[3] = 1;
				}
			}
			if (found[4] == 0)
			{
				if (line.contains("<msDetector")){
					result =  line.substring(line.indexOf("<msDetector"));
					if( (output.detector =		setTagValue(result, INSTRUMENT_LENGTH, "value=\"", "\"")) != null )
						found[4] = 1;
				}
			}

			line = br.readLine();

		} // while

		line = br.readLine();
		storeData = line;

		if (storeData.indexOf('0') > 0 ){
			output.dataprocessing[0] = PROFILE_MODE;
		}
		else if (storeData.indexOf('1') > 0 )
		{
			output.dataprocessing[0] = CENTROIDED_MODE;
		}
		br.close();
		return output; // no data

					   /**********************hanynag ISA************************/

	}
	
}
