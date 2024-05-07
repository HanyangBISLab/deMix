package rds;

import utils.Constants.FILE_TYPE;

public class ReaderFactory {
	public static RawData getRawData(FILE_TYPE file_type){
		String header_n = "acqu";
		switch (file_type)
		{
			
		case FINNIGAN:
			break;
		case MZXMLRAWDATA:
			MZXmlRawData mzxml_raw_data;
			mzxml_raw_data = new MZXmlRawData();
			return mzxml_raw_data;
		default:
			break;
		}
		return null;
	}
}
