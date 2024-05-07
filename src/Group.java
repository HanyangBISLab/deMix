package kr.ac.hanyang.bislab.demix;

import java.util.LinkedHashMap;

public class Group {

	String groupName;
	String ctrl_ms_file;
	LinkedHashMap<String, String> hdx_ms_files;
	
	public Group(String gname, String cname) {
		groupName= gname;
		ctrl_ms_file= cname;
		hdx_ms_files= new LinkedHashMap<String, String>();
	}
	
	public void addHDXFile(String time, String fname) {
		hdx_ms_files.put(time, fname);
	}
	
	
}




