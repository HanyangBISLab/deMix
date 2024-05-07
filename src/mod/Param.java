package mod;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

import utils.Constants;

public class Param {
	public Param() {
		init();
	}

	private void init() {
		destroy();
		Constants._tbl.put("AdvDeconv::MaxMissPeak", "3");
		Constants._tbl.put("DeconvPep::MaxCharge", "10");
		Constants._tbl.put("AdvDeconv::ScanNoModifier", "0");
		Constants._tbl.put("AdvDeconv::MaxAbundancePeak", "3");
		Constants._tbl.put("AdvDeconv::MassErr", "1.0E-05");
		Constants._tbl.put("AdvDeconv::ThClustExt", "-2.0");
		Constants._tbl.put("DeconvPep::ThScore", "-0.1");
		Constants._tbl.put("PeakPicking::Partition", "5");
		Constants._tbl.put("PeakPicking::SNRThreshold", "0.0");
		Constants._tbl.put("PeakPicking::BackgroundRatio", "0.0");
		Constants._tbl.put("PeakPicking::BackgroundRatio2", "0.0");
		Constants._tbl.put("PeakPicking::PrintPeak", "true");
	}

	public String get(String key) {
		String temp = Constants._tbl.get(key);
		if (temp != null) {
			temp = temp.toUpperCase(); 
		}
		return temp;
	}

	public List<String> getList(String key) {
		List<String> temp = Constants._tbl2.get(key);
		return temp;
	}

//	public void get(String key, String val) {
//		map<string, char*>::iterator it = _tbl.find(key);
//		if (it != _tbl.end()) {
//			sscanf(it->second, "%c", &val);
//		}
//	}
//
//	void getList(const char* key, int& val) {
//		map<string, char*>::iterator it = _tbl.find(key);
//		if (it != _tbl.end()) {
//			sscanf(it->second, "%d", &val);
//		}
//	}
//
//	void get(const char* key, double& val) {
//		map<string, char*>::iterator it = _tbl.find(key);
//		if (it != _tbl.end()) {
//			sscanf(it->second, "%lf", &val);
//		}
//	}

	public void readFile(String file_name) {
		try{
			File fp = new File(file_name);
			
			
			BufferedReader br = new BufferedReader(new FileReader(fp));
			
			String line;
			String bk = null, bv = null;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if( line.isEmpty() )
					continue;
				
				String[] split = line.split("\\s+");
				if( split.length == 2 ){			//Read 'parameter' in parameter file.
					bk = new String(split[0]);
					bv = new String(split[1]);
						
//					System.out.println(line);
					Constants._tbl.put(bk, bv);
				}
				else if ( split.length ==  1) {		//Read 'file list' in parameter file.
					List<String> slist = new ArrayList<>(); 
					while( (line = br.readLine()) != null ) {
						line = line.trim();
						if( line.trim().isEmpty() )
							break;
						
//						System.out.println(line);
						bk = new String(split[0]);
						line = line.trim();
						
						slist.add(line);
					}
					Constants._tbl2.put(bk, slist);
				}
			}
			br.close();
		}catch(IOException e){
			System.out.println("Can't open Parameter file.");
			return;
		}

	}
	
	public void dump() {		//print screen
		List<String> i1 = new ArrayList<>(); 
		i1.addAll(Constants._tbl.keySet());
		for (int i = 0 ; i < i1.size(); i++) {
			System.out.println(i1.get(i) + " = "+ Constants._tbl.get(i1.get(i)) + "\n");
		}
		
		List<String> i2 = new ArrayList<>();
		for (int i = 0; i < i2.size(); i++) {
			System.out.println(i2.get(i) + " = \n");
			List<String> l2 = Constants._tbl2.get(i2.get(i));
			for (int j = 0; j < l2.size(); j++) {
				System.out.println("\t" + l2.get(j) + "\n");
			}
		}
		System.out.println();
	}

	void destroy() {
		Constants._tbl = new TreeMap<>();
		Constants._tbl2 = new TreeMap<>();
	}
}
