package kr.ac.hanyang.bislab.demix;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;


public class deMix {

	public static void main(String[] args) throws Exception {
		
	
		System.out.println("************************************************************************************");
		System.out.println("deMix (version 3.20) Hydrogen Deuterium eXchange Analysis");
		System.out.println("Release Date: May 1, 2024");
		System.out.println("Hanyang University, Seoul, Korea");
		System.out.println("************************************************************************************");
		System.out.println();
		
		String paramfile= null, result= null;
		boolean requirement = false;

		for (int i=0; i<args.length; i++){			
			if( args[i].equals("-i") ){
				i++;
				if (i < args.length) {
					requirement = true;
					paramfile= args[i];
				}
			}
			else if( args[i].equals("-o") ){
				i++;
				if (i < args.length) result = args[i];
			}
		}	
		
		if( !requirement ){
			printUsage();
			return;
		}//*/
				
		long startTime= System.currentTimeMillis();
		
		deMix hdx = new deMix(paramfile);
		if( !hdx.validateParams() ) return;

		System.out.println();
		hdx.run();
		
		if( result == null )
			result = hdx.pept_input_file.substring(0, hdx.pept_input_file.lastIndexOf('.'));
		

		
		System.out.println("[deMix] Elapsed Time : " + (System.currentTimeMillis()-startTime)/1000. + " Sec" );
		System.out.println("[deMix] BYE!!");
	}
	
	String project_name;
	
	String pept_input_file;
	String fasta_input_file;

	ArrayList<Group> groupList;

	double PPMTolerance	= 10;
	double DATolerance	= 0.1;
	
	public void run() throws Exception {
		
		String projectfile= pept_input_file.replace('\\', '/');//for project file
		projectfile= projectfile.substring(0, projectfile.lastIndexOf('/')+1)+project_name+".dmxj";
		PrintWriter out = new PrintWriter( new BufferedWriter( new FileWriter(projectfile) ) );
		printParams(out);
		
		for(int gi=0; gi<groupList.size(); gi++) {
			
			System.out.println("Analyzing group data: " + groupList.get(gi).groupName);
			
			HDXAnalysis anal= new HDXAnalysis(this, gi);	
			anal.decode();
			
			System.out.println();		

			out.println("#Group::"+groupList.get(gi).groupName );
			out.println("#HDXProfile::");
			anal.printHDXProfile(out);
			out.println();
			out.println("#DdeuAnal::");
			anal.printDdeuAnalysis(out);
			out.println();
		}
		out.close();
		
	}
	
	
	public void printParams(PrintWriter out) throws IOException {
		out.println("#Params::");
		out.println("Project="+project_name);
		out.println("Peptide="+pept_input_file);

		for( Group g : groupList ) { 
			out.println( String.format("CTRLData[%s]=%s", g.groupName, g.ctrl_ms_file) );
		}
		
		for( Group g : groupList ) { 
			for( Map.Entry<String, String> entry : g.hdx_ms_files.entrySet() ) {
				out.println( String.format("HDXData=%s[%s],%s", entry.getKey(), g.groupName, entry.getValue()) );
			}
		}
		
		if( fasta_input_file != null )
			out.println("Protein="+fasta_input_file);
		out.println("MassTolerance="+(PPMTolerance==0? DATolerance+"da" : PPMTolerance+"ppm") );
		out.println();
	}
	
	public boolean validateParams() throws IOException {
		
		if( pept_input_file == null ) {
			System.out.println("Notice : Peptide list file should be specified in 'Peptide=' parameter.");
			return false;
		}
		else if( groupList.size() == 0 ) {
			System.out.println("Notice : HDX Spectra files should be specified.");
			return false;
		}
		else if( PPMTolerance == DATolerance ) {
			System.out.println("Notice : " + "MassTolerance valus is wrong.");
			return false;
		}
	
		System.out.println("Project name : " + project_name);
		
		System.out.println("Input peptide : " + pept_input_file);
		
		for( Group g : groupList ) { 
		
			System.out.println("Group name : " + g.groupName);
			System.out.println(" -ctrl spectrum : " + g.ctrl_ms_file);
	
			System.out.println(" -hdx spectra : ");
			int hi=1;
			for( Map.Entry<String, String> entry : g.hdx_ms_files.entrySet() ) {
				System.out.println( "  " + (hi++) + ". " + entry.getKey() + ", " + entry.getValue() );
			}
		}
		System.out.println("Input protein : " + fasta_input_file);
		
		System.out.println("Input mass tolerance : " + (PPMTolerance==0? String.format("%.4f da", DATolerance) : String.format("%.4f ppm", PPMTolerance)) );
		
		return true;
	}
	
	public deMix( String paramfile ) throws IOException {
		
		project_name= "deMix_proj_"+System.currentTimeMillis();
		
		System.out.println("Reading parameter.....");
		
		String s;
		
		groupList= new ArrayList<Group>();
		HashMap<String, Integer> groupIndex= new HashMap<String, Integer>();
		
		BufferedReader in = new BufferedReader( new FileReader(paramfile) );				
		while((s = in.readLine()) != null) {
			
			s= s.trim();
			if( s.startsWith("#") ) continue;
			
			String values[]= s.split("=");
			if( values.length != 2 ) continue;
			
			String paramName	= values[0].trim();
			String paramValue	= values[1].trim();

			if( "Project".equalsIgnoreCase(paramName) ){
				project_name = paramValue;			
			}
			
			else if( "Peptide".equalsIgnoreCase(paramName) ){
				pept_input_file = paramValue;			
			}
			
			else if( "Protein".equalsIgnoreCase(paramName) ){
				fasta_input_file = paramValue;
			}

			else if( paramName.startsWith("CTRLData[") ){
			
				String gname= paramName.substring(paramName.lastIndexOf('[')+1, paramName.length()-1);
				groupIndex.put(gname, groupList.size());
				groupList.add( new Group(gname, paramValue) );	
			}
			
			else if( "HDXData".equalsIgnoreCase(paramName) ){
				
				String sub[]= paramValue.split(",");
				if( sub.length != 2 ) continue;
				
				String label= sub[0].trim();
				String time= label.substring(0, label.lastIndexOf('['));
				String gname= label.substring(label.lastIndexOf('[')+1, label.length()-1);
				groupList.get(groupIndex.get(gname)).addHDXFile(time, sub[1].trim());
				
			}

			else if( "MassTolerance".equalsIgnoreCase(paramName) ){
				
				String err= paramValue.toLowerCase();
				if( err.endsWith("ppm") ) {
					PPMTolerance= Double.parseDouble(err.substring(0, err.lastIndexOf("ppm")));			
					DATolerance= 0;
				}
				else if( err.endsWith("da") ) {
					DATolerance= Double.parseDouble(err.substring(0, err.lastIndexOf("da")));
					PPMTolerance= 0;
				}
				else {
					PPMTolerance = DATolerance = 0.;			
				}
			}
		}
		in.close();
	}
	
	protected static void printUsage(){	
		System.out.println("--------------------------------------------------------------------------");
		System.out.println("Usage    : java -jar deMix.jar <options> <file>");
		System.out.println("Options  :");
		System.out.println("  -i <parameter_file> : Config file path for search parameters [Required]");
		System.out.println("  -o <results_name>   : Output file path for search results [Optional]");
		System.out.println();
		System.out.println("Example1 : java -jar deMix.jar -i foo.txt -o hdx_results_name");
		System.out.println("Example2 : java -jar deMix.jar -i foo.txt");			
		System.out.println();
		System.out.println("Citation :");
		System.out.println("deMix: Decoding Deuterated Distributions from Heterogeneous Protein States via HDX-MS.");
		System.out.println("Seungjin Na, Jae-Jin Lee, Jong Wha Joo, Kong-Joo Lee and Eunok Paek.");
		System.out.println("Scientific Reports, 2019, 9, 3176.");	
		System.out.println();
		System.out.println("For feedback, questions and comments,");
		System.out.println("    contact Seungjin Na (sna@hanyang.ac.kr).");
	}
	

}











