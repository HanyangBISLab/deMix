package kr.ac.hanyang.bislab.demix.protein;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class ProxDB extends ArrayList<Prox> {
	
	int sizeOfEntries;
	int sizeOfResidues;	
	String name, orgName;
	
	public int getSizeOfResidues() { return sizeOfResidues; }
	public String getName() { return orgName; }
	
	public void setSizeOfEntries(int x) { sizeOfEntries = x; }
	public void setSizeOfResidues(int x) { sizeOfResidues = x; }
	
	public void readFasta(String fileName) throws Exception
	{
		if( fileName == null ) {
			System.out.println( "The protein DB was not specified." );
			return;
		}
		orgName = fileName;
		int residues= 0;
		try 
		{		
			if( fileName.lastIndexOf('.') > -1 )
				name= fileName.substring(0, fileName.lastIndexOf('.'));
			else name = fileName;
			
			BufferedReader in = new BufferedReader(new FileReader(fileName));
			String s= "<";
			StringBuffer buffer;	
			
			while( s.startsWith(">") == false ) {
				s = in.readLine();
				if( s == null ) break;
			}
		
			while( s != null ) 
			{
				Prox protein= new Prox();				
				protein.setHeader( s.substring(1) );			
				
				buffer = new StringBuffer();
				while( (s = in.readLine()) != null ){
					if( s.startsWith(">") ) break;

					for( int aa=0; aa<s.length(); aa++ ){
						if( Character.isLetter( s.charAt(aa) ) )
							buffer.append( Character.toUpperCase(s.charAt(aa)) );
					}
				}
				if( buffer.length() < 3 ) continue; // check sequence
				
				protein.setSequence( buffer.toString() );
				residues += buffer.length();
				this.add( protein );					
			}
			in.close();						
		} 
		catch (FileNotFoundException e) 
		{
			System.out.println( "Cannot find the protein file, "+fileName );
			e.printStackTrace();
			System.exit(1);
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
			System.exit(1);
		}
		
		sizeOfEntries = this.size();
		sizeOfResidues = residues;	
		System.out.println( sizeOfEntries+" proteins / "  + sizeOfResidues+" residues" );
	}
	
	
	
	
	
	public Match2Protein match2Protein(String peptide) throws Exception {
		
		for( Prox pro : this ){		
			int index= pro.sequence.indexOf(peptide);
			if( -1 < index ) return new Match2Protein(pro.accession, index+1, index+peptide.length());
		}
		return new Match2Protein("~", 0, 0);
	}
	

	
	
}























