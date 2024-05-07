package kr.ac.hanyang.bislab.demix.protein;

public class Match2Protein implements Comparable<Match2Protein> {
	
	String id;
	int start;
	int end;
	
	public Match2Protein(String d, int s, int e) {
		id= d; 
		start= s;
		end= e;
	}
	
	public String toString() {
		return id+"\t"+start+"\t"+end;
	}

	public int compareTo(Match2Protein x) {
		
		if( id.compareTo(x.id) > 0 ) return 1;
		else if( id.compareTo(x.id) < 0 ) return -1;
		
		if( start > x.start ) return 1;
		else if( start < x.start ) return -1;
		
		if( end > x.end ) return 1;
		else if( end < x.end ) return -1;
		
		return 0;
	}
	
	
	
	
}
