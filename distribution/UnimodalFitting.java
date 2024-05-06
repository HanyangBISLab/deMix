package kr.ac.hanyang.bislab.demix.distribution;

import java.util.ArrayList;
import java.util.Collections;

public class UnimodalFitting {
	
	public static ArrayList<DDMatch> run(ArrayList<Double> experimental, ArrayList<Double>[] calcDeuDist, int[] calc_max_index) {
		
		ArrayList<DDMatch> unimodalAnalysis= new ArrayList<DDMatch>();
		
		int exchangeableSites= calcDeuDist.length-1;
		double maxInter= 0.;
		for(int hlevel=0 ; hlevel<=exchangeableSites; hlevel++) {
			
			ArrayList<Double> deuDistT= calcDeuDist[hlevel];
			int c_max_index= calc_max_index[hlevel];
		
			if( experimental.get(c_max_index) < deuDistT.get(c_max_index)*0.1 ) continue;
			
			DDMatch curMatch= Fitting.calculateMPC(experimental, deuDistT, c_max_index);
			curMatch.setDeuterated(hlevel, 1., deuDistT);
			
			if( maxInter < curMatch.interperted ) maxInter= curMatch.interperted;
			unimodalAnalysis.add(curMatch);
		}
		maxInter -= 0.0001;
		for(DDMatch ddm : unimodalAnalysis) {
			if( ddm.interperted > maxInter ) ddm.mpc += 1;
		}
		
		Collections.sort(unimodalAnalysis);
		
		
		return unimodalAnalysis;
	}
	
	
}
