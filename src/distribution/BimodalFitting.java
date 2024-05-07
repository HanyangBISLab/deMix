package kr.ac.hanyang.bislab.demix.distribution;

import java.util.ArrayList;
import java.util.HashSet;

public class BimodalFitting {
	
	private static double error_improvement = 0.1;
	
	public static DDMatch run(ArrayList<DDMatch> unimodalAnalysis, ArrayList<Double> experimental, ArrayList<Double>[] calcDeuDist) {

		DDMatch unimodalResult = unimodalAnalysis.get(0); 
		int HDXSiteLength= calcDeuDist.length-1;
		
		ArrayList<Integer> secHdx = new ArrayList<Integer>();
		int mpc2nd= 0, ev2nd= 0, mem2nd= 0;
		
		for(int ev=1; ev<unimodalAnalysis.size(); ev++){
		
			if( unimodalResult.side != unimodalAnalysis.get(ev).side && 
													mpc2nd <= unimodalAnalysis.get(ev).mpc ) {
				secHdx.add(unimodalAnalysis.get(ev).deuterated);
				mpc2nd = unimodalAnalysis.get(ev).mpc;
				ev2nd = ev;
				mem2nd++;
			}
		}
		
		if( mem2nd == 1 ){
			if( (unimodalResult.deuterated - unimodalAnalysis.get(ev2nd).deuterated) == 1 ){
				if( unimodalAnalysis.get(ev2nd).deuterated > 0 ) secHdx.add(unimodalAnalysis.get(ev2nd).deuterated-1);
			}
			else if( (unimodalResult.deuterated - unimodalAnalysis.get(ev2nd).deuterated) == -1 ){
				if( unimodalAnalysis.get(ev2nd).deuterated < HDXSiteLength ) secHdx.add(unimodalAnalysis.get(ev2nd).deuterated+1);
			}
		}

		DDMatch bestResult = unimodalAnalysis.get(0); 
		int mix1= unimodalResult.deuterated, mix2= unimodalResult.deuterated;		
		
		HashSet<String> hdxpair = new HashSet<String>();		
		for(int sec_h=0; sec_h < secHdx.size(); sec_h++){	
			
			int leftmix = unimodalResult.deuterated, rightmix = secHdx.get(sec_h);
			
			if( rightmix < leftmix ){
				leftmix = rightmix;
				rightmix = unimodalResult.deuterated;
			}
			
			int num_1st= Math.max(leftmix-1, 0);
			int num_2nd= Math.min(rightmix+1, HDXSiteLength);		

			for(int f=num_1st; f<=leftmix; f++){
				for(int s=num_2nd; s>=rightmix; s--){
					if( !hdxpair.add(f+"_"+s) ) continue;
					
					DDMatch tp_eval = mixture_fitting(experimental, calcDeuDist, f, s);
					if( tp_eval.mpc > bestResult.mpc || 
							( tp_eval.mpc == bestResult.mpc && tp_eval.mpc_err < bestResult.mpc_err ) ){					
						bestResult = tp_eval;
						mix1 = f;
						mix2 = s;
					}
				}
			}
		}

		double matPImprove = 2;
		if( mix1 == 0 || mix1 == 1 ) matPImprove = 1.5;
		if( mix2 == 0 || mix2 == 1 ) matPImprove = 1.5;
		
		if( ( bestResult.mpc < unimodalResult.mpc*matPImprove-1 ) || 
				(bestResult.error > unimodalResult.error * error_improvement) || 
					bestResult.mainPortion > 0.9 ){
			bestResult = unimodalResult; 
		}
				

		return bestResult;
	}
	
	
	
	
	
	
	
	
	
	
	
	public static DDMatch mixture_fitting(ArrayList<Double> experimental, ArrayList<Double>[] calcDeuDist, int hdx1st, int hdx2nd) {
		
		//generate bimodal dist
		ArrayList<Double> deuDist1stT= calcDeuDist[hdx1st];
		ArrayList<Double> deuDist2ndT= calcDeuDist[hdx2nd];
		
		DDMatch bestMatch= new DDMatch();
		double mixing = 1.;
		while( mixing > 0 ){
			
			ArrayList<Double> theo_mixed= new ArrayList<Double>();
			int mx_max_index = 0;
			double mx_max = -1;
			
			for(int j=0 ; j<deuDist2ndT.size();j++) {
				
				double mx = mixing*deuDist1stT.get(j) + (1-mixing)*deuDist2ndT.get(j);
				
				theo_mixed.add( mx );
				if( mx_max < mx ){
					mx_max = mx;
					mx_max_index = j;
				}
			}
		/*	for(int j=0 ; j<theo_mixed.size();j++) {
				if( theo_mixed.get(j) < 0.005 )
					theo_mixed.set(j, 0.);
			}//*/
			
			DDMatch curMatch= Fitting.calculateMPC(experimental, theo_mixed, mx_max_index);
			if( 0.5 < mixing ) curMatch.setDeuterated(hdx1st, mixing, hdx2nd, theo_mixed);
			else curMatch.setDeuterated(hdx2nd, 1-mixing, hdx1st, theo_mixed);
			
			if( bestMatch.compareTo(curMatch) == 1 ) {
				bestMatch= curMatch;
			}
			
			mixing -= 0.01;
		}
		
		return bestMatch;
	}
	
	
	
}












