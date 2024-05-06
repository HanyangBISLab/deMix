package kr.ac.hanyang.bislab.demix.distribution;

import java.util.ArrayList;

public class Fitting {
	
	static double scaling_max_ratio = 1.1, scaling_grid_ratio = 0.01, dynamic_match_tol = 0.2;
		
	public static DDMatch calculateMPC(ArrayList<Double> experimental, ArrayList<Double> theoritical, int theo_max_index) {

		
		//normalized by max_peak in theoretical
		double comp_ratio = experimental.get(theo_max_index)/theoritical.get(theo_max_index);				
		ArrayList<Double> normDeuDistT= new ArrayList<Double>();
		for(int j=0 ; j<experimental.size();j++) {
			normDeuDistT.add(theoritical.get(j)*comp_ratio);
		}
		
		double scaling_ratio= scaling_max_ratio;

		int 	max_mpc= 0;
		double 	mpc_error= 100000., mpc_fit_scale= 0.;
		
		int fitting_side= 0;// left, 1: right
		double interpreted = 0., total_error= mpc_error;
				
		while( scaling_ratio > 0 ){
			
			ArrayList<Double> comp_theo_dist = new ArrayList<Double>();
			for(int j=0 ; j<normDeuDistT.size();j++) {
				comp_theo_dist.add( normDeuDistT.get(j)*scaling_ratio );
				
			}
			
			int cur_mpc= 0;
			double cur_mpc_error= 0;
			double cur_interpreted= 0., cur_area= 0., cur_error= 0.;
			
			double left_error = 0, right_error=0;
			int cur_side = 0;
			
			int overlapped= Math.min(experimental.size(), comp_theo_dist.size());
			
			for(int j=0 ; j<overlapped; j++){
				
				double expInten= experimental.get(j), theInten= comp_theo_dist.get(j);
				double temp =  expInten - theInten;
				
				cur_error += temp*temp;	
			
				if( temp > 0 ){
					if( j < theo_max_index ) left_error += temp;
					else if( j > theo_max_index ) right_error += temp;
				}
		
				cur_area += Math.max(expInten, theInten);
				if( expInten < 0.01 || theInten < 0.01 ) continue;
				
				cur_interpreted += Math.min(expInten, theInten);
				
				if( Math.abs(temp) <= theInten*dynamic_match_tol ) {
					cur_mpc++;	
					cur_mpc_error += temp*temp;
				}
			}	
			
			if( left_error > right_error ) cur_side = 1;
			
			if( max_mpc < cur_mpc || 
					( max_mpc == cur_mpc && cur_mpc_error < mpc_error ) ) {				
				max_mpc = cur_mpc;
				mpc_error = cur_mpc_error;
				
				interpreted = cur_interpreted/cur_area;
				total_error = cur_error;
				fitting_side= cur_side;
				mpc_fit_scale= scaling_ratio;				
			}
			scaling_ratio -= scaling_grid_ratio;
		}
		
		return new DDMatch(0, max_mpc, mpc_error, mpc_fit_scale, interpreted, total_error, fitting_side);
		
	}
	
	
	
	
	
	
}





