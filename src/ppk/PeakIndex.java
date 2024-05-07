package ppk;

import java.util.List;

import utils.Constants;

/**
 * porting form C++ code to java code
 * @author user_Lee
 *
 */
public class PeakIndex {
	/**
	 * 
	 * @param vec
	 * @param mz_val is the value we are looking for.
	 * @param start_index index around which we are looking for the specified point.
	 * @return returns the index of the point that is closest to the specified value.
	 */
	public static int getNearest(List<Peak> vec, double mz_val, int start_index)
	{
		// we're going to use continuity here, look at the difference 
		// between consecutive points and estimate how much further we have to 
		// go and start there. 
		int num_pts = (int)vec.size() - 1;

		if (mz_val >= vec.get(num_pts).mdbl_mz)
			return num_pts;
		if (mz_val < vec.get(0).mdbl_mz)
			return 0;

		double distance_to_go = mz_val - vec.get(start_index).mdbl_mz;
		double step;
		if (start_index < num_pts)
			step = vec.get(start_index + 1).mdbl_mz - vec.get(start_index).mdbl_mz;
		else
			step = vec.get(start_index).mdbl_mz - vec.get(start_index - 1).mdbl_mz;


		int move_by = (int)(distance_to_go / step);
		int next_index = start_index + move_by;

		if (next_index < 0)
			next_index = 0;
		if (next_index > num_pts)
			next_index = num_pts - 1;

		if (mz_val >= vec.get(next_index).mdbl_mz)
			return lookRight(vec, mz_val, next_index);
		else
			return lookLeft(vec, mz_val, next_index);
	}
	
	/**
	 * 
	 * @param vec
	 * @param mz_val is the value we are looking for.
	 * @param start_index minimum index of the point.
	 * @param stop_index maximum index of the point.
	 * @return returns the index of the point that is closest to the specified value.
	 */
	public static int getNearestBinary(List<Peak> vec, double mz_val, int start_index, int stop_index)
	{
		double min_val, max_val, mid_val, mid_next_val;
//		System.out.println(vec.size());
		if (vec.get(start_index).mdbl_mz > mz_val)
			return start_index;
		if (vec.get(stop_index).mdbl_mz < mz_val)
			return stop_index;

		int mid_index;
		while (true)
		{
			min_val = vec.get(start_index).mdbl_mz;
			max_val = vec.get(stop_index).mdbl_mz;

			if (Math.abs(stop_index - start_index) <= 1 && mz_val >= min_val && mz_val <= max_val)
			{
				if (Constants.absolute(min_val - mz_val) < Constants.absolute(max_val - mz_val))
					return start_index;
				return stop_index;
			}

			double ratio = ((max_val - mz_val) * 1.0) / (max_val - min_val);
			mid_index = (int)(start_index * ratio + stop_index * (1 - ratio) + 0.5);
			if (mid_index == start_index)
				mid_index = start_index + 1;
			else if (mid_index == stop_index)
				mid_index = stop_index - 1;

			mid_val = vec.get(mid_index).mdbl_mz;
			if (mid_val >= mz_val)
			{
				stop_index = mid_index;
			}
			else if (mid_index + 1 == stop_index)
			{
				if (Constants.absolute(mid_val - mz_val) < Constants.absolute(max_val - mz_val))
					return mid_index;
				return stop_index;
			}
			else
			{
				mid_next_val = vec.get(mid_index + 1).mdbl_mz;
				if (mz_val >= mid_val && mz_val <= mid_next_val)
				{
					if (mz_val - mid_val < mid_next_val - mid_val)
						return mid_index;
					return mid_index + 1;
				}
				start_index = mid_index + 1;
			}
		}
	}
	
	/**
	 * does a search for the given value by doing a linear scan to the left of the given index
	 * @param vec 
	 * @param mz_val is the value we are looking for.
	 * @param start_index index of the peak to the left of which we are scanning.
	 * @return returns the index of the point that is closest to the specified value.
	 */
	private static int lookLeft(List<Peak> vec, double mz_val, int start_index){
		// mv_val <= vec[start_index] so start moving index further left.
		int nearest_index = start_index;
		int next_index = start_index;

		if (next_index == 0)
			return 0;

		double next_val = vec.get(next_index).mdbl_mz;
		double best_distance = Constants.absolute(mz_val - next_val);

		while (next_val > mz_val)
		{
			next_index--;
			next_val = vec.get(next_index).mdbl_mz;
			double dist = Constants.absolute(next_val - mz_val);
			if (dist < best_distance)
			{
				best_distance = dist;
				nearest_index = next_index;
			}
			if (next_index == 0)
				break;
		}
		return nearest_index;
	}
	
	/**
	 * 
	 * @param vec
	 * @param mz_val is the value we are looking for.
	 * @param start_index index of the peak to the right of which we are scanning.
	 * @return returns the index of the point that is closest to the specified value.
	 */
	private static int lookRight(List<Peak> vec, double mz_val, int start_index)
	{
		// mv_val >= vec[start_index] so start moving index further right.
		int nearest_index = start_index;
		int next_index = start_index;
		int num_pts = (int)vec.size();

		if (next_index >= num_pts - 1)
			return num_pts - 1;

		double next_val = vec.get(next_index).mdbl_mz;
		double best_distance = Constants.absolute(mz_val - next_val);

		// we've gone back too far, posibly. Move pas the mz_val and return that value. 
		while (next_val < mz_val)
		{
			next_index++;

			next_val = vec.get(next_index).mdbl_mz;
			double dist = Constants.absolute(next_val - mz_val);
			if (dist < best_distance)
			{
				best_distance = dist;
				nearest_index = next_index;
			}

			if (next_index == num_pts - 1)
				break;
		}
		return nearest_index;
	}
}
