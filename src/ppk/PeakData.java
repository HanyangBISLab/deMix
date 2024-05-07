package ppk;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class PeakData {
	public PeakData() {
		mptr_vect_peaks = new ArrayList<>();
		mvect_peak_tops = new ArrayList<>();
		mmap_pk_intensity_index = new TreeMap<>();
		mmap_pk_mz_index = new TreeMap<Double, Integer>();
		mmap_pk_mz_index_all = new TreeMap<Double, Integer>();
	}
	
//	public PeakData(PeakData a) {
//		mdbl_tic = a.mdbl_tic;
//		mvect_peak_tops.addAll(0, a.mvect_peak_tops);
//	}
	
	public void clear()
	{
		mptr_vect_peaks = new ArrayList<>();
		mvect_peak_tops = new ArrayList<>();
		mmap_pk_intensity_index = new TreeMap<>();
		mmap_pk_mz_index = new TreeMap<Double, Integer>();
		mmap_pk_mz_index_all = new TreeMap<Double, Integer>();

	}
	//! pointer to std::vector of mz's in the raw data
//	std::vector <double> *mptr_vect_mzs ; 
	//! pointer to std::vector of intensities in the raw data.  
//	std::vector<double> *mptr_vect_intensities 
	List<Peak> mptr_vect_peaks;
	//! pointer to the std::vector of temporary peaks that are used during the processing. 
	List<Peak> mvect_temp_peak_tops;
	//! multimap of indices of unprocessed peaks in PeakData::mvect_temp_peak_tops sorted in descending intensity. This helps fast retrieval of the next most intense unprocessed peak.
	/*!
	\remarks While the intensity of the peaks might actually be double, for the map, we only used the integral values.
	*/
	//, std::greater<int> 
	TreeMap<Integer, Integer> mmap_pk_intensity_index;
	//! multimap of indices of unprocessed peaks in PeakData::mvect_temp_peak_tops sorted in ascending m/z. This helps fast retrieval when looking for an unprocessed peak around an approximate m/z.
	Map<Double, Integer> mmap_pk_mz_index;
	//! multimap of indices of all peaks in PeakData::mvect_temp_peak_tops sorted in ascending m/z. This helps fast retrieval when looking for a peak(not only unprocessed ones) around an approximate m/z.
	Map<Double, Integer> mmap_pk_mz_index_all;
	//! The TIC for current scan. 
	double mdbl_tic;
	//! helps searching for specific values in any sorted std::vector. 
	PeakIndex mobj_pk_index;
	//! std::vector of peaks found in the data. It is recommended that this object not be touched
	// by calling functions. 
	public List<Peak> mvect_peak_tops;
	//! Gets the index^th peak in the std::vector PeakData::mvect_peak_tops.
	/*!
	\param index is the index of the peak in the std::vector PeakData::mvect_peak_tops.
	\return pk is assigned the peak which is at the index^th position in std::vector PeakData::mvect_peak_tops.
	*/
	public void getPeak(int index, Peak pk){
		pk = mvect_peak_tops.get(index);
	}
	
	//! Adds a peak to the std::vector of peaks PeakData::mvect_peak_tops
	/*!
	\param pk is the peak that we want to add to our std::vector of peaks.
	*/
	public void addPeak(Peak pk){
		mvect_peak_tops.add(pk);
	}
	
	//! Returns number of peaks.
	/*!
	\return number of peaks found.
	*/
	public int getNumPeaks(){
		return mvect_peak_tops.size();
	}
	
	//! Gets the most intense peak(whether or not it is processed) in the m/z range (mz - width to mz + width). The intensity returned is the intensity in the original raw data std::vector PeakData::mptr_vect_intensities
	/*!
	\param start_mz minimum m/z of the Peak.
	\param stop_mz mimum m/z of the Peak.
	\param pk is set to the peak that was found.
	\param exclude_mass is the mass we need to exclude in this search.
	\return returns true if a peak was found in the window (mz - width to mz + width) and false if not found.
	\note The returned peak has the intensity in the original raw data std::vector PeakData::mptr_vect_intensities.
	*/
//	public boolean getPeakFromAllOriginalIntensity(double start_mz, double stop_mz, Peak pk, double exclude_mass){
//		pk.mdbl_intensity = -10;
//		pk.mdbl_mz = 0;
//		boolean found = false;
//		
//		List<Double> mzs = new ArrayList<>();
//		mzs.addAll(mmap_pk_mz_index_all.keySet());
//		Collections.sort(mzs);
//		Iterator<Double> iter_mz = mzs.iterator();
//		double tmp = iter_mz.next();
//		while( tmp < start_mz ) tmp = iter_mz.next();
//		
//		while(true){
//			int peak_index = mmap_pk_mz_index_all.get(tmp);
//			double mz_val = tmp;
//			if (mz_val == exclude_mass)
//				continue;
//			if (mz_val > stop_mz)
//				return found;
//			int data_index = mvect_peak_tops.get(peak_index).mint_data_index;
//			if (mvect_peak_tops.get(data_index).in() >= pk.mdbl_intensity && mz_val >= start_mz)
//			{
//				double this_mz = mvect_peak_tops.get(peak_index).mdbl_mz;
//				pk = mvect_peak_tops.get(peak_index);
//				found = true;
//			}
//			
//			if( !iter_mz.hasNext() )
//				break;
//			else
//				tmp = iter_mz.next();
//		}
//		
//		return found;
//	}
	
	//! Gets the most intense peak(whether or not it is processed) in the m/z range (mz - width to mz + width). The intensity returned is the intensity in the original raw data std::vector PeakData::mptr_vect_intensities
	/*!
	\param start_mz minimum m/z of the Peak.
	\param stop_mz mimum m/z of the Peak.
	\param pk is set to the peak that was found.
	\return returns true if a peak was found in the window (mz - width to mz + width) and false if not found.
	\note The returned peak has the intensity in the original raw data std::vector PeakData::mptr_vect_intensities.
	*/
	public boolean getPeakFromAllOriginalIntensity(double start_mz, double stop_mz, Peak pk){
		pk.mdbl_intensity = -10;
		pk.mdbl_mz = 0;
		boolean found = false;

		List<Double> mzs = new ArrayList<>();
		mzs.addAll(mmap_pk_mz_index_all.keySet());
		Collections.sort(mzs);
		int i;
		double tmp = Double.MAX_VALUE;
		for( i = mzs.size()-1;;i-- ){
			if( tmp < stop_mz ){
				break;
			}
			
			tmp = mzs.get(i);
		}
		i++;
		int index= mzs.size()-1;
		while(true){
			double mz_val = mzs.get(index);
			int peak_index = mmap_pk_mz_index_all.get(mz_val);
			if (mz_val < start_mz)
				return found;
			if (mvect_peak_tops.get(peak_index).mdbl_intensity > pk.mdbl_intensity && mz_val >= start_mz && mz_val <= stop_mz)
			{
				double this_mz = mvect_peak_tops.get(peak_index).mdbl_mz;
				pk = mvect_peak_tops.get(peak_index);
				found = true;
			}
			
			if (index == 0)
				break;
			
			
			index--;
		}
		
		return found;
	}
}
