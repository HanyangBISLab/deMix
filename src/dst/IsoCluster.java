package dst;

import java.util.ArrayList;
import java.util.List;

import ppk.Peak;
import utils.Constants;

public class IsoCluster implements Comparable<IsoCluster>{
	
	public int peakNum;
	public int abundant;
	public int charge;
	public double mass;
	public double score;
	public double sum;
	public double avg;
	
	public List<Peak> plist;
	public IsoCluster() {
		plist = new ArrayList<Peak>();
	}
	
	static double AveragineMonoMass = 111.05429999675;
	static double AveragineAvgMass = 111.1257144233;
	static double Averagine[] = { 4.9384, 7.7583, 1.3577, 1.4773, 0.0417 };
	static double AtomMonoMass[] = {
		12.0,
		1.0078246,
		14.0030732,
		15.9949141,
		31.972070
	};

	static double AtomAvgMass[] = {
		12.011107,
		1.007976,
		14.006724,
		15.999370,
		32.064387
	};
	

	public double getMass() {
		double ints = plist.get(0).mdbl_intensity;
		int idx = 0;
		for (int i = 1; i < plist.size(); ++i) {
			if (plist.get(i).mdbl_intensity > ints) {
				ints = plist.get(i).mdbl_intensity;
				idx = i;
			}
		}

		double tmp = plist.get(0).mass(charge) - peakNum*Constants.MASS_ONE;
		return plist.get(idx).mass(charge) - diffMass(peakNum + idx, tmp);
	}

	/*double CIsoCluster::getMassSulfur(int numSulfur) const {
	double tmp = List[0].mass(Charge) - PeakNum*MassOne;
	return List[this->Abundant].mass(Charge) - diffMassSulfur(PeakNum + this->Abundant, tmp, numSulfur);
	}*/

	public double getMassAveragine() {
		double m = getMass();

		int nC = (int)(m / AveragineMonoMass * Averagine[0] + 0.5);
		int nN = (int)(m / AveragineMonoMass * Averagine[2] + 0.5);
		int nO = (int)(m / AveragineMonoMass * Averagine[3] + 0.5);
		int nS = (int)(m / AveragineMonoMass * Averagine[4] + 0.5);

		double mono_mw = AtomMonoMass[0] * nC + AtomMonoMass[2] * nN + AtomMonoMass[3] * nO + AtomMonoMass[4] * nS;
		int nH = (int)((m - mono_mw) / AtomMonoMass[1] + 0.5);
		mono_mw += nH*AtomMonoMass[1];

		return mono_mw;
	}


	double diffMass(int pknum, double mass) {
		switch (pknum) {
		case 0:
			return 0;
		case 1:
			return 1.002858;
		case 2:
			return (3.636027E-3 + 110.998 / mass*9.845597E-3) / (1.812832E-3 + 110.998 / mass*4.920315E-3);
		case 3:
			return (8.026887E-5 + 110.998 / mass*8.899539E-4) / (2.667061E-5 + 110.998 / mass*2.962693E-4);
		default:
			return pknum*Constants.MASS_ONE;
		}
	}

	@Override
	public int compareTo(IsoCluster iso) {
		if (Constants.ResultOrder.compareTo("M/Z") == 0) {
			// if (a.plist.get(a.abundant).mz() < b.plist.get(b.abundant).mz()) return 1;
			// else return -1;
			if (plist.get(abundant).mz() - iso.plist.get(iso.abundant).mz() < 0)
				return 1;
			else if (plist.get(abundant).mz() - iso.plist.get(iso.abundant).mz() > 0)
				return -1;
			return 0;
		} else if (Constants.ResultOrder.compareTo("SCORE") == 0) {
			if (score > iso.score)
				return 1;
			else if (score < iso.score)
				return -1;
			return 0;
		} else if (Constants.ResultOrder.compareTo("MASS") == 0) {
			if (mass < iso.mass)
				return 1;
			else if (mass > iso.mass)
				return -1;
			return 0;
		}
		if (sum < iso.sum)
			return 1;
		else if (sum > iso.sum)
			return -1;
		return 0;
	}
}
