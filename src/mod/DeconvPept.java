package mod;

import java.util.ArrayList;
import java.util.List;

import dst.IsoCluster;
import dst.NextPair;
import ppk.Peak;
import utils.Constants;

public class DeconvPept {
	List<Peak> pkl;
	List<Peak> pks;
//	int peakNum;
//	int abundant;
//	int charge;
//	double mass;
//	double score;
//	double sum;
//	double avg;
	double lowIn;
	
	public DeconvPept() {
		pkl = null;
		pks = null;
	}
	
	public DeconvPept(List<Peak> pkl, List<Peak> pks) {
		this.pkl = pkl;
		this.pks = pks;
	}
	
	public List<IsoCluster> process(List<Peak> pkl, List<Peak> pks) {
		this.pkl = pkl;
		this.pks = pks;
		List<IsoCluster> result = new ArrayList<IsoCluster>(), clusts2 = new ArrayList<IsoCluster>();

		lowIn = 1000000;
//		System.out.println("Process start\t" + (System.currentTimeMillis()-Deconv.stime)/1000);
		for (int i = 0; i < pkl.size(); ++i) {
//		for (int i = 214; i < pkl.size(); ++i) {
//			System.out.println((i+1) + "th loof: start\t" + (System.currentTimeMillis()-Deconv.stime)/1000);
			if (pkl.get(i).in() > 0) {
				List<IsoCluster> clusts = new ArrayList<IsoCluster>();
				for (int z = 1; z <= Constants.maxCharge; ++z) {
					//if (getFirstTriples(clusts, i, z) == 0)
					clusts.addAll(getFirstTriples(i, z));
					clusts.addAll(getDouble(i, z));
//					System.out.println("Check Points");
				}
				clusts.addAll(getSingle(i, 1));
				clusts2.addAll(extendCluster(clusts, i));
			}
//			System.out.println((i+1) + "th loof: end\t" + (System.currentTimeMillis()-Deconv.stime)/1000);
		}
		
//		System.out.println("Recalculate score\t" + (System.currentTimeMillis()-Deconv.stime)/1000);
		clusts2 = recalcScore(clusts2);
//		System.out.println("Remove redundants\t" + (System.currentTimeMillis()-Deconv.stime)/1000);
		result = removeDuplicate(result, clusts2);

		
//		System.out.println("Process End\t" + (System.currentTimeMillis()-Deconv.stime)/1000);
		return result;
	}
	
//	public List<IsoCluster> process() {
//		List<IsoCluster> result = new ArrayList<IsoCluster>(), clusts2 = new ArrayList<IsoCluster>();
//
//		lowIn = 1000000;
//
//		for (int i = 0; i < plist.size(); ++i) {
//			if (plist.get(i).in() > 0) {
//				List<IsoCluster> clusts = new ArrayList<IsoCluster>();
//				for (int z = 1; z <= Constants.maxCharge; ++z) {
//					//if (getFirstTriples(clusts, i, z) == 0)
//					clusts = getFirstTriples(i, z);
//					clusts = getDouble(i, z);
//				}
//				clusts = getSingle(i, 1);
//				clusts2 = extendCluster(clusts, i);
//			}
//		}
//		
//		clusts2 = recalcScore(clusts2);
//		result = removeDuplicate(result, clusts2);
//
//		return result;
//	}
	
	void printScore(List<IsoCluster> clusts) {
		for( int i = 0 ; i < clusts.size(); i++ )
			System.out.println(clusts.get(i).score);
	}
	
	public List<IsoCluster> extendCluster(List<IsoCluster> triples, int curr) {
		List<IsoCluster> clusts = new ArrayList<IsoCluster>();
		int cnt = 0;
		int n, maxIdx;
		NextPair next_p = new NextPair();
		int shift = 0;
		double nextScore, maxScore;
		for (int i = 0; i<triples.size(); i++) {
			if (triples.get(i).plist.size() == 2) {
				clusts.add(triples.get(i));
				continue;
			}
			if (triples.get(i).plist.size() == 1) {		//by slee
				clusts.add(triples.get(i));
				continue;
			}

			IsoCluster c = triples.get(i);
			int z = c.charge;
			Peak p = c.plist.get(c.plist.size() - 1);
			getNextRange(next_p, curr,
				p.mz() + (Constants.MASS_ONE - Constants.massErr*p.mass(z)) / z,
				p.mz() + (Constants.MASS_ONE + Constants.massErr*p.mass(z)) / z);

			maxIdx = curr;

//			System.out.println(c.score);
			while (maxIdx < next_p.next_l && next_p.next_l <= next_p.next_r && c.score > Constants.thClustExt) {	// best next peak
				maxScore = Constants.NEGINF;
				n = c.plist.size();
//				System.out.println("Iso clu: " + n);
				for (int j = next_p.next_l; j <= next_p.next_r; ++j) {
					nextScore = evalNextTriple(c.peakNum + n - 2, c.mass, c.plist.get(n - 2), c.plist.get(n - 1), pkl.get(j));
					if (nextScore > maxScore) {
						maxScore = nextScore;
						maxIdx = j;
					}
				}
				if (maxScore < 0.0 && IntsRangeFunc.avgRatio(c.peakNum + n - 1, c.mass) <= 1.0) {
					IsoCluster tmp = new IsoCluster(); 
					tmp.score = c.score;
					tmp.peakNum = c.peakNum;
					tmp.abundant = c.abundant;
					tmp.avg = c.avg;
					tmp.charge = c.charge;
					tmp.mass = c.mass;
					tmp.plist = new ArrayList<>();
					tmp.plist.addAll(c.plist);
					tmp.sum = c.sum;
					
					Peak p1 = c.plist.get(0);
					Peak p2 = c.plist.get(n - 1);
					double prePeak, postPeak;
					if (Constants.truncated) {	// find intensity of smallest peak
						prePeak = c.plist.get(0).in();
						for (int k = 1; k<c.plist.size(); k++) {
							if (prePeak > c.plist.get(k).in()) prePeak = c.plist.get(k).in();
						}
						postPeak = prePeak = prePeak / 2;
					}
					else {	// find noise in raw data
						prePeak = highNoise(
							p1.mz() - (Constants.MASS_ONE + Constants.massErr*p1.mass(z)) / z,
							p1.mz() - (Constants.MASS_ONE - Constants.massErr*p1.mass(z)) / z);
						postPeak = highNoise(
							p2.mz() + (Constants.MASS_ONE - Constants.massErr*p2.mass(z)) / z,
							p2.mz() + (Constants.MASS_ONE + Constants.massErr*p2.mass(z)) / z);
					}
					if (c.peakNum > 0) {
						if (p1.in() / IntsRangeFunc.maxRatio(c.peakNum - 1, c.mass) > prePeak) {
							double x, l;
							x = p1.in() / prePeak - IntsRangeFunc.avgRatio(c.peakNum - 1, c.mass);
							l = IntsRangeFunc.maxRatio(c.peakNum - 1, c.mass) - IntsRangeFunc.avgRatio(c.peakNum - 1, c.mass);
							assert(l > 0);
							assert(x / l >= 0);
							tmp.score += 2.0 * (1.0 - x / l);
						}
					}
					if (p2.in()*IntsRangeFunc.minRatio(c.peakNum + n - 1, c.mass) > postPeak) {
						double x, l;
						x = postPeak / p2.in() - IntsRangeFunc.avgRatio(c.peakNum + n - 1, c.mass);
						l = IntsRangeFunc.minRatio(c.peakNum + n - 1, c.mass) - IntsRangeFunc.avgRatio(c.peakNum + n - 1, c.mass);
						assert(l < 0);
						assert(x / l >= 0);
						tmp.score += 2.0 * (1.0 - x / l);
					}
					//tmp.score /= tmp.mass;
					tmp.score /= -1.4633125034E-07*tmp.mass*tmp.mass + 3.4111225678E-03*tmp.mass + 7.4696036696E-01;
					
					if (tmp.score > Constants.thScore/*-0.1*/)
						clusts.add(tmp);
				}
				c.plist.add(pkl.get(maxIdx));
				c.score += maxScore;
				getNextRange(next_p, maxIdx,
					pkl.get(maxIdx).mz() + (Constants.MASS_ONE - Constants.massErr*pkl.get(maxIdx).mass(z)) / z,
					pkl.get(maxIdx).mz() + (Constants.MASS_ONE + Constants.massErr*pkl.get(maxIdx).mass(z)) / z);
			}
			n = c.plist.size();
			if (IntsRangeFunc.avgRatio(c.peakNum + n - 1, c.mass) <= 1.0) {
				Peak p1 = c.plist.get(0);
				Peak p2 = c.plist.get(n - 1);
				double prePeak, postPeak;
				if (Constants.truncated) {	// find intensity of smallest peak
					prePeak = c.plist.get(0).in();
					for (int k = 1; k<c.plist.size(); k++) {
						if (prePeak > c.plist.get(k).in()) prePeak = c.plist.get(k).in();
					}
					postPeak = prePeak = prePeak / 2;
				}
				else {	// find noise in raw data
					prePeak = highNoise(
						p1.mz() - (Constants.MASS_ONE + Constants.massErr*p1.mass(z)) / z,
						p1.mz() - (Constants.MASS_ONE - Constants.massErr*p1.mass(z)) / z);
					postPeak = highNoise(
						p2.mz() + (Constants.MASS_ONE - Constants.massErr*p2.mass(z)) / z,
						p2.mz() + (Constants.MASS_ONE + Constants.massErr*p2.mass(z)) / z);
				}
				if (c.peakNum > 0) {
					if (p1.in() / IntsRangeFunc.maxRatio(c.peakNum - 1, c.mass) > prePeak) {
						double x, l;
						x = p1.in() / prePeak - IntsRangeFunc.avgRatio(c.peakNum - 1, c.mass);
						l = IntsRangeFunc.maxRatio(c.peakNum - 1, c.mass) - IntsRangeFunc.avgRatio(c.peakNum - 1, c.mass);
						assert(l > 0);
						assert(x / l >= 0);
						c.score += 2.0 * (1.0 - x / l);
					}
				}
				if (p2.in()*IntsRangeFunc.minRatio(c.peakNum + n - 1, c.mass) > postPeak) {
					double x, l;
					x = postPeak / p2.in() - IntsRangeFunc.avgRatio(c.peakNum + n - 1, c.mass);
					l = IntsRangeFunc.minRatio(c.peakNum + n - 1, c.mass) - IntsRangeFunc.avgRatio(c.peakNum + n - 1, c.mass);
					assert(l < 0);
					assert(x / l >= 0);
					c.score += 2.0 * (1.0 - x / l);
				}
				//c.score /= c.mass;
				c.score /= -1.4633125034E-07*c.mass*c.mass + 3.4111225678E-03*c.mass + 7.4696036696E-01;
				if (c.score > Constants.thScore) {
					clusts.add(c);
				}
			}
		}

		return clusts;
	}
	
	private List<IsoCluster> getFirstTriples(int curr, int z) {
		List<IsoCluster> clusts = new ArrayList<IsoCluster>();
		int cnt1, cnt2;
		int j2;
		int j3;

		NextPair next_p2 = new NextPair(), next_p3 = new NextPair(); 
		
		cnt1 = cnt2 = 0;

		//find second peaks
		getNextRange(next_p2, curr,
			pkl.get(curr).mz() + (Constants.MASS_ONE - Constants.massErr*pkl.get(curr).mass(z)) / z,
			pkl.get(curr).mz() + (Constants.MASS_ONE + Constants.massErr*pkl.get(curr).mass(z)) / z);

		if (curr < next_p2.next_l && next_p2.next_l <= next_p2.next_r) {
			for (j2 = next_p2.next_l; j2 <= next_p2.next_r; ++j2) {

				//find third peaks
				getNextRange(next_p3, j2,
						pkl.get(j2).mz() + (Constants.MASS_ONE - Constants.massErr*pkl.get(j2).mass(z)) / z,
						pkl.get(j2).mz() + (Constants.MASS_ONE + Constants.massErr*pkl.get(j2).mass(z)) / z);
				if (j2 < next_p3.next_l && next_p3.next_l <= next_p3.next_r) {
					for (j3 = next_p3.next_l; j3 <= next_p3.next_r; ++j3) {
						//create a first triple
						++cnt1;

						//scoring each peaknum
						for (int k = 0; k <= Constants.maxMissPeak; ++k) {
							IsoCluster c = new IsoCluster();
							c.charge = z;
							c.plist.add(pkl.get(curr));
							c.plist.add(pkl.get(j2));
							c.plist.add(pkl.get(j3));
							c.peakNum = k;
							c.mass = c.getMass();
							c.score = evalTriple(k, c.mass, z, pkl.get(curr), pkl.get(j2), pkl.get(j3));
//							System.out.println(c.score);
							//miss peak이 존재할때 penalty를 부여 (by slee)
							if (k == 1)
								c.score = 0.8*c.score;
							else if (k == 2)
								c.score = 0.75*c.score;
							else if (k == 3)
								c.score = 0.7*c.score;

							if (c.score > Constants.thClustExt/* -2 */ && c.mass > Constants.thMassAbn1/* 300 */) {
								++cnt2;
								clusts.add(c);
								
								//by slee							
								if (pkl.get(curr).in() < lowIn)
									lowIn = pkl.get(curr).in();

							}
						}
					}
				}
			}
		}

		return clusts;
	}
	
	public List<IsoCluster> getDouble(int curr, int z) {
		List<IsoCluster> clusts = new ArrayList<IsoCluster>();
		int cnt;
		int j2;
		NextPair next_p2 = new NextPair();
		
		cnt = 0;
		//find second peaks
		getNextRange(next_p2, curr,
				pkl.get(curr).mz() + (Constants.MASS_ONE - Constants.massErr/*0.00001*/*pkl.get(curr).mass(z)) / z,
				pkl.get(curr).mz() + (Constants.MASS_ONE + Constants.massErr*pkl.get(curr).mass(z)) / z);


		if (curr < next_p2.next_l && next_p2.next_l <= next_p2.next_r) {
			for (j2 = next_p2.next_l; j2 <= next_p2.next_r; ++j2) {
				//create a double

				//scoring each peaknum
				for (int k = 0; k <= Constants.maxMissPeak; ++k) {
					IsoCluster c = new IsoCluster();
					c.charge = z;
					c.plist.add(pkl.get(curr));
					c.plist.add(pkl.get(j2));
					c.peakNum = k;
					c.mass = c.getMass();
					c.score = evalDouble(k, c.mass, z, pkl.get(curr), pkl.get(j2));

					//c.score /= c.mass;
					c.score /= -1.4633125034E-07*c.mass*c.mass + 3.4111225678E-03*c.mass + 7.4696036696E-01;

					if (c.score > Constants.thScore/*0*/ /*&& c.mass > this->ThMassAbn1/*300*/ && c.mass < 4000.0) {
						++cnt;
						clusts.add(c);
						
						//by slee
						if (pkl.get(curr).in() < lowIn)
							lowIn = pkl.get(curr).in();
					}
				}
			}
		}
		return clusts;
	}
	
	public List<IsoCluster> getSingle(int curr, int z) {
		List<IsoCluster> clusts = new ArrayList<IsoCluster>();
		int cnt = 0;
		NextPair next_p2 = new NextPair();
		//find second peaks
		getNextRange(next_p2, curr,
				pkl.get(curr).mz() + (Constants.MASS_ONE - Constants.massErr/*0.00001*/*pkl.get(curr).mass(z)) / z,
				pkl.get(curr).mz() + (Constants.MASS_ONE + Constants.massErr*pkl.get(curr).mass(z)) / z);

		if (0 < curr && curr < pkl.size() - 1) {
			if (next_p2.next_r == curr || curr == pkl.size() - 1) {		//현재 peak의 오른쪽으로 약1.01 da안에 다른 peak이 없을경우
																	//if(plist.get(curr).mz()-1.01 < _pkl[curr-1].mz())	//현재 peak의 왼쪽으로 1.01 da안에 다른 peak이 없을경우
																	//	return cnt;
				IsoCluster c = new IsoCluster();
				c.charge = z;
				c.plist.add(pkl.get(curr));
				c.peakNum = 0;
				c.mass = c.getMass();
				c.score = Constants.thScore + 0.001;
//				System.out.println(c.score);
				clusts.add(c);
				++cnt;


			}
		}
		return clusts;
	}
	
	//return 값이 -2보다 커야지 살아남음
	double evalTriple(int pknum, double mass, int charge, Peak pk1, Peak pk2, Peak pk3) {
		//miss peak이 1개이고 mass < 1839.60이면 버림
		//miss peak이 2개이고 mass < 3372.78이면 버림
		//miss peak이 3개이고 mass < 4975.32이면 버림	
		
		assert( pk1.mz() != pk2.mz() );
		assert( pk2.mz() != pk3.mz() );
		
		if (pknum > 0 && IntsRangeFunc.avgRatio(pknum - 1, mass) < 1.0) return Constants.NEGINF;
		
		double score = 0.0;
		double pr1 = pk2.in() / pk1.in();
		double pr2 = pk3.in() / pk2.in();
		double pr3 = pk1.in()*pk3.in() / pk2.in() / pk2.in();

//		System.out.println("pr1 "+ pr1);
//		System.out.println("pr2 "+ pr2);
//		System.out.println("pr3 "+ pr3);
		
		
		
		if (mass <= Constants.thMassAbn1/*300*/) {
			if (pr1 > pr2 && pr1 > pr3) return 1.0;
			return 0;
		}

		double x, l;
		x = pr1 - IntsRangeFunc.avgRatio(pknum, mass);
		if (x > 0) {
			l = IntsRangeFunc.maxRatio(pknum, mass) - IntsRangeFunc.avgRatio(pknum, mass);
			assert(l > 0);
		}
		else {
			l = IntsRangeFunc.minRatio(pknum, mass) - IntsRangeFunc.avgRatio(pknum, mass);
			assert(l < 0);
		}
		assert(x / l >= 0);
		score += (1.0 - x / l);
		x = pr2 - IntsRangeFunc.avgRatio(pknum + 1, mass);
		if (x > 0) {
			l = IntsRangeFunc.maxRatio(pknum + 1, mass) - IntsRangeFunc.avgRatio(pknum + 1, mass);
		}
		else {
			l = IntsRangeFunc.minRatio(pknum + 1, mass) - IntsRangeFunc.avgRatio(pknum + 1, mass);
		}
		assert(x / l >= 0);
		score += (1.0 - x / l);

		if (mass > Constants.thMassAbn2) {
			x = pr3 - IntsRangeFunc.avgRatioProduct(pknum, mass);
			if (x > 0) {
				l = IntsRangeFunc.maxRatioProduct(pknum, mass) - IntsRangeFunc.avgRatioProduct(pknum, mass);
				assert(l > 0);
			}
			else {
				l = IntsRangeFunc.minRatioProduct(pknum, mass) - IntsRangeFunc.avgRatioProduct(pknum, mass);
				assert(l < 0);
			}
			assert(x / l >= 0);
			score += (1.0 - x / l);
		}
		
		return score;
	}
	
	double evalDouble(int pknum, double mass, int charge, Peak pk1, Peak pk2) {
		double score = 0.0;
		double x, l;

		if (pknum > 0 && IntsRangeFunc.avgRatio(pknum - 1, mass) < 1.0) return Constants.NEGINF;
		if (IntsRangeFunc.avgRatio(pknum + 1, mass) > 1.0) return Constants.NEGINF;

		x = pk2.in() / pk1.in() - IntsRangeFunc.avgRatio(pknum, mass);
		if (x > 0) {
			l = IntsRangeFunc.maxRatio(pknum, mass) - IntsRangeFunc.avgRatio(pknum, mass);
		}
		else {
			l = IntsRangeFunc.minRatio(pknum, mass) - IntsRangeFunc.avgRatio(pknum, mass);
		}
		assert(x / l >= 0);
		score += (1.0 - x / l);

		double prePeak, postPeak;
		if (Constants.truncated) {	// find intensity of smallest peak
			postPeak = prePeak = (pk1.in() < pk2.in()) ? pk1.in() / 2 : pk2.in() / 2;
		}
		else {	// find noise in raw data
			prePeak = highNoise(
				pk1.mz() - (Constants.MASS_ONE + Constants.massErr*pk1.mass(charge)) / charge,
				pk1.mz() - (Constants.MASS_ONE - Constants.massErr*pk1.mass(charge)) / charge);
			postPeak = highNoise(
				pk2.mz() + (Constants.MASS_ONE - Constants.massErr*pk2.mass(charge)) / charge,
				pk2.mz() + (Constants.MASS_ONE + Constants.massErr*pk2.mass(charge)) / charge);
		}

		if (pknum > 0) {
			if (pk1.in() / IntsRangeFunc.maxRatio(pknum - 1, mass) > prePeak) {
				x = pk1.in() / prePeak - IntsRangeFunc.avgRatio(pknum - 1, mass);
				l = IntsRangeFunc.maxRatio(pknum - 1, mass) - IntsRangeFunc.avgRatio(pknum - 1, mass);
				assert(l > 0);
				assert(x / l >= 0);
				score += 2.0 * (1.0 - x / l);
			}
		}
		if (pk2.in()*IntsRangeFunc.minRatio(pknum + 1, mass) > postPeak) {
			x = postPeak / pk2.in() - IntsRangeFunc.avgRatio(pknum + 1, mass);
			l = IntsRangeFunc.minRatio(pknum + 1, mass) - IntsRangeFunc.avgRatio(pknum + 1, mass);
			assert(l < 0);
			assert(x / l >= 0);
			score += 2.0 * (1.0 - x / l);
		}

		return score;
	}
	
	double evalNextTriple(int pknum, double mass, Peak pk1, Peak pk2, Peak pk3) {
		double score = 0.0;
		double x, l;

		x = pk3.in() / pk2.in() - IntsRangeFunc.avgRatio(pknum + 1, mass);
		if (x > 0) {
			l = IntsRangeFunc.maxRatio(pknum + 1, mass) - IntsRangeFunc.avgRatio(pknum + 1, mass);
		}
		else {
			l = IntsRangeFunc.minRatio(pknum + 1, mass) - IntsRangeFunc.avgRatio(pknum + 1, mass);
		}
		assert(x / l >= 0);
		score += (1.0 - x / l);
		if (mass > Constants.thMassAbn2) {
			x = pk1.in()*pk3.in() / pk2.in() / pk2.in() - IntsRangeFunc.avgRatioProduct(pknum, mass);
			if (x > 0) {
				l = IntsRangeFunc.maxRatioProduct(pknum, mass) - IntsRangeFunc.avgRatioProduct(pknum, mass);
				assert(l >= 0);
			}
			else {
				l = IntsRangeFunc.minRatioProduct(pknum, mass) - IntsRangeFunc.avgRatioProduct(pknum, mass);
				assert(l <= 0);
			}
			assert(x / l >= 0);
			score += (1.0 - x / l);
		}
		return score;
	}
	
	List<IsoCluster> recalcScore(List<IsoCluster> clusts) {
		//	double maxI, minI, minI2;
		//	double x, l;
		//	double bias;
		double tempSum;
		for (int i = 0; i < clusts.size(); ++i) {
			clusts.get(i).abundant = 0;
			clusts.get(i).sum = 0.0;
			tempSum = 0.0;
			for (int j = 0; j < clusts.get(i).plist.size(); ++j) {
				if (clusts.get(i).plist.get(clusts.get(i).abundant).in() < clusts.get(i).plist.get(j).in()) clusts.get(i).abundant = j;
			}
			for (int j = 0; j < clusts.get(i).plist.size(); ++j) {
				clusts.get(i).sum += clusts.get(i).plist.get(j).in();
				tempSum += clusts.get(i).plist.get(j).in() * clusts.get(i).plist.get(j).mass(clusts.get(i).charge);
			}
			clusts.get(i).avg = tempSum / clusts.get(i).sum;
			if (clusts.get(i).plist.size() > Constants.maxAbundancePeak) {
				int i_s, i_e;
				i_s = i_e = clusts.get(i).abundant;
				for (int j = 0; j < Constants.maxAbundancePeak - 1; ++j) {
					if (i_s == 0) i_e++;
					else if (i_e == clusts.get(i).plist.size() - 1) i_s--;
					else if (clusts.get(i).plist.get(i_s - 1).in() > clusts.get(i).plist.get(i_e + 1).in()) i_s--;
					else i_e++;
				}
				clusts.get(i).sum = 0.0;
				for (int j = i_s; j <= i_e; ++j) {
					clusts.get(i).sum += clusts.get(i).plist.get(j).in();
				}
			}

		}
		return clusts;
	}
	
	List<IsoCluster> removeDuplicate(List<IsoCluster> result, List<IsoCluster> clusts) {
		int cnt = 0;
		int i, j;

		for (i = 0; i < clusts.size(); ++i) {
			//if( clusts.get(i).score < ThScore ) continue;		//extend에서 수행함.
			for (j = 0; j < clusts.size(); ++j) {
				if (clusts.get(j).plist.get(clusts.get(j).plist.size() - 1).mz() < clusts.get(i).plist.get(0).mz()) continue;
				if (clusts.get(i).plist.get(clusts.get(i).plist.size() - 1).mz() < clusts.get(j).plist.get(0).mz()) {	//clust가 안겹치면 넣어줌
					result.add(clusts.get(i));
					++cnt;
					break;
				}
				
				//clust i, j 가 겹칠때

				/*	if(clusts.get(i).List[0].mz() > 897.22 && clusts.get(i).List[0].mz() < 897.24 && clusts.get(i).charge == 4 && clusts.get(j).charge == 4){
				printf("i = %lf, %lf, %d, %d\n",clusts.get(i).Mass, clusts.get(i).score, clusts.get(i).charge, clusts.get(i).PeakNum);
				printf("j = %lf, %lf, %d, %d\n\n",clusts.get(j).Mass, clusts.get(j).score, clusts.get(j).charge, clusts.get(j).PeakNum);
				}*/

				if (clusts.get(i).plist.get(clusts.get(i).abundant).in() > clusts.get(j).plist.get(clusts.get(j).abundant).in() //abundant의 intensity가 더 클때
					|| (clusts.get(i).charge % clusts.get(j).charge != 0 && clusts.get(j).charge % clusts.get(i).charge != 0)) continue;


				//monoisotopic의 intensity가 더 작고 두 클러스트의 charge가 배수관계일때
				if (isIncluded(clusts.get(i).plist, clusts.get(j).plist) /*isIncluded(clusts.get(i).List, clusts.get(j).List[clusts.get(j).abundant].mz())*/) {
					if (clusts.get(i).plist.get(clusts.get(i).abundant).in() < clusts.get(j).plist.get(clusts.get(j).abundant).in()) break;
					if (clusts.get(i).charge < clusts.get(j).charge) break;
					if (clusts.get(i).charge == clusts.get(j).charge && clusts.get(i).score < clusts.get(j).score) break;
				}
			}
			if (j == clusts.size()) {
				result.add(clusts.get(i));
				++cnt;
			}
		}

		return result;
	}

	boolean isIncluded(List<Peak> list, double mz) {
		for (int i = 0; i < list.size(); ++i) {
			if (list.get(i).mz() == mz) return true;
		}
		return false;
	}

	boolean isIncluded(List<Peak> list1, List<Peak> list2) {
		int j = 0;
		//printf("\t\t\tlist size = %d, %d\n", list1.size(), list2.size() );
		for (int i = 0; i < list1.size(); ++i) {
			j = 0;//이거없으면....안될듯한데???
			while (list1.get(i).mz() > list2.get(j).mz() && j < list2.size() - 1) j++;
			if (list1.get(i).mz() == list2.get(j).mz()) return true;
		}
		return false;
	}

	double highNoise(double start, double end) {
		int s = 0, e = pks.size() - 1;
		while (s < e) {
			if (pks.get((s + e) / 2).mdbl_mz == start) s = e = (s + e) / 2;
			else if (pks.get((s + e) / 2).mdbl_mz < start) s = (s + e) / 2 + 1;
			else e = (s + e) / 2 - 1;
		}
		if (pks.get(s).mdbl_mz > end) return Constants.backgroundIntensity/*0*/;
		double max = pks.get(s).mdbl_intensity;
		while (pks.get(s).mdbl_mz < end && s < pks.size() - 1) {
			s++;
			if (pks.get(s).mdbl_intensity > max) max = pks.get(s).mdbl_intensity;
		}
		return max;
	}
	
	void getNextRange(NextPair next_p, int curr, double mz_l, double mz_r) {
		int i, j;
		next_p.next_l = -1;	//범위안에 드는 첫번째 peak
		next_p.next_r = -1;	//범위안에 드는 마지막 peak

		if (pkl.get(curr).mz() < mz_l) {
			for (i = curr + 1; i < pkl.size(); ++i) {
				if (pkl.get(i).mz() >= mz_l) break;
			}
			if (i < pkl.size()) {
				next_p.next_l = i;
				for (j = i; j < pkl.size(); ++j) {
					if (pkl.get(j).mz() > mz_r) break;
				}
				next_p.next_r = j - 1;
			}

		}
		else if (mz_r < pkl.get(curr).mz()) {
			for (i = curr - 1; i >= 0; --i) {
				if (pkl.get(i).mz() <= mz_r) break;
			}
			if (i >= 0) {
				next_p.next_r = i;
				for (j = i; j >= 0; --j) {
					if (pkl.get(j).mz() < mz_l) break;
				}
				next_p.next_l = j + 1;
			}
		}
		else {
			for (i = curr - 1; i >= 0; --i) {
				if (pkl.get(i).mz() < mz_l) break;
			}
			if (i >= 0) {
				next_p.next_l = i + 1;
				for (j = i; j < pkl.size(); ++j) {
					if (pkl.get(j).mz() > mz_r) break;
				}
				next_p.next_r = j - 1;
			}
		}

	} 
}
