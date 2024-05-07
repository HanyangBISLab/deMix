package mod;

import utils.Constants;

public class IntsRangeFunc {
	public static double FctRange = 0.5;

	// 비율, 비율곱 범위 설정 배율, default 1.0
	private static double _r_int;

	// 질량 mass에 따른 이상적인 idx번째 비율((idx+1)/idx) 리턴
	public static double avgRatio(int idx, double mass) {
		if (mass <= Constants.T_MASS2/*1800*/) {
			return avgR1(idx, mass);
		}
		else if (mass <= Constants.T_MASS3 /*4000*/) {
			return avgR2(idx, mass);
		}
		return avgR3(idx, mass);
	}
	
	public static double maxRatio(int idx, double mass) {
		double avg, max = 0;

		if (mass <= Constants.T_MASS2) {
			avg = avgR1(idx, mass);
			max = maxR1(idx, mass);
		}
		else if (mass <= Constants.T_MASS3) {
			avg = avgR2(idx, mass);
			max = maxR2(idx, mass);
		}
		else {
			avg = avgR3(idx, mass);
			max = maxR3(idx, mass);
		}

		if (max > 1.0) return max / (1.0 - FctRange);
		else return max * (1.0 + FctRange);
	}
	
	public static double minRatio(int idx, double mass) {
		double avg, min = 999999;

		if (mass <= Constants.T_MASS2) {
			avg = avgR1(idx, mass);
			min = minR1(idx, mass);
		}
		else if (mass <= Constants.T_MASS3) {
			avg = avgR2(idx, mass);
			min = minR2(idx, mass);
		}
		else {
			avg = avgR3(idx, mass);
			min = minR3(idx, mass);
		}

		if (min < 1.0) return min * (1.0 - FctRange);
		else return min / (1.0 + FctRange);
	}

	// 질량 mass에 따른 이상적인 idx번째 비율곱(idx idx+2 / idx+1 ^2) 리턴
	public static double avgRatioProduct(int idx, double mass) {
		if (mass <= Constants.T_MASS2) {
			return avgRP1(idx, mass);
		}
		else if (mass <= Constants.T_MASS3) {
			return avgRP2(idx, mass);
		}
		return avgRP3(idx, mass);
	}
	
	public static double maxRatioProduct(int idx, double mass) {
		double avg, max = 0;

		if (mass <= Constants.T_MASS2) {
			avg = avgRP1(idx, mass);
			max = maxRP1(idx, mass);
		}
		else if (mass <= Constants.T_MASS3) {
			avg = avgRP2(idx, mass);
			max = maxRP2(idx, mass);
		}
		else {
			avg = avgRP3(idx, mass);
			max = maxRP3(idx, mass);
		}
		return max * (1.0 + FctRange);
	}

	public static double minRatioProduct(int idx, double mass) {
		double avg, min = 999999;

		if (mass <= Constants.T_MASS2) {
			avg = avgRP1(idx, mass);
			min = minRP1(idx, mass);
		}
		else if (mass <= Constants.T_MASS3) {
			avg = avgRP2(idx, mass);
			min = minRP2(idx, mass);
		}
		else {
			avg = avgRP3(idx, mass);
			min = minRP3(idx, mass);
		}

		return min * (1.0 - FctRange);
	}


	public static double avgRS(int idx, double mass, int sulfur) {
		if (mass < Constants.T_MASS1) {
			return 0;
		}
		else if (mass < 1000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return 5.1156765634E-04*mass - 1.4884062127E-02;
				case 1:
					return 1.4056685185E-07*mass*mass + 3.2497226322E-05*mass + 1.4910215008E-01;
				case 2:
					return 1.9600517826E-04*mass + 4.2625179836E-02;
				case 3:
					return 3.3436122859E-08*mass*mass + 1.0163352381E-04*mass + 6.2206027491E-02;
				case 4:
					return 1.3182588545E-04*mass + 3.6488072794E-02;
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return 5.1390100036E-04*mass - 2.3888543901E-02;
				case 1:
					return 7.3391359384E-12*mass*mass*mass*mass - 2.2864148498E-08*mass*mass*mass + 2.6577770588E-05*mass*mass - 1.3565883937E-02*mass + 2.9214291833E+00;
				case 2:
					return -1.6676364717E-07*mass*mass + 4.6333949123E-04*mass - 1.0807247825E-02;
				case 3:
					return 2.5143592613E-12*mass*mass*mass*mass - 7.6805574353E-09*mass*mass*mass + 8.6629493239E-06*mass*mass - 4.1321717680E-03*mass + 8.7715548566E-01;
				case 4:
					return -4.0992269673E-08*mass*mass + 2.0100080329E-04*mass + 4.1766398246E-02;
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return 5.0661088174E-04*mass - 2.8630182526E-02;
				case 1:
					return 9.3025607596E-12*mass*mass*mass*mass - 2.9531742070E-08*mass*mass*mass + 3.5163215151E-05*mass*mass - 1.8690132667E-02*mass + 4.2543402487E+00;
				case 2:
					return -2.3098689991E-07*mass*mass + 6.0393043472E-04*mass - 5.7206042697E-02;
				case 3:
					return 4.3952342834E-12*mass*mass*mass*mass - 1.3763988053E-08*mass*mass*mass + 1.6080281849E-05*mass*mass - 8.2349888797E-03*mass + 1.8106984267E+00;
				case 4:
					return -1.4217866407E-07*mass*mass + 3.6766379356E-04*mass + 6.4286826385E-03;
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return 5.0620701896E-04*mass - 3.9991518444E-02;
				case 1:
					return 1.3421847449E-06*mass*mass - 2.4666492623E-03*mass + 1.7262391637E+00;
				case 2:
					return 2.9608194958E-04*mass + 4.8269670271E-02;
				case 3:
					return 4.5425863870E-07*mass*mass - 7.5391715107E-04*mass + 6.3554738808E-01;
				case 4:
					return 1.8335603608E-04*mass + 8.1081143039E-02;
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return 4.8458368554E-04*mass - 3.4280504615E-02;
				case 1:
					return 2.0062789369E-06*mass*mass - 3.7390941798E-03*mass + 2.4457835598E+00;
				case 2:
					return -2.9006888789E-07*mass*mass + 7.6808888000E-04*mass - 1.3458175242E-01;
				case 3:
					return -1.0990391745E-09*mass*mass*mass + 3.2972465742E-06*mass*mass - 3.2645900272E-03*mass + 1.4498785139E+00;
				case 4:
					return -2.5882492332E-07*mass*mass + 6.1699277860E-04*mass - 8.3291360158E-02;
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 2000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return 4.8636049445E-04*mass + 6.5925583258E-03;
				case 1:
					return 2.4705896250E-04*mass + 7.1387486256E-02;
				case 2:
					return 1.7545106571E-04*mass + 6.3522431440E-02;
				case 3:
					return 1.3754016982E-04*mass + 5.9523655275E-02;
				case 4:
					return 1.1514168197E-04*mass + 5.4141147048E-02;
				case 5:
					return 9.9892462046E-05*mass + 4.9762466931E-02;
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return 4.8761736394E-04*mass - 1.7189959046E-03;
				case 1:
					return 2.0442374857E-04*mass + 1.9378110714E-01;
				case 2:
					return 1.6225987879E-04*mass + 1.2809140182E-01;
				case 3:
					return 1.2767456886E-04*mass + 1.1329847671E-01;
				case 4:
					return 1.1076184938E-04*mass + 9.2947165923E-02;
				case 5:
					return 9.7480166217E-05*mass + 8.0841020675E-02;
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return 4.8906818302E-04*mass - 1.3408284527E-02;
				case 1:
					return 6.0305525631E-08*mass*mass - 2.2486768311E-05*mass + 4.5337328187E-01;
				case 2:
					return 1.6339067027E-04*mass + 1.6009513460E-01;
				case 3:
					return 1.1477283348E-04*mass + 1.6864735977E-01;
				case 4:
					return 1.0569832103E-04*mass + 1.3046619845E-01;
				case 5:
					return 9.1713379538E-05*mass + 1.1681458069E-01;
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return 4.9388920732E-04*mass - 3.3587670981E-02;
				case 1:
					return 1.0353808544E-07*mass*mass - 2.0417169227E-04*mass + 6.8725372519E-01;
				case 2:
					return 1.7314677588E-04*mass + 1.7040335618E-01;
				case 3:
					return 2.2395404246E-08*mass*mass + 3.0380451528E-05*mass + 2.7766938610E-01;
				case 4:
					return 1.0498082268E-04*mass + 1.5729205598E-01;
				case 5:
					return 8.3683167117E-05*mass + 1.5496646932E-01;
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return 4.9044371317E-04*mass - 3.9980329757E-02;
				case 1:
					return 7.2412086340E-14*mass*mass*mass*mass - 5.2477339171E-10*mass*mass*mass + 1.5085495688E-06*mass*mass - 1.8946245510E-03*mass + 1.5280724165E+00;
				case 2:
					return -4.8710968276E-08*mass*mass + 3.3123386552E-04*mass + 6.6265751767E-02;
				case 3:
					return 3.9296969013E-08*mass*mass - 4.0332770224E-05*mass + 3.7884024511E-01;
				case 4:
					return -2.4727066344E-08*mass*mass + 1.8235372228E-04*mass + 1.2103085865E-01;
				case 5:
					return 9.7217576500E-09*mass*mass + 4.4864072385E-05*mass + 2.1591049125E-01;
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 3000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return 4.8696202100E-04*mass + 6.8728674948E-03;
				case 1:
					return 2.4529567117E-04*mass + 7.5066742082E-02;
				case 2:
					return 1.6832553015E-04*mass + 7.7304134617E-02;
				case 3:
					return 1.2985959856E-04*mass + 7.4364093600E-02;
				case 4:
					return 1.0686059586E-04*mass + 7.0106940351E-02;
				case 5:
					return 9.1495059706E-05*mass + 6.5946800735E-02;
				case 6:
					return 8.0492042604E-05*mass + 6.2086263489E-02;
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return 4.9029561356E-04*mass - 7.1179539033E-03;
				case 1:
					return 2.3170983138E-04*mass + 1.4201784154E-01;
				case 2:
					return 1.5919027123E-04*mass + 1.3334281193E-01;
				case 3:
					return 1.2306245905E-04*mass + 1.2201063239E-01;
				case 4:
					return 1.0212366757E-04*mass + 1.0938793108E-01;
				case 5:
					return 8.8180322232E-05*mass + 9.8707544419E-02;
				case 6:
					return 7.8214260426E-05*mass + 8.9642318935E-02;
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return 4.9285343234E-04*mass - 2.0902430858E-02;
				case 1:
					return 2.1714852805E-04*mass + 2.1130341987E-01;
				case 2:
					return 1.5309700870E-04*mass + 1.7823322825E-01;
				case 3:
					return 1.1758014701E-04*mass + 1.6347038246E-01;
				case 4:
					return 9.8022777678E-05*mass + 1.4478678581E-01;
				case 5:
					return 8.4811543354E-05*mass + 1.3003675161E-01;
				case 6:
					return 7.5541848185E-05*mass + 1.1718415954E-01;
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return 4.9450341734E-04*mass - 3.6359993354E-02;
				case 1:
					return 2.0131719952E-04*mass + 2.8340619120E-01;
				case 2:
					return 1.4962778758E-04*mass + 2.1305749105E-01;
				case 3:
					return 1.1260081908E-04*mass + 2.0141816230E-01;
				case 4:
					return 9.4732403840E-05*mass + 1.7609332040E-01;
				case 5:
					return 8.1746997900E-05*mass + 1.5907924525E-01;
				case 6:
					return 7.3025630363E-05*mass + 1.4318223409E-01;
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return 4.9128805641E-04*mass - 4.3499321014E-02;
				case 1:
					return 1.8312973237E-04*mass + 3.6138043907E-01;
				case 2:
					return 1.4751434023E-04*mass + 2.4117590941E-01;
				case 3:
					return 1.0703413321E-04*mass + 2.3910967697E-01;
				case 4:
					return 9.1837045305E-05*mass + 2.0449558674E-01;
				case 5:
					return 7.8543575758E-05*mass + 1.8707421061E-01;
				case 6:
					return 7.0522686469E-05*mass + 1.6794373383E-01;
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 4000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return 4.8536852145E-04*mass + 1.3467724918E-02;
				case 1:
					return 2.4370454845E-04*mass + 8.0514700408E-02;
				case 2:
					return 1.6522394719E-04*mass + 8.6887584818E-02;
				case 3:
					return 1.2630357036E-04*mass + 8.5157943750E-02;
				case 4:
					return 1.0306650405E-04*mass + 8.1529745824E-02;
				case 5:
					return 8.7604946295E-05*mass + 7.7607737969E-02;
				case 6:
					return 7.6564976298E-05*mass + 7.3825472958E-02;
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return 4.8855036664E-04*mass - 3.0840322828E-04;
				case 1:
					return 2.3786123432E-04*mass + 1.2454024987E-01;
				case 2:
					return 1.5995347195E-04*mass + 1.3132119818E-01;
				case 3:
					return 1.2196385299E-04*mass + 1.2538292455E-01;
				case 4:
					return 9.9646844884E-05*mass + 1.1676329290E-01;
				case 5:
					return 8.4949634563E-05*mass + 1.0828351903E-01;
				case 6:
					return 7.4525809781E-05*mass + 1.0055636577E-01;
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return 4.9660580738E-04*mass - 3.4350765833E-02;
				case 1:
					return 2.3325403960E-04*mass + 1.6333035139E-01;
				case 2:
					return 1.5640458206E-04*mass + 1.6787516280E-01;
				case 3:
					return 1.1888773477E-04*mass + 1.5932986829E-01;
				case 4:
					return 9.7094583858E-05*mass + 1.4731739650E-01;
				case 5:
					return 8.2801603960E-05*mass + 1.3584549614E-01;
				case 6:
					return 7.2721826583E-05*mass + 1.2543066696E-01;
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return 4.9553403780E-04*mass - 3.9529472313E-02;
				case 1:
					return 2.2502082088E-04*mass + 2.1387670691E-01;
				case 2:
					return 1.5189011821E-04*mass + 2.0590876626E-01;
				case 3:
					return 1.1520663546E-04*mass + 1.9366449588E-01;
				case 4:
					return 9.4201632163E-05*mass + 1.7750824743E-01;
				case 5:
					return 8.0388549655E-05*mass + 1.6301297621E-01;
				case 6:
					return 7.0691861386E-05*mass + 1.4999601515E-01;
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return 4.9920963756E-04*mass - 6.5530361473E-02;
				case 1:
					return 2.1781985599E-04*mass + 2.5996622405E-01;
				case 2:
					return 1.4906200600E-04*mass + 2.3623480900E-01;
				case 3:
					return 1.1240803180E-04*mass + 2.2348705869E-01;
				case 4:
					return 9.2039884788E-05*mass + 2.0385620324E-01;
				case 5:
					return 7.8470383438E-05*mass + 1.8737663797E-01;
				case 6:
					return 6.9021534353E-05*mass + 1.7240555976E-01;
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		return 0.0;
	}
	
	public static double maxRS(int idx, double mass, int sulfur) {
		if (mass < Constants.T_MASS1) {
			return Constants.POSINF;
		}
		else if (mass < 1000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return 6.2005168323E-04*mass - 3.3350589217E-02;
				case 1:
					return 1.2548519225E-07*mass*mass + 2.0863921013E-05*mass + 1.8640479082E-01;
				case 2:
					return 2.0691981297E-04*mass + 4.3016497954E-02;
				case 3:
					return 1.5699720481E-08*mass*mass + 1.1748196221E-04*mass + 7.2990669481E-02;
				case 4:
					return 1.3722908825E-04*mass + 4.1736713121E-02;
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return 6.3869794624E-04*mass - 6.3273442624E-02;
				case 1:
					return 1.3947731479E-12*mass*mass*mass*mass - 5.1628291692E-09*mass*mass*mass + 7.2421716200E-06*mass*mass - 4.4623958410E-03*mass + 1.4103646625E+00;
				case 2:
					return -2.3930629795E-07*mass*mass + 5.7568364456E-04*mass - 4.0446074237E-02;
				case 3:
					return 1.2473281569E-12*mass*mass*mass*mass - 4.1758608009E-09*mass*mass*mass + 5.1702521843E-06*mass*mass - 2.6799375806E-03*mass + 6.8536695981E-01;
				case 4:
					return -6.9340619161E-08*mass*mass + 2.4868001125E-04*mass + 2.8313988711E-02;
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return 6.1368701057E-04*mass - 6.2491828828E-02;
				case 1:
					return 2.1949794107E-12*mass*mass*mass*mass - 7.5518195880E-09*mass*mass*mass + 1.0268994044E-05*mass*mass - 6.5662465128E-03*mass + 2.1836683401E+00;
				case 2:
					return -3.2012784410E-07*mass*mass + 7.4770765148E-04*mass - 9.5944548394E-02;
				case 3:
					return 1.7033984779E-12*mass*mass*mass*mass - 5.5833214024E-09*mass*mass*mass + 6.9928272087E-06*mass*mass - 3.9140904538E-03*mass + 1.0979291981E+00;
				case 4:
					return -1.5201222693E-07*mass*mass + 3.8158259264E-04*mass + 7.7573174564E-03;
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return 5.9776436170E-04*mass - 7.3585821809E-02;
				case 1:
					return 9.6322660424E-07*mass*mass - 1.9247813160E-03*mass + 1.6039278137E+00;
				case 2:
					return 3.1774632285E-04*mass + 4.8587445999E-02;
				case 3:
					return 4.1125522318E-07*mass*mass - 7.0780320423E-04*mass + 6.4765504906E-01;
				case 4:
					return 1.8439929463E-04*mass + 8.8014369392E-02;
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return 5.2143187040E-04*mass - 4.1324822846E-02;
				case 1:
					return 2.2215876313E-06*mass*mass - 4.0707064490E-03*mass + 2.6075531578E+00;
				case 2:
					return -2.2221896028E-07*mass*mass + 6.8781216007E-04*mass - 1.0459431282E-01;
				case 3:
					return -1.4071715684E-09*mass*mass*mass + 4.1675929813E-06*mass*mass - 4.0659123555E-03*mass + 1.7032205870E+00;
				case 4:
					return -2.1394605579E-07*mass*mass + 5.4902323513E-04*mass - 5.2353851881E-02;
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 2000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return 5.1453228503E-04*mass + 5.3429942457E-02;
				case 1:
					return 2.6420140714E-04*mass + 6.9793999289E-02;
				case 2:
					return 1.8276993579E-04*mass + 6.5152556310E-02;
				case 3:
					return 1.3896974437E-04*mass + 6.6077573438E-02;
				case 4:
					return 1.1372192619E-04*mass + 6.3722090711E-02;
				case 5:
					return 9.6875024902E-05*mass + 6.1587892646E-02;
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return 5.1763279388E-04*mass + 3.5567719181E-02;
				case 1:
					return 2.2138223102E-04*mass + 1.8589267347E-01;
				case 2:
					return 1.6620687549E-04*mass + 1.3518610677E-01;
				case 3:
					return 1.3110075368E-04*mass + 1.1601500949E-01;
				case 4:
					return 1.1237588959E-04*mass + 9.6825805617E-02;
				case 5:
					return 9.7115611761E-05*mass + 8.7402292358E-02;
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return 5.3498749415E-04*mass - 1.7289912241E-03;
				case 1:
					return 1.1894330550E-07*mass*mass - 1.9797750364E-04*mass + 5.9644828717E-01;
				case 2:
					return 1.6638114131E-04*mass + 1.7027710803E-01;
				case 3:
					return 1.1739218362E-04*mass + 1.7220566457E-01;
				case 4:
					return 1.0809757456E-04*mass + 1.3282359816E-01;
				case 5:
					return 9.1783217522E-05*mass + 1.2196835372E-01;
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return 5.5287670687E-04*mass - 4.4728500306E-02;
				case 1:
					return 1.7127425173E-07*mass*mass - 4.2275926946E-04*mass + 8.7488608075E-01;
				case 2:
					return 1.7531949895E-04*mass + 1.8395896158E-01;
				case 3:
					return 4.6014891645E-08*mass*mass - 3.9788049673E-05*mass + 3.3521653746E-01;
				case 4:
					return 1.0738247285E-04*mass + 1.6008826073E-01;
				case 5:
					return 8.4579691569E-05*mass + 1.5844316265E-01;
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return 5.7465161956E-04*mass - 9.6609929343E-02;
				case 1:
					return 1.5701996550E-07*mass*mass - 4.5980293041E-04*mass + 1.0357381046E+00;
				case 2:
					return -6.0615266088E-08*mass*mass + 3.7528938302E-04*mass + 4.5161981210E-02;
				case 3:
					return 5.6452411568E-08*mass*mass - 9.4249026482E-05*mass + 4.2689388650E-01;
				case 4:
					return -2.5912312180E-08*mass*mass + 1.8910453630E-04*mass + 1.2041317137E-01;
				case 5:
					return 1.6966763663E-08*mass*mass + 2.3841812821E-05*mass + 2.3499342694E-01;
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 3000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return 5.1927873027E-04*mass + 5.5734994317E-02;
				case 1:
					return 2.6356993699E-04*mass + 7.5080657516E-02;
				case 2:
					return 1.7942455026E-04*mass + 7.3207434362E-02;
				case 3:
					return 1.3705983378E-04*mass + 7.0561025542E-02;
				case 4:
					return 1.1117583798E-04*mass + 6.8834545040E-02;
				case 5:
					return 9.3431323732E-05*mass + 6.8159460669E-02;
				case 6:
					return 8.0485041704E-05*mass + 6.8080055740E-02;
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return 5.3543753015E-04*mass + 5.3017546175E-03;
				case 1:
					return 2.5658013201E-04*mass + 1.2063626997E-01;
				case 2:
					return 1.7200148035E-04*mass + 1.2358360913E-01;
				case 3:
					return 1.3106543054E-04*mass + 1.1580990364E-01;
				case 4:
					return 1.0700651965E-04*mass + 1.0670281087E-01;
				case 5:
					return 9.1108828083E-05*mass + 9.8450289793E-02;
				case 6:
					return 7.9442307831E-05*mass + 9.2334890423E-02;
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return 5.5314706631E-04*mass - 5.1325805767E-02;
				case 1:
					return 2.4745051365E-04*mass + 1.7133438587E-01;
				case 2:
					return 1.6691783738E-04*mass + 1.6486760654E-01;
				case 3:
					return 1.2678673747E-04*mass + 1.5359085632E-01;
				case 4:
					return 1.0372003840E-04*mass + 1.3977068399E-01;
				case 5:
					return 8.8460769877E-05*mass + 1.2778494531E-01;
				case 6:
					return 7.7547306931E-05*mass + 1.1772226267E-01;
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return 5.4942059826E-04*mass - 5.5511045650E-02;
				case 1:
					return 2.3092153495E-04*mass + 2.4124276262E-01;
				case 2:
					return 1.6097201620E-04*mass + 2.0608042950E-01;
				case 3:
					return 1.2145356676E-04*mass + 1.9224377311E-01;
				case 4:
					return 1.0006362856E-04*mass + 1.7223755859E-01;
				case 5:
					return 8.5628115812E-05*mass + 1.5633508047E-01;
				case 6:
					return 7.5596403240E-05*mass + 1.4225025190E-01;
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return 5.4768212481E-04*mass - 6.8121942031E-02;
				case 1:
					return 2.1412308071E-04*mass + 3.1117909823E-01;
				case 2:
					return 1.5773096250E-04*mass + 2.3720268376E-01;
				case 3:
					return 1.1682107771E-04*mass + 2.2684396573E-01;
				case 4:
					return 9.7405411941E-05*mass + 2.0001938015E-01;
				case 5:
					return 8.3144699070E-05*mass + 1.8243228233E-01;
				case 6:
					return 7.3663788779E-05*mass + 1.6560334805E-01;
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 4000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return 5.0059769397E-04*mass + 1.3399304111E-01;
				case 1:
					return 2.5503492229E-04*mass + 1.0855085198E-01;
				case 2:
					return 1.7333931713E-04*mass + 9.5133370039E-02;
				case 3:
					return 1.3256893849E-04*mass + 8.5933975272E-02;
				case 4:
					return 1.0793133153E-04*mass + 7.9705619634E-02;
				case 5:
					return 9.1230825683E-05*mass + 7.5665160111E-02;
				case 6:
					return 7.8851927393E-05*mass + 7.3896904125E-02;
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return 5.2112366937E-04*mass + 5.2774347216E-02;
				case 1:
					return 2.5573380138E-04*mass + 1.2501825517E-01;
				case 2:
					return 1.7057071407E-04*mass + 1.2867025075E-01;
				case 3:
					return 1.2916618602E-04*mass + 1.2190697893E-01;
				case 4:
					return 1.0477056586E-04*mass + 1.1365407950E-01;
				case 5:
					return 8.8770575458E-05*mass + 1.0561174590E-01;
				case 6:
					return 7.7399944194E-05*mass + 9.8591965320E-02;
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return 5.2200894269E-04*mass + 4.9446230570E-02;
				case 1:
					return 2.4862834623E-04*mass + 1.7259745818E-01;
				case 2:
					return 1.6433059226E-04*mass + 1.7526276709E-01;
				case 3:
					return 1.2383047405E-04*mass + 1.6442897083E-01;
				case 4:
					return 1.0025374317E-04*mass + 1.5175318889E-01;
				case 5:
					return 8.5000448245E-05*mass + 1.3938982114E-01;
				case 6:
					return 7.4268778278E-05*mass + 1.2847324603E-01;
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return 5.6305474587E-04*mass - 1.0259293056E-01;
				case 1:
					return 2.5893088809E-04*mass + 1.5595702169E-01;
				case 2:
					return 1.7029024422E-04*mass + 1.7612966521E-01;
				case 3:
					return 1.2754224242E-04*mass + 1.7277530152E-01;
				case 4:
					return 1.0279041824E-04*mass + 1.6314225615E-01;
				case 5:
					return 8.6659248725E-05*mass + 1.5263790946E-01;
				case 6:
					return 7.5444787879E-05*mass + 1.4222331242E-01;
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return 5.5656848485E-04*mass - 9.0193796970E-02;
				case 1:
					return 2.4833644284E-04*mass + 2.1220453005E-01;
				case 2:
					return 1.6419441704E-04*mass + 2.1811290035E-01;
				case 3:
					return 1.2309139514E-04*mass + 2.0858307701E-01;
				case 4:
					return 9.9543716772E-05*mass + 1.9365901130E-01;
				case 5:
					return 8.4181039904E-05*mass + 1.7933083034E-01;
				case 6:
					return 7.3320990099E-05*mass + 1.6654889465E-01;
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		return 0.0;
	}
	
	public static double minRS(int idx, double mass, int sulfur) {
		if (mass < Constants.T_MASS1) {
			return Constants.NEGINF;
		}
		else if (mass < 1000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return 4.3061180012E-04*mass - 8.0095224722E-03;
				case 1:
					return 1.6617710152E-07*mass*mass - 5.9150959659E-06*mass + 1.4462164639E-01;
				case 2:
					return 1.7947419053E-04*mass + 4.3023322794E-02;
				case 3:
					return 9.4060035834E-08*mass*mass + 6.2053553190E-06*mass + 8.5937667032E-02;
				case 4:
					return 1.1897106955E-04*mass + 3.5302650205E-02;
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return 4.2429234097E-04*mass - 2.9762724730E-03;
				case 1:
					return 7.1513188877E-12*mass*mass*mass*mass - 2.2056290354E-08*mass*mass*mass + 2.5402528595E-05*mass*mass - 1.2844025000E-02*mass + 2.7416360524E+00;
				case 2:
					return -1.2509853051E-07*mass*mass + 4.0318928429E-04*mass - 2.3237934550E-03;
				case 3:
					return 1.4240750502E-12*mass*mass*mass*mass - 4.3862443546E-09*mass*mass*mass + 5.0416224808E-06*mass*mass - 2.4241036270E-03*mass + 5.7695910166E-01;
				case 4:
					return 5.0298199140E-09*mass*mass + 1.2439421259E-04*mass + 6.4438206800E-02;
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return 4.1837313334E-04*mass - 9.3303789302E-03;
				case 1:
					return 6.7434108365E-12*mass*mass*mass*mass - 2.2711417009E-08*mass*mass*mass + 2.8754396664E-05*mass*mass - 1.6189219184E-02*mass + 3.8801657503E+00;
				case 2:
					return -1.5789024465E-07*mass*mass + 4.9002358299E-04*mass - 3.3812392150E-02;
				case 3:
					return 3.8987917517E-12*mass*mass*mass*mass - 1.2201382958E-08*mass*mass*mass + 1.4274657314E-05*mass*mass - 7.3270528342E-03*mass + 1.6337769947E+00;
				case 4:
					return -1.3598864184E-07*mass*mass + 3.6164397284E-04*mass + 2.1862247749E-05;
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return 4.1490521508E-04*mass - 5.0958937905E-03;
				case 1:
					return 9.6567999679E-07*mass*mass - 1.7870631402E-03*mass + 1.3918716291E+00;
				case 2:
					return 2.7698840194E-04*mass + 4.4342829961E-02;
				case 3:
					return 3.2231749545E-07*mass*mass - 5.2266899885E-04*mass + 5.2571452422E-01;
				case 4:
					return 1.8257257169E-04*mass + 7.3040110315E-02;
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return 4.7569281011E-04*mass - 5.0259818609E-02;
				case 1:
					return 1.7148376554E-06*mass*mass - 3.2490402194E-03*mass + 2.2120814948E+00;
				case 2:
					return -3.5357809611E-07*mass*mass + 8.4012182393E-04*mass - 1.6049514656E-01;
				case 3:
					return -6.8871731127E-10*mass*mass*mass + 2.2264913064E-06*mass*mass - 2.3333002599E-03*mass + 1.1698641549E+00;
				case 4:
					return -3.1451474049E-07*mass*mass + 7.0209771857E-04*mass - 1.2072132850E-01;
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 2000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return 4.5017388899E-04*mass - 2.1782093483E-02;
				case 1:
					return 2.3619792919E-04*mass + 6.7550016211E-02;
				case 2:
					return 1.7340598920E-04*mass + 5.4221966202E-02;
				case 3:
					return 1.3912415122E-04*mass + 4.7461053177E-02;
				case 4:
					return 1.1788453045E-04*mass + 4.1411684320E-02;
				case 5:
					return 1.0313143114E-04*mass + 3.6689933285E-02;
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return 4.4181872607E-04*mass - 9.4862391089E-03;
				case 1:
					return 1.9529723792E-04*mass + 1.9263815311E-01;
				case 2:
					return 1.5787790099E-04*mass + 1.2302790851E-01;
				case 3:
					return 1.2749614761E-04*mass + 1.0518887858E-01;
				case 4:
					return 1.1145835164E-04*mass + 8.4540482547E-02;
				case 5:
					return 9.9095182118E-05*mass + 7.1251236823E-02;
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return 4.4987005701E-04*mass - 2.8683525509E-02;
				case 1:
					return 3.8603195013E-08*mass*mass + 4.2166017720E-05*mass + 3.9438937342E-01;
				case 2:
					return 1.6448535494E-04*mass + 1.4483070760E-01;
				case 3:
					return 1.1440471107E-04*mass + 1.6170198339E-01;
				case 4:
					return 1.0505241824E-04*mass + 1.2508014264E-01;
				case 5:
					return 9.2745316532E-05*mass + 1.0911832520E-01;
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return 4.3635550435E-04*mass - 1.8982616526E-02;
				case 1:
					return 6.7814680688E-08*mass*mass - 8.7222426102E-05*mass + 5.8247099934E-01;
				case 2:
					return 1.7279145395E-04*mass + 1.5468655908E-01;
				case 3:
					return 1.4209932398E-08*mass*mass + 5.4534125999E-05*mass + 2.5379582049E-01;
				case 4:
					return 1.0405268947E-04*mass + 1.5236181580E-01;
				case 5:
					return 8.3493306331E-05*mass + 1.4985997050E-01;
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return 4.1342646205E-04*mass + 5.8666369307E-03;
				case 1:
					return 2.4710312166E-13*mass*mass*mass*mass - 1.6641397976E-09*mass*mass*mass + 4.2251760483E-06*mass*mass - 4.6820144870E-03*mass + 2.5455261574E+00;
				case 2:
					return -3.8860617875E-08*mass*mass + 2.9514650769E-04*mass + 7.8448589669E-02;
				case 3:
					return 3.4191986245E-08*mass*mass - 2.5120537794E-05*mass + 3.6166276039E-01;
				case 4:
					return -2.2845464919E-08*mass*mass + 1.7509468639E-04*mass + 1.2043487485E-01;
				case 5:
					return 1.0274099739E-08*mass*mass + 4.1744064420E-05*mass + 2.1453929960E-01;
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 3000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return 4.4884928473E-04*mass - 2.2605361821E-02;
				case 1:
					return 2.2513850105E-04*mass + 8.7453327375E-02;
				case 2:
					return 1.5847866007E-04*mass + 8.2274939835E-02;
				case 3:
					return 1.2520697030E-04*mass + 7.3822094257E-02;
				case 4:
					return 1.0518066247E-04*mass + 6.5220543834E-02;
				case 5:
					return 9.1288032403E-05*mass + 5.8822238992E-02;
				case 6:
					return 8.1205981398E-05*mass + 5.3444806505E-02;
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return 4.5501123492E-04*mass - 3.2858897309E-02;
				case 1:
					return 2.1051442364E-04*mass + 1.6362283089E-01;
				case 2:
					return 1.4996416682E-04*mass + 1.3855740296E-01;
				case 3:
					return 1.1787902610E-04*mass + 1.2370248474E-01;
				case 4:
					return 9.9665167717E-05*mass + 1.0716696071E-01;
				case 5:
					return 8.7538001200E-05*mass + 9.3465927000E-02;
				case 6:
					return 7.8606487249E-05*mass + 8.2639791878E-02;
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return 4.5965309631E-04*mass - 4.5654590774E-02;
				case 1:
					return 1.9708726253E-04*mass + 2.3500303368E-01;
				case 2:
					return 1.4656205761E-04*mass + 1.7731756599E-01;
				case 3:
					return 1.1351812361E-04*mass + 1.6312000097E-01;
				case 4:
					return 9.6111788179E-05*mass + 1.4172059955E-01;
				case 5:
					return 8.3941594359E-05*mass + 1.2593826410E-01;
				case 6:
					return 7.5404033003E-05*mass + 1.1203126749E-01;
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return 4.5906205281E-04*mass - 5.5603322013E-02;
				case 1:
					return 1.7972459826E-04*mass + 3.1445057435E-01;
				case 2:
					return 1.4350564356E-04*mass + 2.1040218109E-01;
				case 3:
					return 1.0722200780E-04*mass + 2.0462114050E-01;
				case 4:
					return 9.2015864386E-05*mass + 1.7494979903E-01;
				case 5:
					return 7.9965317732E-05*mass + 1.5735596567E-01;
				case 6:
					return 7.1935927393E-05*mass + 1.4063257152E-01;
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return 4.4670362916E-04*mass - 3.5638652907E-02;
				case 1:
					return 1.5975297630E-04*mass + 4.0055552926E-01;
				case 2:
					return 1.4071686049E-04*mass + 2.3968164878E-01;
				case 3:
					return 9.9764673867E-05*mass + 2.4760722533E-01;
				case 4:
					return 8.7929921392E-05*mass + 2.0624229652E-01;
				case 5:
					return 7.5576970297E-05*mass + 1.8837652426E-01;
				case 6:
					return 6.8447720372E-05*mass + 1.6787257907E-01;
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 4000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return 4.8225837504E-04*mass - 1.3994605263E-01;
				case 1:
					return 2.3320278428E-04*mass + 5.8832575026E-02;
				case 2:
					return 1.5726981278E-04*mass + 8.3958625266E-02;
				case 3:
					return 1.2067329133E-04*mass + 8.6036100348E-02;
				case 4:
					return 9.9040943894E-05*mass + 8.2774986370E-02;
				case 5:
					return 8.5074882688E-05*mass + 7.7016370591E-02;
				case 6:
					return 7.4924243024E-05*mass + 7.2105759415E-02;
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return 4.4783733933E-04*mass - 7.1631676687E-03;
				case 1:
					return 2.1702178038E-04*mass + 1.4673559868E-01;
				case 2:
					return 1.4902125593E-04*mass + 1.4242949426E-01;
				case 3:
					return 1.1506840144E-04*mass + 1.3278562496E-01;
				case 4:
					return 9.5341057306E-05*mass + 1.2025553943E-01;
				case 5:
					return 8.2478056406E-05*mass + 1.0826005258E-01;
				case 6:
					return 7.3255730573E-05*mass + 9.8013092994E-02;
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return 4.5365866127E-04*mass - 1.7427964431E-02;
				case 1:
					return 2.0834594659E-04*mass + 2.0420396692E-01;
				case 2:
					return 1.4370230183E-04*mass + 1.8639354359E-01;
				case 3:
					return 1.1062621482E-04*mass + 1.7193111812E-01;
				case 4:
					return 9.1604548455E-05*mass + 1.5509706041E-01;
				case 5:
					return 7.9100546055E-05*mass + 1.4024896881E-01;
				case 6:
					return 7.0403781578E-05*mass + 1.2672824448E-01;
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return 4.5457539034E-04*mass - 3.6523366187E-02;
				case 1:
					return 2.0218234083E-04*mass + 2.4986229708E-01;
				case 2:
					return 1.4118230963E-04*mass + 2.1771374629E-01;
				case 3:
					return 1.0768142574E-04*mass + 2.0398838990E-01;
				case 4:
					return 8.9153267327E-05*mass + 1.8387220436E-01;
				case 5:
					return 7.6886370837E-05*mass + 1.6687771207E-01;
				case 6:
					return 6.8230100810E-05*mass + 1.5199593716E-01;
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return 4.7480187399E-04*mass - 1.2118920896E-01;
				case 1:
					return 2.0060948995E-04*mass + 2.7934325518E-01;
				case 2:
					return 1.4341463786E-04*mass + 2.3067348748E-01;
				case 3:
					return 1.0765826223E-04*mass + 2.2450934221E-01;
				case 4:
					return 8.8952922892E-05*mass + 2.0340154988E-01;
				case 5:
					return 7.5959720972E-05*mass + 1.8776380660E-01;
				case 6:
					return 6.7136463846E-05*mass + 1.7235836654E-01;
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		return 0.0;
	}

	
	public static double avgPkd(int idx, double mass, int sulfur) {
		if (mass < Constants.T_MASS1) {
			return 0;
		}
		else if (mass < 1000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return (1.00293368805513*mass + 0.935708718497934) / (mass + 0.991843047748524);
				case 1:
					return (2.00562314455784*mass - 221.656776224854) / (mass - 110.407196702121);
				case 2:
					return (3.00823137914166*mass + 1.01682121818915) / (mass + 0.464049378904185);
				case 3:
					return (4.01080139863879*mass + 1.14405930657068) / (mass + 0.417163979934148);
				case 4:
					return (5.01324722247508*mass + 1.13386191758042) / (mass + 0.34406742758602);
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return (1.00291099785934*mass + 0.984785545300121) / (mass + 1.07281404422701);
				case 1:
					return (2.00803289720052*mass + 731.334823380727) / (mass + 367.515169996758);
				case 2:
					return (3.01709199956793*mass + 6211.06236068832) / (mass + 2071.55796224543);
				case 3:
					return (4.0140851572608*mass + 3684.46726347861) / (mass + 921.95081597202);
				case 4:
					return (5.01985425060732*mass + 8577.45783212504) / (mass + 1715.60988661398);
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return (1.00291637541836*mass + 0.928853804338038) / (mass + 1.08147448332841);
				case 1:
					return (2.01088446404143*mass + 2053.32124280217) / (mass + 1030.4258103677);
				case 2:
					return (3.01551397580277*mass + 7103.0620404595) / (mass + 2369.91106603995);
				case 3:
					return (4.01272484274554*mass + 2586.19822109) / (mass + 648.605569010246);
				case 4:
					return (5.0217971975202*mass + 9167.8697617209) / (mass + 1835.95930332742);
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return (1.00288625044951*mass + 0.926908248509697) / (mass + 1.11296410908056);
				case 1:
					return (2.01337509765593*mass + 3565.73960699026) / (mass + 1788.78805910417);
				case 2:
					return (3.01122865955959*mass + 6271.43771568732) / (mass + 2092.61554881518);
				case 3:
					return (4.01519969384884*mass + 4406.56372430441) / (mass + 1105.05103602182);
				case 4:
					return (5.01991498112207*mass + 10921.2156880335) / (mass + 2187.48723052913);
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return (1.0030321936735*mass + 0.737050658349466) / (mass + 1.12883574674179);
				case 1:
					return (2.00407817184208*mass + 410.115264013997) / (mass + 207.29217595888);
				case 2:
					return (3.00332740006045*mass + 1.07671551308492) / (mass + 1.18247440425356);
				case 3:
					return (4.00478649697988*mass + 1.12899575033535) / (mass + 1.69452199570285);
				case 4:
					return (5.00965848278663*mass + 4048.41296495178) / (mass + 811.731474375441);
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 2000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return (1.00295950840363*mass + 0.961220787174525) / (mass + 1.0432566509175);
				case 1:
					return (2.00582467881989*mass + 1.15374246709954) / (mass + 0.800596318050647);
				case 2:
					return (3.0089236837415*mass + 3269.10662204862) / (mass + 1087.2042321278);
				case 3:
					return (4.01179609257744*mass + 4924.94934254691) / (mass + 1228.45714761347);
				case 4:
					return (5.01464700270306*mass + 7062.1229917312) / (mass + 1409.24510511963);
				case 5:
					return (6.01744041514622*mass + 9112.87507338323) / (mass + 1515.41110829221);
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return (1.00295391404271*mass + 0.936833366445414) / (mass + 1.06602995116468);
				case 1:
					return (2.00635440122245*mass - 189.492094071571) / (mass - 92.9923993436242);
				case 2:
					return (3.00979900966767*mass + 1163.23995850531) / (mass + 388.997415270093);
				case 3:
					return (4.01299212577052*mass + 2560.2981780313) / (mass + 641.034069885366);
				case 4:
					return (5.01624231676911*mass + 4848.70792374594) / (mass + 970.202393390242);
				case 5:
					return (6.01919476509136*mass + 7144.1941173926) / (mass + 1190.73060128619);
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return (1.00296256122887*mass + 0.896819931357756) / (mass + 1.09795950658375);
				case 1:
					return (2.00694894057729*mass + 237.911429443904) / (mass + 121.519283577764);
				case 2:
					return (3.01112830009987*mass + 3024.43230146041) / (mass + 1010.09485470343);
				case 3:
					return (4.01326486350479*mass + 2764.78303634486) / (mass + 693.360262343827);
				case 4:
					return (5.01682607964317*mass + 5167.43675442151) / (mass + 1035.44258460745);
				case 5:
					return (6.01942133732949*mass + 6506.50241997006) / (mass + 1086.19597107412);
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return (1.00294844020763*mass + 0.886586382192758) / (mass + 1.13914252828926);
				case 1:
					return (2.00792290741722*mass + 963.327994739202) / (mass + 485.232160483669);
				case 2:
					return (3.01407118304633*mass + 6834.10757976449) / (mass + 2280.97656077705);
				case 3:
					return (4.01398569007967*mass + 3728.26147537136) / (mass + 935.20928502631);
				case 4:
					return (5.01976415736369*mass + 9092.58590659403) / (mass + 1821.78058232433);
				case 5:
					return (6.02058575768338*mass + 8275.99105265916) / (mass + 1382.31038260688);
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return (1.00294975458523*mass + 0.795270993463695) / (mass + 1.12058198075075);
				case 1:
					return (2.00829863859034*mass + 1348.69158661254) / (mass + 678.746948940858);
				case 2:
					return (3.01380023238372*mass + 7563.69877808966) / (mass + 2524.81290637155);
				case 3:
					return (4.01307744357889*mass + 3389.60777686819) / (mass + 851.043027508048);
				case 4:
					return (5.01934776906142*mass + 9917.04844256749) / (mass + 1987.48494150336);
				case 5:
					return (6.01833105585862*mass + 6826.91276199546) / (mass + 1141.29058866015);
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 3000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return (1.00297140085314*mass + 0.95034341462094) / (mass + 1.04977017899761);
				case 1:
					return (2.00590021033957*mass + 1.08483648507666) / (mass + 0.831484129872426);
				case 2:
					return (3.00873622055941*mass + 1.09850592888673) / (mass + 0.713464047008087);
				case 3:
					return (4.01151276645993*mass + 1.09348714183758) / (mass + 0.650036099876211);
				case 4:
					return (5.01423865186348*mass + 1.08784891914275) / (mass + 0.607863374850521);
				case 5:
					return (6.01692357276294*mass + 1.08317681684485) / (mass + 0.5765590498714);
				case 6:
					return (7.01957541173439*mass + 1.080249139226) / (mass + 0.551977239096394);
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return (1.00297459032596*mass + 0.962637397483761) / (mass + 1.1255468283858);
				case 1:
					return (2.0062400511767*mass - 259.653740470464) / (mass - 128.110083487467);
				case 2:
					return (3.00929202431592*mass + 0.697974607174444) / (mass + 1.99847979538262);
				case 3:
					return (4.01207767986914*mass + 0.82410656367374) / (mass + 2.04962080421739);
				case 4:
					return (5.01595158293566*mass + 4220.81708231031) / (mass + 844.754187791304);
				case 5:
					return (6.01914034963153*mass + 6910.77707384658) / (mass + 1151.87752964792);
				case 6:
					return (7.02224690721867*mass + 9998.00824812674) / (mass + 1427.86247325871);
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return (1.00297704572907*mass + 0.828570678572173) / (mass + 1.05894515179825);
				case 1:
					return (2.00625061697481*mass - 831.861219352272) / (mass - 412.953389074961);
				case 2:
					return (3.0095899306176*mass + 391.79843701314) / (mass + 133.125294325674);
				case 3:
					return (4.01274416965519*mass + 1416.61218580688) / (mass + 356.6181837333);
				case 4:
					return (5.01589393114694*mass + 3013.17909069383) / (mass + 604.899385173173);
				case 5:
					return (6.01889806471087*mass + 4714.5143571775) / (mass + 787.813467814725);
				case 6:
					return (7.02186375661506*mass + 6802.56015896404) / (mass + 973.602331073011);
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return (1.00298434194247*mass + 0.927977091477208) / (mass + 1.24546855596251);
				case 1:
					return (2.00667413799264*mass - 224.793445993497) / (mass - 109.023943221837);
				case 2:
					return (3.01069460062557*mass + 2581.87768722357) / (mass + 863.457611915027);
				case 3:
					return (4.0139480021359*mass + 3700.31051629477) / (mass + 928.21165187131);
				case 4:
					return (5.01760377038444*mass + 6404.67985819644) / (mass + 1283.99307375765);
				case 5:
					return (6.02086362993975*mass + 8677.03311020124) / (mass + 1449.16458292866);
				case 6:
					return (7.02415472685762*mass + 11523.376587455) / (mass + 1649.0385153186);
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return (1.0029718956106*mass + 0.964028589393642) / (mass + 1.32986840065261);
				case 1:
					return (2.00708666662462*mass + 325.84806064546) / (mass + 166.852815112311);
				case 2:
					return (3.01193446136126*mass + 5014.50896164342) / (mass + 1674.87278215358);
				case 3:
					return (4.01491802297074*mass + 5543.06841822767) / (mass + 1389.75289767779);
				case 4:
					return (5.01922815690291*mass + 9685.33766490196) / (mass + 1941.11356597271);
				case 5:
					return (6.02237363200899*mass + 11846.4998010726) / (mass + 1978.49928794925);
				case 6:
					return (7.02595872202835*mass + 15414.1823392104) / (mass + 2206.05284697366);
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 4000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return (1.00296823886383*mass + 0.958626307824402) / (mass + 1.03943968145155);
				case 1:
					return (2.0059148433284*mass + 1.06451265197153) / (mass + 0.833885772112208);
				case 2:
					return (3.00879893276037*mass + 1.07963130337461) / (mass + 0.759735458311488);
				case 3:
					return (4.01163238245427*mass + 1.06877614548068) / (mass + 0.723272929555923);
				case 4:
					return (5.01442121899282*mass + 1.05944403126076) / (mass + 0.701424424561989);
				case 5:
					return (6.01717260269282*mass + 1.05236448749887) / (mass + 0.685588002501871);
				case 6:
					return (7.01989196458963*mass + 1.04699725819895) / (mass + 0.672669281902025);
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return (1.00297435514572*mass + 0.934262842144595) / (mass + 1.08949432946512);
				case 1:
					return (2.00615439272962*mass + 0.865452890131085) / (mass + 1.6701548596557);
				case 2:
					return (3.00922404951199*mass + 1.22185688318408) / (mass + 2.09587032177023);
				case 3:
					return (4.01304874878666*mass + 5525.43184421457) / (mass + 1380.58891932832);
				case 4:
					return (5.0149961815169*mass + 0.791674871150801) / (mass + 2.1429669517967);
				case 5:
					return (6.01771216641173*mass + 0.839058104351687) / (mass + 2.13532378661159);
				case 6:
					return (7.02033652357364*mass + 0.879827172324812) / (mass + 2.09089302941654);
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return (1.00298997140931*mass + 0.869642399263411) / (mass + 1.13597761193287);
				case 1:
					return (2.00615904192065*mass - 1312.54250250541) / (mass - 652.838154693951);
				case 2:
					return (3.00953001940095*mass + 0.172811477399175) / (mass + 2.81384110096259);
				case 3:
					return (4.01249021516468*mass + 1.07593142277282) / (mass + 3.29006712559875);
				case 4:
					return (5.0159938929733*mass + 2741.7837952117) / (mass + 550.785907604928);
				case 5:
					return (6.01793724914704*mass + 0.661719558950418) / (mass + 3.22554468450629);
				case 6:
					return (7.02047793439992*mass + 0.735200769449895) / (mass + 3.17160826502371);
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return (1.00299272284475*mass + 0.812562066566496) / (mass + 1.15662451932822);
				case 1:
					return (2.00620216775292*mass - 1480.03158227583) / (mass - 735.903038839805);
				case 2:
					return (3.00940966219523*mass - 791.098670837569) / (mass - 259.855795010345);
				case 3:
					return (4.012541042136*mass + 0.545661107405565) / (mass + 3.94992730339765);
				case 4:
					return (5.0158262966805*mass + 1808.16271067829) / (mass + 365.243030043943);
				case 5:
					return (6.01779277548327*mass + 0.542167332170347) / (mass + 3.97951179751125);
				case 6:
					return (7.02022449534525*mass + 0.638546461458168) / (mass + 3.92320915955648);
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return (1.00299940857136*mass + 0.769304143039377) / (mass + 1.21003144750379);
				case 1:
					return (2.00641938919981*mass - 1011.13550924436) / (mass - 501.225352629905);
				case 2:
					return (3.01015237774343*mass + 1318.49865388763) / (mass + 443.343353022124);
				case 3:
					return (4.01325448124215*mass + 1986.59814881933) / (mass + 500.847783353272);
				case 4:
					return (5.01662123281798*mass + 4056.51878915689) / (mass + 815.494008217453);
				case 5:
					return (6.0197921777858*mass + 6018.01280266545) / (mass + 1007.18169097412);
				case 6:
					return (7.01984786316731*mass + 0.617668673230637) / (mass + 4.51373751614008);
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		return 0.0;
	}

	public static double maxPkd(int idx, double mass, int sulfur) {
		if (mass < Constants.T_MASS1) {
			return Constants.POSINF;
		}
		else if (mass < 1000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return (1.00323465009313*mass - 101.911088369617) / (mass - 101.486938597329);
				case 1:
					return (2.00603412691658*mass - 549.380179306983) / (mass - 273.796404071801);
				case 2:
					return (3.00887783994786*mass - 589.081953970283) / (mass - 195.684920408319);
				case 3:
					return (4.01158422185472*mass - 808.313473535027) / (mass - 201.400808362753);
				case 4:
					return (5.01428420172806*mass - 863.546362846799) / (mass - 172.115037473941);
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return (1.00318338360925*mass - 120.207453299689) / (mass - 119.713629153107);
				case 1:
					return (2.00702626236606*mass + 46.5872125267643) / (mass + 24.8714019470776);
				case 2:
					return (3.01261334148514*mass + 2731.58093331645) / (mass + 911.564702777241);
				case 3:
					return (4.01124787256745*mass + 1011.39355431633) / (mass + 253.695834016275);
				case 4:
					return (5.01502143472786*mass + 3086.86184215871) / (mass + 617.838569401118);
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return (1.00311010427832*mass - 127.192266361806) / (mass - 126.673577057391);
				case 1:
					return (2.00954485538437*mass + 1087.93081355998) / (mass + 546.97165845165);
				case 2:
					return (3.01453984239629*mass + 5429.39153358859) / (mass + 1811.72199676215);
				case 3:
					return (4.01138762371718*mass + 1359.82242819306) / (mass + 341.65183142468);
				case 4:
					return (5.01944500729388*mass + 6629.20117400348) / (mass + 1327.80028566915);
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return (1.00310540349736*mass - 173.537189442841) / (mass - 172.821047756357);
				case 1:
					return (2.01692211105366*mass + 3692.45867567624) / (mass + 1853.08563314816);
				case 2:
					return (3.0131505109971*mass + 6034.09324200222) / (mass + 2013.75873399406);
				case 3:
					return (4.01648474689816*mass + 3819.7737289407) / (mass + 958.356741891123);
				case 4:
					return (5.02183632313415*mass + 10088.0776162089) / (mass + 2020.92657555153);
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return (1.00316322579353*mass - 281.984337002833) / (mass - 280.814205585822);
				case 1:
					return (2.00595228581962*mass + 677.397108974326) / (mass + 341.695858463498);
				case 2:
					return (3.00402027608075*mass + 0.969423881734885) / (mass + 1.29992760793369);
				case 3:
					return (4.00708103746794*mass + 448.866194483685) / (mass + 114.111417914335);
				case 4:
					return (5.00843223890757*mass + 2161.95765465304) / (mass + 434.085459576455);
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 2000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return (1.00302145240148*mass + 1.02791661191195) / (mass + 0.966054637287335);
				case 1:
					return (2.00594887203368*mass + 1.20354318543046) / (mass + 0.680507716038747);
				case 2:
					return (3.00876474452885*mass + 1.25215336872159) / (mass + 0.531918420354397);
				case 3:
					return (4.01151434281288*mass + 1.27976885126513) / (mass + 0.452656842699228);
				case 4:
					return (5.01421178440333*mass + 1.30507459490297) / (mass + 0.402590738466823);
				case 5:
					return (6.01686546043182*mass + 1.32729362331526) / (mass + 0.36655966505376);
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return (1.00304030219015*mass + 0.987266673475433) / (mass + 1.01031566018996);
				case 1:
					return (2.00630304719678*mass - 429.030162344813) / (mass - 212.819877631674);
				case 2:
					return (3.00952599393135*mass + 331.91297025599) / (mass + 112.000696254764);
				case 3:
					return (4.01264268132983*mass + 1431.13974032514) / (mass + 358.833172177393);
				case 4:
					return (5.01574194063236*mass + 3137.96046935861) / (mass + 628.209851400963);
				case 5:
					return (6.01862377453611*mass + 4978.74785167351) / (mass + 830.012477421433);
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return (1.00305807953947*mass + 0.9415027883532) / (mass + 1.05602907023557);
				case 1:
					return (2.0068728242289*mass - 90.7199121718033) / (mass - 42.9363526275577);
				case 2:
					return (3.01096491256242*mass + 2136.91240239059) / (mass + 714.222663335546);
				case 3:
					return (4.01345349432759*mass + 2268.14670789103) / (mass + 569.102060391229);
				case 4:
					return (5.01683598124571*mass + 4161.68383747292) / (mass + 834.218877389377);
				case 5:
					return (6.01956689135986*mass + 5529.75539038818) / (mass + 923.329868336394);
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return (1.00309492341683*mass + 0.907174939359879) / (mass + 1.1498491485713);
				case 1:
					return (2.00863442990557*mass + 1061.69428212719) / (mass + 534.396829522044);
				case 2:
					return (3.01439131399044*mass + 5845.13687124386) / (mass + 1951.29792209599);
				case 3:
					return (4.01660163959572*mass + 5095.66595458407) / (mass + 1277.32266000289);
				case 4:
					return (5.02093292530504*mass + 8571.42426076109) / (mass + 1717.56917134946);
				case 5:
					return (6.02487348808437*mass + 10985.6277241118) / (mass + 1834.26077516118);
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return (1.00305051163165*mass - 468.584359443262) / (mass - 466.983744015617);
				case 1:
					return (2.00768890204288*mass + 410.668742514545) / (mass + 209.258313423963);
				case 2:
					return (3.0141646505511*mass + 6258.78185051894) / (mass + 2089.89247552913);
				case 3:
					return (4.01296752983917*mass + 2223.84430904944) / (mass + 559.406890227186);
				case 4:
					return (5.01868684339485*mass + 7293.8263955581) / (mass + 1462.55553474879);
				case 5:
					return (6.0191897066367*mass + 5938.90467951422) / (mass + 993.286884487851);
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 3000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return (1.00301505916509*mass + 1.04910963517566) / (mass + 0.951003579478781);
				case 1:
					return (2.00605586125415*mass + 1.13254662579167) / (mass + 0.731375529772049);
				case 2:
					return (3.00900586210842*mass + 1.11846352093606) / (mass + 0.626386144902348);
				case 3:
					return (4.01189325956952*mass + 1.09561816818437) / (mass + 0.573787698989272);
				case 4:
					return (5.01473031057534*mass + 1.07561522869956) / (mass + 0.54207216932026);
				case 5:
					return (6.01752444216044*mass + 1.05936764342531) / (mass + 0.519394548585614);
				case 6:
					return (7.02028045076763*mass + 1.04656672135212) / (mass + 0.501119090776945);
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return (1.00304219750861*mass + 0.98787522950027) / (mass + 1.01183285388182);
				case 1:
					return (2.00644430353369*mass + 0.666393628233783) / (mass + 1.58876826964424);
				case 2:
					return (3.00960483724371*mass + 0.662185732564852) / (mass + 1.87382260898857);
				case 3:
					return (4.01254786259615*mass + 0.74031748140113) / (mass + 1.96830227828186);
				case 4:
					return (5.01529433455163*mass + 0.838186881159386) / (mass + 1.95087251147434);
				case 5:
					return (6.01789104555109*mass + 0.926440565057262) / (mass + 1.87849316222549);
				case 6:
					return (7.02038104200768*mass + 1.00019157236017) / (mass + 1.78349149216224);
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return (1.00306696325623*mass + 0.926190753557292) / (mass + 1.07526276473253);
				case 1:
					return (2.00683634706047*mass + 0.28069657269397) / (mass + 2.46731714149232);
				case 2:
					return (3.00993814851697*mass + 0.573967495533016) / (mass + 2.88455321183324);
				case 3:
					return (4.01280033781308*mass + 0.872273460047192) / (mass + 3.02866331618853);
				case 4:
					return (5.01903786133895*mass + 7771.74528176492) / (mass + 1555.90880077518);
				case 5:
					return (6.02297228967683*mass + 11311.2145046748) / (mass + 1886.48074627407);
				case 6:
					return (7.02675962631248*mass + 15220.8288794449) / (mass + 2175.3486066783);
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return (1.00307318948725*mass + 0.888187227580665) / (mass + 1.11524195761261);
				case 1:
					return (2.00702529895715*mass + 0.0183965449294023) / (mass + 3.1063429207308);
				case 2:
					return (3.01208530470937*mass + 4099.06480341464) / (mass + 1368.79529713624);
				case 3:
					return (4.01613023949351*mass + 6302.52174981496) / (mass + 1578.44312910256);
				case 4:
					return (5.02032306138408*mass + 9576.69627507416) / (mass + 1918.17533359631);
				case 5:
					return (6.02405616324884*mass + 12505.4409960939) / (mass + 2087.04915790294);
				case 6:
					return (7.02771588303639*mass + 15864.2543499862) / (mass + 2268.99997026786);
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return (1.00308328214053*mass + 0.839304695949488) / (mass + 1.15917347414996);
				case 1:
					return (2.0071936449972*mass - 0.419878734727797) / (mass + 3.62321348435638);
				case 2:
					return (3.00968086409406*mass + 0.383494024973809) / (mass + 3.80595543540441);
				case 3:
					return (4.01230582202874*mass + 0.715495817882188) / (mass + 4.01174450447243);
				case 4:
					return (5.01741448112024*mass + 4625.05028756207) / (mass + 928.912950872199);
				case 5:
					return (6.02039857041521*mass + 6119.03213000374) / (mass + 1023.73554043879);
				case 6:
					return (7.02348749262967*mass + 8225.94059101681) / (mass + 1178.96099077579);
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 4000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return (1.00301407310622*mass + 1.10961192510235) / (mass + 0.975043010510933);
				case 1:
					return (2.00603837764656*mass + 1.14921820004958) / (mass + 0.697384840109725);
				case 2:
					return (3.00898640970256*mass + 1.13160245108093) / (mass + 0.591872667769032);
				case 3:
					return (4.01188280486344*mass + 1.1066027130129) / (mass + 0.54529803604492);
				case 4:
					return (5.01473614775442*mass + 1.08589954166806) / (mass + 0.52167786281039);
				case 5:
					return (6.01755689273087*mass + 1.06922284653966) / (mass + 0.510029643252435);
				case 6:
					return (7.02034427528611*mass + 1.05617200148107) / (mass + 0.501602757510202);
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return (1.00301879699861*mass + 1.03110355480469) / (mass + 0.968780767014032);
				case 1:
					return (2.00618451683076*mass + 0.861037167932796) / (mass + 1.27783400804699);
				case 2:
					return (3.00930655099496*mass + 0.800765793742524) / (mass + 1.59802894995137);
				case 3:
					return (4.01233347050461*mass + 0.804046678977369) / (mass + 1.79590362498352);
				case 4:
					return (5.01524947384032*mass + 0.829150383026894) / (mass + 1.89204112768305);
				case 5:
					return (6.0180582092141*mass + 0.860332306297885) / (mass + 1.919411239374);
				case 6:
					return (7.02077105686772*mass + 0.890873299958481) / (mass + 1.90281654572774);
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return (1.00300281458096*mass + 1.03458905640391) / (mass + 0.965399031375029);
				case 1:
					return (2.00636494397111*mass - 0.197203728103424) / (mass + 1.4927522086572);
				case 2:
					return (3.00957556688493*mass + 0.503339117983076) / (mass + 2.45088972333953);
				case 3:
					return (4.01262643114986*mass + 0.547460979628313) / (mass + 2.76429638093308);
				case 4:
					return (5.01551414179175*mass + 0.613424298802826) / (mass + 2.90647453887888);
				case 5:
					return (6.01826234105047*mass + 0.678138120737821) / (mass + 2.94954040708378);
				case 6:
					return (7.02089128591202*mass + 0.736381772696561) / (mass + 2.93239725888199);
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return (1.00308571561896*mass + 0.88778075046291) / (mass + 1.15961862587008);
				case 1:
					return (2.00681088209776*mass + 0.268335929186727) / (mass + 2.94160691746875);
				case 2:
					return (3.01016852195219*mass + 2.33725063476955) / (mass + 4.35845876715403);
				case 3:
					return (4.01329981569193*mass + 0.304160171715496) / (mass + 4.0037123213406);
				case 4:
					return (5.01622085552407*mass + 0.460961045843221) / (mass + 4.15249541265207);
				case 5:
					return (6.01897204927991*mass + 0.596453694547288) / (mass + 4.1828878348407);
				case 6:
					return (7.02158100619135*mass + 0.710600363175461) / (mass + 4.14481410955902);
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return (1.00307739766197*mass + 0.828578743624013) / (mass + 1.13721980515349);
				case 1:
					return (2.0064155743967*mass - 1821.11001724222) / (mass - 905.764734270153);
				case 2:
					return (3.01003301825426*mass - 508.747552261909) / (mass - 165.230677348243);
				case 3:
					return (4.01325107553331*mass + 0.0743829075235624) / (mass + 4.54119888308982);
				case 4:
					return (5.01604508126254*mass + 0.245330966915111) / (mass + 4.67756639572166);
				case 5:
					return (6.01866589271048*mass + 0.386993071719538) / (mass + 4.70027185463661);
				case 6:
					return (7.02113909992155*mass + 0.50570120043085) / (mass + 4.65607918112508);
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		return 0.0;
	}
	
	public static double minPkd(int idx, double mass, int sulfur) {
		if (mass < Constants.T_MASS1) {
			return Constants.NEGINF;
		}
		else if (mass < 1000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return (1.00247519086927*mass - 363.635280706548) / (mass - 362.742196679545);
				case 1:
					return (2.00505038778732*mass + 1331.72816955264) / (mass + 664.408059642483);
				case 2:
					return (3.0071967921405*mass + 1564.169729515) / (mass + 520.219552188265);
				case 3:
					return (4.00940136098405*mass + 1.18076766815339) / (mass + 0.32505205151551);
				case 4:
					return (5.01161863770218*mass + 1.16212584883309) / (mass + 0.255212910293626);
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return (1.00245585149465*mass + 1.00990425257206) / (mass + 0.987844635633999);
				case 1:
					return (2.01101696516683*mass + 2246.3158985131) / (mass + 1125.86541834185);
				case 2:
					return (3.0191830678755*mass + 9251.41387991941) / (mass + 3085.29216899982);
				case 3:
					return (4.02069677261752*mass + 10154.3807939149) / (mass + 2539.57805265345);
				case 4:
					return (5.02078063865046*mass + 12642.9367997232) / (mass + 2528.4898625405);
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return (1.00241138291841*mass + 1.00772009385201) / (mass + 0.991924106345017);
				case 1:
					return (2.01302911871099*mass + 3547.06741885911) / (mass + 1778.54492372203);
				case 2:
					return (3.01004947174719*mass + 5488.91834918521) / (mass + 1831.2198792425);
				case 3:
					return (4.01298903627605*mass + 3872.15433104658) / (mass + 970.364792734504);
				case 4:
					return (5.01679005639199*mass + 8757.90868886194) / (mass + 1753.62201524984);
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return (1.00232416970916*mass + 0.981874560809136) / (mass + 0.905624160374287);
				case 1:
					return (2.01282575039249*mass + 4492.81246051024) / (mass + 2252.74220802513);
				case 2:
					return (3.00623165868734*mass + 4493.07908964022) / (mass + 1498.9720734543);
				case 3:
					return (4.01474328957685*mass + 5740.35954480459) / (mass + 1438.77737952656);
				case 4:
					return (5.0116124254603*mass + 8131.02435164025) / (mass + 1628.38565260837);
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return (1.00267609186131*mass + 0.887658465616967) / (mass + 1.04915787098515);
				case 1:
					return (2.00241999892306*mass + 0.591657772538667) / (mass + 1.77913404882978);
				case 2:
					return (3.00263863676066*mass + 0.999331613722848) / (mass + 1.00267304259035);
				case 3:
					return (4.0038687090921*mass + 0.863261708332881) / (mass + 1.47547607743452);
				case 4:
					return (5.00314077661724*mass + 1.00854004435314) / (mass + 0.972787519158438);
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 2000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return (1.00293894899701*mass + 0.780335877900164) / (mass + 1.16285996898145);
				case 1:
					return (2.00568990117336*mass + 1.11251159337178) / (mass + 0.971706921645214);
				case 2:
					return (3.00832111160642*mass + 1.24535767902664) / (mass + 0.80159270057117);
				case 3:
					return (4.01089286157339*mass + 1.33903953371736) / (mass + 0.702287302126615);
				case 4:
					return (5.01341638694903*mass + 1.39530993977548) / (mass + 0.627661813570923);
				case 5:
					return (6.01590735019695*mass + 1.44366036779752) / (mass + 0.573679458705331);
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return (1.00287516555065*mass + 0.866687431949152) / (mass + 1.18760369509669);
				case 1:
					return (2.00624269436416*mass + 0.430671688167438) / (mass + 2.01097322205534);
				case 2:
					return (3.00941656331905*mass + 1089.28415664509) / (mass + 364.539831770761);
				case 3:
					return (4.01248753524962*mass + 2403.38894609511) / (mass + 602.034566795384);
				case 4:
					return (5.01551469610518*mass + 4328.90995607133) / (mass + 866.565689411989);
				case 5:
					return (6.01836738415094*mass + 6430.42826154745) / (mass + 1072.13335106877);
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return (1.00285357319135*mass + 0.92988510110626) / (mass + 1.27643906488388);
				case 1:
					return (2.00739028476144*mass + 826.59295936924) / (mass + 416.165177265183);
				case 2:
					return (3.0125684750019*mass + 4977.56424255715) / (mass + 1661.49784902856);
				case 3:
					return (4.01467725622738*mass + 4706.32478508556) / (mass + 1179.11011254176);
				case 4:
					return (5.01896742492162*mass + 8123.99458262906) / (mass + 1627.161686528);
				case 5:
					return (6.02217744618157*mass + 10452.742230323) / (mass + 1744.28241897821);
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return (1.00278483346136*mass + 0.850525557734518) / (mass + 1.16026667965768);
				case 1:
					return (2.00791238517161*mass + 1392.07749884338) / (mass + 699.861157921907);
				case 2:
					return (3.01277448326267*mass + 6331.34864065179) / (mass + 2113.51037412827);
				case 3:
					return (4.01355646237528*mass + 4153.34226026309) / (mass + 1041.63185077373);
				case 4:
					return (5.01798699717837*mass + 8255.89613943105) / (mass + 1654.48538421565);
				case 5:
					return (6.01896918624442*mass + 7664.02932128535) / (mass + 1280.35817561694);
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return (1.00266086308364*mass + 0.97356498733149) / (mass + 1.15099426520165);
				case 1:
					return (2.00865206075014*mass + 2320.81138513174) / (mass + 1165.2906256898);
				case 2:
					return (3.01229480263002*mass + 7762.21037392657) / (mass + 2590.77354062871);
				case 3:
					return (4.01367095749154*mass + 5229.22723849497) / (mass + 1311.35921562771);
				case 4:
					return (5.0193812446659*mass + 12267.1141148875) / (mass + 2457.78357665112);
				case 5:
					return (6.01999239814899*mass + 10554.332941259) / (mass + 1763.09751616537);
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 3000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return (1.00293212654061*mass + 0.796501365099407) / (mass + 1.20420013872846);
				case 1:
					return (2.00577291826055*mass + 0.994445812260939) / (mass + 1.01026063806221);
				case 2:
					return (3.00851155286873*mass + 1.0443819525622) / (mass + 0.871538620830198);
				case 3:
					return (4.01118928148316*mass + 1.05618111855945) / (mass + 0.786100631833679);
				case 4:
					return (5.01382101453159*mass + 1.05731117865914) / (mass + 0.726518940245762);
				case 5:
					return (6.01641442855417*mass + 1.05509417791109) / (mass + 0.681236462151582);
				case 6:
					return (7.0189779795808*mass + 1.05214344692937) / (mass + 0.645381190801582);
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return (1.00292263745227*mass + 0.79903995306618) / (mass + 1.20482600816886);
				case 1:
					return (2.00616856702034*mass + 0.483778922984113) / (mass + 1.93588620814729);
				case 2:
					return (3.00901007618355*mass + 0.614979575011007) / (mass + 2.10764736768077);
				case 3:
					return (4.01168316162039*mass + 0.740320070610099) / (mass + 2.13026201473617);
				case 4:
					return (5.01417787544454*mass + 0.865856177996661) / (mass + 2.05727009620231);
				case 5:
					return (6.0165581880339*mass + 0.967043674437596) / (mass + 1.95169927168619);
				case 6:
					return (7.01885384507897*mass + 1.04814605171985) / (mass + 1.83500956703396);
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return (1.00293597701617*mass + 0.755284841473897) / (mass + 1.25421618889217);
				case 1:
					return (2.00650713403873*mass + 0.127177059919506) / (mass + 2.83248754770818);
				case 2:
					return (3.00916778524052*mass + 0.547594141117172) / (mass + 3.00415510524145);
				case 3:
					return (4.0117619957278*mass + 0.822086298589229) / (mass + 3.08566001687468);
				case 4:
					return (5.01667065382035*mass + 5445.09988430536) / (mass + 1091.25620672233);
				case 5:
					return (6.01973135381839*mass + 7465.64232878279) / (mass + 1246.33392547954);
				case 6:
					return (7.02282975972314*mass + 10082.6861272135) / (mass + 1442.18544758555);
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return (1.00289231518375*mass + 0.769511660153019) / (mass + 1.2309009789379);
				case 1:
					return (2.00658998153334*mass - 0.0991151140417467) / (mass + 3.46331482121619);
				case 2:
					return (3.01080103868897*mass + 3609.53069472763) / (mass + 1206.06798247675);
				case 3:
					return (4.01393840439701*mass + 4639.56575557761) / (mass + 1163.18127564152);
				case 4:
					return (5.01688910157483*mass + 6308.54146316528) / (mass + 1264.94444373522);
				case 5:
					return (6.01955131764776*mass + 7584.11920321323) / (mass + 1267.19977644073);
				case 6:
					return (7.02215344732951*mass + 9190.64838969138) / (mass + 1315.98699022132);
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return (1.00281919594529*mass + 0.817795139604754) / (mass + 1.19163470672453);
				case 1:
					return (2.0069998520841*mass + 727.249142573273) / (mass + 367.652478391595);
				case 2:
					return (3.01242766063648*mass + 7261.45775080779) / (mass + 2423.87580757177);
				case 3:
					return (4.01542041242909*mass + 7830.54862893935) / (mass + 1961.82973643875);
				case 4:
					return (5.02006087699485*mass + 13183.3098842564) / (mass + 2640.95215663264);
				case 5:
					return (6.02310531876644*mass + 15465.7608371081) / (mass + 2581.96064305708);
				case 6:
					return (7.02658783119541*mass + 19331.6513520466) / (mass + 2765.92374219067);
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		else if (mass < 4000.0) {
			switch (sulfur) {
			case 0:
				switch (idx) {
				case 0:
					return (1.00290506761367*mass + 0.861759290088543) / (mass + 1.19342269510149);
				case 1:
					return (2.00576517661949*mass + 0.999903082920408) / (mass + 1.00019435134877);
				case 2:
					return (3.00859642288693*mass + 1.02247055495623) / (mass + 0.951654902734284);
				case 3:
					return (4.01134737785316*mass + 1.03407392100934) / (mass + 0.900687148233851);
				case 4:
					return (5.01405424164148*mass + 1.03919309307976) / (mass + 0.864413871217093);
				case 5:
					return (6.01672340069931*mass + 1.04204143453997) / (mass + 0.834977323856892);
				case 6:
					return (7.01936285511274*mass + 1.04383976216282) / (mass + 0.810685643038795);
				default:
					return 0.0;
				}
			case 1:
				switch (idx) {
				case 0:
					return (1.00293295174458*mass + 0.893349411519525) / (mass + 1.33718472219975);
				case 1:
					return (2.00608690686286*mass + 0.557308974379463) / (mass + 1.85344827832154);
				case 2:
					return (3.00908351891954*mass + 0.607844417411841) / (mass + 2.17297863575655);
				case 3:
					return (4.01194615699251*mass + 0.674299421185616) / (mass + 2.30203123575115);
				case 4:
					return (5.01467796635756*mass + 0.740724520477367) / (mass + 2.32132423619774);
				case 5:
					return (6.01730363920719*mass + 0.798170322325643) / (mass + 2.28432313969019);
				case 6:
					return (7.0198428492951*mass + 0.846120263111408) / (mass + 2.21769382307467);
				default:
					return 0.0;
				}
			case 2:
				switch (idx) {
				case 0:
					return (1.00293300095178*mass + 0.759913004149397) / (mass + 1.2565486926105);
				case 1:
					return (2.00621543554579*mass - 0.654919585581755) / (mass + 1.99722548508934);
				case 2:
					return (3.00918002398859*mass - 0.616699799461956) / (mass + 2.60825240549106);
				case 3:
					return (4.01199311483407*mass + 0.457200087478591) / (mass + 3.1461008830423);
				case 4:
					return (5.01466030319039*mass + 0.556280596641475) / (mass + 3.2115528517585);
				case 5:
					return (6.01720935532251*mass + 0.63775031189316) / (mass + 3.1985473508617);
				case 6:
					return (7.01965692614781*mass + 0.70555519707735) / (mass + 3.13963062984322);
				default:
					return 0.0;
				}
			case 3:
				switch (idx) {
				case 0:
					return (1.00292210071565*mass + 0.69909138316283) / (mass + 1.21314257583115);
				case 1:
					return (2.00647743443781*mass - 2.04731818018074) / (mass + 2.28872295952506);
				case 2:
					return (3.00937441884854*mass - 1.28403970101672) / (mass + 3.28082670764372);
				case 3:
					return (4.01214578525685*mass + 0.245661364135073) / (mass + 3.98816635692297);
				case 4:
					return (5.01470441207058*mass + 0.401200861171199) / (mass + 4.03136684612641);
				case 5:
					return (6.01713997467128*mass + 0.520076859993294) / (mass + 4.00083574106646);
				case 6:
					return (7.01946371995006*mass + 0.616702302589006) / (mass + 3.92315219567947);
				default:
					return 0.0;
				}
			case 4:
				switch (idx) {
				case 0:
					return (1.00295179732532*mass + 0.516323424561711) / (mass + 1.21784374686008);
				case 1:
					return (2.0062133708259*mass - 1428.27219103666) / (mass - 709.232458840463);
				case 2:
					return (3.00946353863036*mass - 2.0503055749893) / (mass + 3.7542429320458);
				case 3:
					return (4.01218829322677*mass - 1.37444555453438) / (mass + 4.3177550301271);
				case 4:
					return (5.0146221731438*mass + 0.264623618433069) / (mass + 4.68440343056405);
				case 5:
					return (6.01696008583559*mass + 0.403521815710386) / (mass + 4.64632790640854);
				case 6:
					return (7.01917934163197*mass + 0.511953925392437) / (mass + 4.55969640768029);
				default:
					return 0.0;
				}
			default:
				return 0.0;
			}
		}
		return 0.0;
	}

	// DB에서 구한 원래 비율, 비율곱 범위, 질량 작은 경우
	private static double avgR1(int idx, double mass) {	//idx=0 앞에 miss peak이 하나있음  idx=1 misspeak 가 두개있음
		if (mass < Constants.T_MASS1) {
			return 0;
		}

		switch (idx) {
		case 0:
			return 0.000542415752193563*mass;
		case 1:
			return 0.00027125*mass + 0.0817923096774193;
		case 2:
			return (0.000180833333333333*mass*mass + 0.163421198156682*mass) / (mass + 301.23723162522);
		case 3:
			return (0.000135625*mass*mass*mass + 0.245131797235023*mass*mass + 36.9214119911966*mass + 743.606111889506) / (mass*mass + 904.615406570537*mass);
		case 4:
			return (0.0001085*mass*mass*mass*mass + 0.326842396313364*mass*mass*mass + 147.685647964786*mass*mass + 972.168014092919*mass) / (mass*mass*mass + 1807.42338975132*mass*mass + 272231.60915168*mass + 21953171.4064928);
		default:
			return 0;
		}
	}

	private static  double maxR1(int idx, double mass) {
		if (mass < Constants.T_MASS1) {
			return Constants.POSINF;
		}

		switch (idx) {
		case 0:
			return 0.000676239300871878*mass;
			//		return 0.000542415752193563*mass + 0.240882388;
			//		return 0.000565221085438767*mass + 0.164821233808595;
		case 1:
			return 0.000212750085327369*mass + 0.395175590871491;
		case 2:
			return (0.000210681739144816*mass*mass + 0.329016729075631*mass) / (mass + 301.294523673496);
		case 3:
			return (0.000161445207597051*mass*mass*mass + 0.166871551561399*mass*mass - 81.1104742069731*mass + 992.650069249618) / (mass*mass - 333.506395960027*mass);
		case 4:
			return (0.000127724276064344*mass*mass*mass*mass + 0.157365554050182*mass*mass*mass - 182.292929202306*mass*mass + 53800.0533800914*mass) / (mass*mass*mass - 590.066392879796*mass*mass - 962.712351142465*mass + 47535811.9487792);
		default:
			return 0;
		}

	}
	
	private static  double minR1(int idx, double mass) {
		if (mass < Constants.T_MASS1) {
			return Constants.NEGINF;
		}

		switch (idx) {
		case 0:
			return 0.00041149412414058*mass;
			//		return 0.000542415752193563*mass - 0.23565893;
			//		return 0.000519610418948359*mass - 0.155883125684508;
		case 1:
			return 0.000298048896263896*mass - 0.15013290934429;
		case 2:
			return (0.00017459437706253*mass*mass - 0.0198268692884906*mass) / (mass - 0.0198268692884906);
		case 3:
			return (0.00015174741917834*mass*mass*mass - 0.144477133760431*mass*mass + 28.4337348520276*mass + 2092.30324452523) / (mass*mass - 387.671789769344*mass);
		case 4:
			return (9.65617327501916e-05*mass*mass*mass*mass - 0.0469207552687451*mass*mass*mass - 75.3894322769834*mass*mass + 30985.3976564314*mass) / (mass*mass*mass - 488.334703879274*mass*mass + 20081.8835459434*mass + 8764389.73660916);
		default:
			return 0;
		}

	}

	
	private static  double avgRP1(int idx, double mass) {
		if (mass < Constants.T_MASS1) {
			return 0;
		}

		switch (idx) {
		case 0:
			return 0.511834778766633 + 113.655447911276 / mass;
		case 1:
			return 0.735889734963173 - 21.5652918316741 / (mass + 1934.1606296315);
		case 2:
			return 0.75562558469661 + 150.296854120924 / (mass + 1807.45158568401);
		case 3:
			return 0.823973968392781 + 25.7700994896072 / (mass - 77.0997141396938);
		default:
			return 0;
		}

	}
	
	private static  double maxRP1(int idx, double mass) {
		if (mass < Constants.T_MASS1) {
			return Constants.POSINF;
		}

		switch (idx) {
		case 0:
			return 0.407752590742579 + 783.603778796224 / mass;
		case 1:
			return 0.736604624444898 + 150.078546067355 / (mass - 51.5379499203973);
		case 2:
			return 0.722434726036907 + 365.212069675127 / (mass - 100.673006438236);
		case 3:
			return 0.322175486611715 + 1162.33684117372 / (mass - 162.769375072638);
		default:
			return 0;
		}

	}
	
	private static  double minRP1(int idx, double mass) {
		if (mass < Constants.T_MASS1) {
			return Constants.NEGINF;
		}

		switch (idx) {
		case 0:
			return 0.553607822946901 - 332.903209475669 / mass;
		case 1:
			return 0.737558485389733 - 162.46842805394 / (mass - 115.45975496947);
		case 2:
			return 0.872825854542323 - 361.31510810905 / (mass - 48.5894639880709);
		case 3:
			return 1.27799490706067 - 1080.27549645868 / (mass - 113.959973644233);
		default:
			return 0;
		}
	}

		// DB에서 구한 원래 비율, 비율곱 범위, 질량 중간 경우
	private static  double avgR2(int idx, double mass) {
		if (mass < Constants.T_MASS2) {
			return 0;
		}

		switch (idx) {
		case 0:
			return 0.000543595809751539*mass;
		case 1:
			return 0.000278886759073851*mass + 0.0593755682707423;
		case 2:
			return 0.00018616539663102*mass + 0.0741090194786749;
		case 3:
			return 0.000142636534821345*mass + 0.0718221265784284;
		case 4:
			return 0.000116671719376751*mass + 0.0674157215378921;
		case 5:
			return 9.91297857715994e-05*mass + 0.0637691457326418;
		case 6:
			return 8.6824087413594e-05*mass + 0.0593198002510331;
		default:
			return 0;
		}

	}
	
	private static  double maxR2(int idx, double mass) {
		if (mass < Constants.T_MASS2) {
			return 0;
		}

		switch (idx) {
		case 0:
			return 0.000656788437554653*mass;
		case 1:
			return 0.000293139943018318*mass + 0.224435336534894;
		case 2:
			return 0.000193944573048274*mass + 0.243258705422461;
		case 3:
			return 0.000145701702790161*mass + 0.245110728040748;
		case 4:
			return 0.000120471924198261*mass + 0.228652538033455;
		case 5:
			return 9.82480127793669e-05*mass + 0.23847678742629;
		case 6:
			return 7.02828087937578e-05*mass + 0.299470983579012;
		default:
			return 0;
		}

	}
	
	private static  double minR2(int idx, double mass) {
		if (mass < Constants.T_MASS2) {
			return 0;
		}

		switch (idx) {
		case 0:
			return 0.000433008682456966*mass;
		case 1:
			return 0.000245358559231721*mass - 0.0403204503971093;
		case 2:
			return 0.00017059058721391*mass - 0.0453494797851088;
		case 3:
			return 0.000132365050071466*mass - 0.050980579045075;
		case 4:
			return 0.000110631299550308*mass - 0.0627008308617091;
		case 5:
			return 0.000101204060652847*mass - 0.0966918302086797;
		case 6:
			return 0.000107958655292456*mass - 0.184954335888997;
		default:
			return 0;
		}

	}

	
	private static  double avgRP2(int idx, double mass) {
		if (mass < Constants.T_MASS2) {
			return 0;
		}

		switch (idx) {
		case 0:
			return 0.503245602637237 + 129.125338925119 / mass;
		case 1:
			return 0.665832736897384 + 231.16379922434 / (mass + 2062.611998139);
		case 2:
			return 0.758771780232274 + 150.530364451505 / (mass + 1807.4310188011);
		case 3:
			return 0.816721604267822 + 73.9230055844914 / (mass + 966.467555305382);
		case 4:
			return 0.849528126443098 + 53.5360327282436 / (mass + 312.299788149686);
		case 5:
			return 0.872011305090549 + 49.3969441138668 / (mass + 89.6675445851847);
		case 6:
			return 0.887140619837095 + 31.9081710443378 / (mass - 1528.85159380718);
		default:
			return 0;
		}

	}
	
	private static  double maxRP2(int idx, double mass) {
		if (mass < Constants.T_MASS2) {
			return Constants.POSINF;
		}

		switch (idx) {
		case 0:
			return 0.507672931038912 + 545.507772005401 / mass;
		case 1:
			return 0.731678227122246 + 380.374840226825 / (mass + 1853.77725850901);
		case 2:
			return 0.824834840599112 + 219.756975349203 / (mass - 16.5263909709524);
		case 3:
			return 0.86223629835813 + 277.895690728042 / (mass - 741.863306487537);
		case 4:
			return 0.848737456991001 + 538.882009064593 / (mass - 926.114948283104);
		case 5:
			return 0.798166653681199 + 950.23439327722 / (mass - 1180.38864177168);
		case 6:
			return 0.617933668390925 + 1792.67686764809 / (mass - 1433.11134010978);
		default:
			return 0;
		}

	}
	
	private static  double minRP2(int idx, double mass) {
		if (mass < Constants.T_MASS2) {
			return Constants.NEGINF;
		}

		switch (idx) {
		case 0:
			return 0.477206733750061 - 151.940839267185 / mass;
		case 1:
			return 0.60889965165977 + 40.1296575225364 / (mass + 2402.86189363569);
		case 2:
			return 0.70496680034191 - 40.7169197391889 / (mass - 1079.58316897807);
		case 3:
			return 0.776106440733451 - 181.480055670946 / (mass - 946.548040748999);
		case 4:
			return 0.856343729650243 - 478.778616744075 / (mass - 917.443449320424);
		case 5:
			return 0.954063111383699 - 911.36008179514 / (mass - 1113.55910634638);
		case 6:
			return 1.11114026642057 - 1628.43895088082 / (mass - 1342.82592799463);
		default:
			return 0;
		}
	}

		// DB에서 구한 원래 비율, 비율곱 범위, 질량 큰 경우
	private static  double avgR3(int idx, double mass) {
		if (mass < Constants.T_MASS3) {
			return 0;
		}

		switch (idx) {
		case 0:
			return 0.00054755009064469 *mass + 1.001e-30;
		case 1:
			return 0.000276648501600684*mass + 0.0555469243485302;
		case 2:
			return 0.000185273724175025*mass + 0.0782028105854292;
		case 3:
			return 0.000139957881463058*mass + 0.084364128158295;
		case 4:
			return 0.000112963445181452*mass + 0.0847186536097926;
		case 5:
			return 9.49946004647695e-05*mass + 0.0834046798928754;
		case 6:
			return 8.22950358601902e-05*mass + 0.080553780647878;
		case 7:
			return 7.27177475060235e-05*mass + 0.0779919065508693;
		case 8:
			return 6.45462284105087e-05*mass + 0.0814692833459906;
		case 9:
			return 5.92312799517667e-05*mass + 0.0736250327801014;
		case 10:
			return 5.3586471443843e-05*mass + 0.0777919536474791;
		case 11:
			return 4.94313780163201e-05*mass + 0.0763190813133117;
		case 12:
			return 4.54910475093714e-05*mass + 0.0787451515118554;
		case 13:
			return 4.09628815622464e-05*mass + 0.0911521684923799;
		default:
			return 0;
		}

	}
	
	private static  double maxR3(int idx, double mass) {
		if (mass < Constants.T_MASS3) {
			return 0;
		}

		switch (idx) {
		case 0:
			return 0.000723114855391446*mass + 1.001e-30;
		case 1:
			return 0.000363595652519864*mass - 0.120052571475427;
		case 2:
			return 0.000223910383023828*mass + 0.0916721110261479;
		case 3:
			return 0.000161302792805347*mass + 0.166388568760553;
		case 4:
			return 0.000127434496451832*mass + 0.188516637654909;
		case 5:
			return 0.00010534430897877*mass + 0.201270973319086;
		case 6:
			return 8.97006567498902e-05*mass + 0.20917612335033;
		case 7:
			return 7.66834967469582e-05*mass + 0.227185459201798;
		case 8:
			return 6.94419323111995e-05*mass + 0.218751804620972;
		case 9:
			return 6.22776028965714e-05*mass + 0.222888675513929;
		case 10:
			return 5.3586471443843e-05*mass + 0.278991236695385;
		case 11:
			return 4.94313780163201e-05*mass + 0.25841100765211;
		case 12:
			return 4.54910475093714e-05*mass + 0.26725660229666;
		case 13:
			return 4.09628815622464e-05*mass + 0.304690897164158;
		default:
			return 0;
		}

	}
	
	private static  double minR3(int idx, double mass) {
		if (mass < Constants.T_MASS3) {
			return 0;
		}

		switch (idx) {
		case 0:
			return 0.000380162468528693*mass + 1.001e-30;
		case 1:
			return 0.000188956406213499*mass + 0.240140521205437;
		case 2:
			return 0.000141903740345654*mass + 0.0989740414518865;
		case 3:
			return 0.000113295041163779*mass + 0.0449636542114608;
		case 4:
			return 9.39667925372497e-05*mass + 0.0213551492083562;
		case 5:
			return 8.13775320919289e-05*mass - 0.00073171152409789;
		case 6:
			return 7.28126474514515e-05*mass - 0.0219128936834298;
		case 7:
			return 6.80017847852428e-05*mass - 0.054905744583198;
		case 8:
			return 5.94130703056629e-05*mass - 0.0437513731601505;
		case 9:
			return 5.46720340352568e-05*mass - 0.0538423709493245;
		case 10:
			return 5.3586471443843e-05*mass - 0.125167796604039;
		case 11:
			return 4.94313780163201e-05*mass - 0.097761230788494;
		case 12:
			return 4.54910475093714e-05*mass - 0.104762706615447;
		case 13:
			return 4.09628815622464e-05*mass - 0.118988297132509;
		default:
			return 0;
		}

	}

	
	private static  double avgRP3(int idx, double mass) {
		if (mass < Constants.T_MASS3) {
			return 0;
		}

		switch (idx) {
		case 0:
			return 0.568316339968907 - 1.04849096806338e-05*mass + 5.29583465838217e-10*mass*mass;
		case 1:
			return 0.733460256217147 - 9.18658035762221e-06*mass + 4.37558731480206e-10*mass*mass;
		case 2:
			return 0.807432670246971 - 6.85989338754518e-06*mass + 3.00321107090024e-10*mass*mass;
		case 3:
			return 0.848792400341667 - 5.14537679442748e-06*mass + 2.10372854185365e-10*mass*mass;
		case 4:
			return 0.876992946428059 - 4.3071501285348e-06*mass + 1.69700038175106e-10*mass*mass;
		case 5:
			return 0.89636433825948 - 3.76396347682756e-06*mass + 1.56661804621463e-10*mass*mass;
		case 6:
			return 0.910752542118436 - 3.07313102723661e-06*mass + 1.08929701630143e-10*mass*mass;
		case 7:
			return 0.93714576129582 - 6.83631892958637e-06*mass + 3.64168663718938e-10*mass*mass;
		case 8:
			return 0.975519353370565 - 1.38758467603088e-05*mass + 7.99504253392897e-10*mass*mass;
		case 9:
			return 0.971648408579592 - 9.6833411343062e-06*mass + 4.91036874403108e-10*mass*mass;
		case 10:
			return 1.01828916680805 - 1.81953781058948e-05*mass + 9.49446803287465e-10*mass*mass;
		case 11:
			return 1.22464703658373 - 6.30297296836898e-05*mass + 3.45920274200409e-09*mass*mass;
		case 12:
			return 1.95099166705112 - 0.00022690814544934*mass + 1.27469499988468e-08*mass*mass;
		default:
			return 0;
		}

	}
	
	private static  double maxRP3(int idx, double mass) {
		if (mass < Constants.T_MASS3) {
			return Constants.POSINF;
		}

		switch (idx) {
		case 0:
			return 0.831789180281756 - 7.58841159585501e-05*mass + 6.96175358036141e-09*mass*mass;
		case 1:
			return 0.911006329452626 - 4.37405672538403e-05*mass + 3.85614214011361e-09*mass*mass;
		case 2:
			return 0.955433550006254 - 2.82605144216244e-05*mass + 2.28449217614987e-09*mass*mass;
		case 3:
			return 1.05368830283444 - 3.70698615821216e-05*mass + 2.50662767003313e-09*mass*mass;
		case 4:
			return 1.20062429798949 - 6.01943236228298e-05*mass + 3.62831083577672e-09*mass*mass;
		case 5:
			return 1.44399802858061 - 0.000106730429512693*mass + 6.12437216302458e-09*mass*mass;
		case 6:
			return 1.91072149779034 - 0.000206131813363327*mass + 1.17494847573272e-08*mass*mass;
		case 7:
			return 2.79046256024063 - 0.000403308332304838*mass + 2.32120785677021e-08*mass*mass;
		case 8:
			return 4.65846765998948 - 0.000839070280470542*mass + 4.90063556168516e-08*mass*mass;
		case 9:
			return 4.18096532811616 - 0.000629081651027411*mass + 3.2463063725138e-08*mass*mass;
		case 10:
			return 6.88267435457757 - 0.00118183867657941*mass + 6.16138421163546e-08*mass*mass;
		case 11:
			return 14.1170311988883 - 0.00274223529011024*mass + 1.47048870410186e-07*mass*mass;
		case 12:
			return 33.7513537151654 - 0.00709538770527136*mass + 3.89856998765122e-07*mass*mass;
		default:
			return 0;
		}

	}
	
	private static  double minRP3(int idx, double mass) {
		if (mass < Constants.T_MASS3) {
			return Constants.NEGINF;
		}

		switch (idx) {
		case 0:
			return 0.351772029522468 + 4.24661648161992e-05*mass - 5.01100223965691e-09*mass*mass;
		case 1:
			return 0.569044013675839 + 2.22700172959916e-05*mass - 2.750217995896e-09*mass*mass;
		case 2:
			return 0.659298266157043 + 1.51271254590081e-05*mass - 1.71291415520162e-09*mass*mass;
		case 3:
			return 0.64558089531114 + 2.66278677299093e-05*mass - 2.06949282469354e-09*mass*mass;
		case 4:
			return 0.558565659037797 + 5.05135007239057e-05*mass - 3.21245117741485e-09*mass*mass;
		case 5:
			return 0.355579394158034 + 9.78999844395257e-05*mass - 5.73494844723805e-09*mass*mass;
		case 6:
			return -0.0391476292908423 + 0.00018794679683782*mass - 1.08025486966718e-08*mass*mass;
		case 7:
			return -0.768092961597599 + 0.000353114455127164*mass - 2.02372437621952e-08*mass*mass;
		case 8:
			return -2.28535389754734 + 0.000706525395348056*mass - 4.09716723699235e-08*mass*mass;
		case 9:
			return -1.98870581260309 + 0.000555831475738199*mass - 2.85191205080274e-08*mass*mass;
		case 10:
			return -4.08003922853999 + 0.000978126099246765*mass - 5.05006331331271e-08*mass*mass;
		case 11:
			return -8.88988096274267 + 0.00199187524209113*mass - 1.04986006885241e-07*mass*mass;
		case 12:
			return -20.2621512269365 + 0.00446301323883152*mass - 2.40891206944064e-07*mass*mass;
		default:
			return 0;
		}

	}
	
	private static  int fac(int n) {
		if (n <= 1) return 1;
		return n*fac(n - 1);
	}
	
	private static  double calcIn(int k, double c, double h, double n, double o, double s) {
		double in = 0.0;
		for (int k4 = 0; k4 <= k / 4; k4++) {
			for (int k2 = 0; k2 * 2 + k4 * 4 <= k; k2++) {
				int k1 = k - 4 * k4 - 2 * k2;
				double t1 = c*0.011070 / 0.988930 + h*0.00015 / 0.99985 + n*0.003663 / 0.996337 + o*0.000374 / 0.997590 + s*0.0075 / 0.9502;
				double t2 = o*0.002036 / 0.997590 + s*0.0421 / 0.9502;
				double t4 = s*0.0002 / 0.9502;
				in += Math.pow(t1, k1) * Math.pow(t2, k2) * Math.pow(t4, k4) / fac(k1) / fac(k2) / fac(k4);
			}
		}
		return in;
	}

}
