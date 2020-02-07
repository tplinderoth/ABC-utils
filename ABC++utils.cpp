/*
 * ABC++utils.cpp
 *
 */

/*
 * functions to implement
 *
 * 1) rescale summary stats to range 0 to 1: z_i = [x_i - min(x)] / [max(x) - min(x)]
 *
 * 2) distance rejection function - in ABCToolbox, which gives choice of standardizing stats prior to rejection
 * weighted euclidean distance: d(s, s_obs) = [ sum_i=1_m ( (s_i - s_obs,i)/delta_i )^2 ]^(1/2)
 * delta_i = empirical standard deviation of the simulated s_i values (Prangle 2017)
 *
 */

#include "ABC++utils.h"
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>

AlleleFreq::AlleleFreq ()
	: allele(4),
      freq(0)
{}

PopSeq::PopSeq (const char* filename)
	: npops(0),
	  nsnp(0)
{
	// preallocate some memory to be a bit faster
	int pre_n = 2000;
	popmap.reserve(pre_n);
	seq.reserve(pre_n);
	sampsz.reserve(4);

	if (filename) {
		parseArp(filename);
	}
}

unsigned int PopSeq::nseq () {
	return seq.size();
}

int PopSeq::popn () {
	return npops;
}

unsigned int PopSeq::snpn () {
	return nsnp;
}

int PopSeq::parseArp (const char* filename) {

	std::ifstream arpfs(filename);
	std::string arpline;
	int popn = -1;
	int nseq_total = 0;
	int n = 0;

	while (getline(arpfs, arpline)) {

		if (nsnp == 0 && arpline.find("#Total number of polymorphic sites:") != std::string::npos) {
			// find number of SNPs
			nsnp = stoi(arpline.substr(arpline.rfind(": ") + 2));
		} else {
			if (arpline[0] != '#' && arpline.find("SampleSize") != std::string::npos) {
				// get sequences
				int nseq = stoi(arpline.substr(arpline.rfind("=")+1));
				sampsz.push_back(nseq);
				nseq_total += nseq;
				++popn;
				arpfs.ignore(15, '\n');
				for (int i=0; i<nseq; ++i) {
					getline(arpfs, arpline);
					seq.push_back(arpline.substr(arpline.length()-nsnp));
					popmap.push_back(popn);
					++n;
				}
			}
		}
	}

	npops = popn+1;

	return n;
}

int popInput (int argc, int i, char** argv, std::vector<int*> &pops, int maxn, std::vector<std::string> &popid, int exclude = -1) {
	int npops = 0;

	while (i < argc) {
		if (argv[i][0] != '-') {
			int* vals = new int[maxn];
			char* p = strtok(argv[i], ",");
			int n = 0;
			while (p != NULL) {
				if (n >= maxn) return -1;
				vals[n] = atoi(p);
				if (n != exclude) {
					popid.push_back(p);
				}
				p = strtok(NULL, ",");
				++n;
			}
			pops.push_back(vals);
			++npops;
		} else {
			break;
		}
		++i;
	}

	return npops > 0 ? i:-1;
}

void statsInfo () {
	int w=12;

	std::cerr << "\nCalculate summary statistics from Fastsimcoal2 DNA sequence\n";

	std::cerr << "\nABC++utils doStats [ARGUMENTS]\n";

	std::cerr << "\nARGUMENTS:\n\n"
	<< std::setw(w) << std::left << "-arpfile" << "Name of Arlequin style DNA sequence file written by fastsimcoal2\n"
	<< std::setw(w) << std::left << "-out" << "Name of output file, default: standard output\n"
	<< std::setw(w) << std::left << "-segsites" << "Number segregating sites, takes list of population IDs\n"
	<< std::setw(w) << std::left << "-sfs" << "SFS bins for a list of <number derived alleles,pop ID>\n"
	<< std::setw(w) << std::left << "-pi" << "Average number of pairwise sequence differences, takes list of population IDs\n"
	<< std::setw(w) << std::left << "-thetaw" << "Watterson's theta, takes list of population IDs\n"
	<< std::setw(w) << std::left << "-tajimaD" << "Tajima's D, takes list of population IDs\n"
	<< std::setw(w) << std::left << "-doubleton" << "Probability of sharing doubleton alleles, takes list of <1st allele pop ID,2nd allele pop ID>\n"
	<< std::setw(w) << std::left << "-fst" << "Pairwise Fst, takes list of <pop1 ID,pop2 ID>\n"
	<< std::setw(w) << std::left << "-fold" << "Folds the frequency spectrum, default: unfolded\n"
	<< std::setw(w) << std::left << "-sfsprob" << "Convert SFS bin counts to probabilities, default: counts\n"
	<< std::setw(w) << std::left << "-seqlen" << "An INT specifying the total sequence length to convert stats to per-site estimates\n"
	<< std::setw(w) << std::left << "-noheader" << "Suppress writing the header to output, default: header is written\n"
	<< std::setw(w) << std::left << "-verbose" << "Print all messages to screen\n";

	std::cerr << "\nNOTES:\n"
	<< "* Input to arguments requiring lists should be space-delimited, e.g. -segsites 0 1 2 -doubleton 0,0 0,1 1,1\n"
	<< "* Statistics are written in the same order as they are called in the argument list\n\n";
}

int statArgs (int argc, char** argv, std::string &arpname, std::ofstream &outfstream, std::ostream &os, std::vector< std::vector<int*> > &pops, std::vector<std::string> &popid,
		std::vector<int> &writeorder, int &fold, int &sfsprob, unsigned int &seqlen, int &writeheader, int &verbose) {

	if (argc < 6 || (argc > 2 && strcmp(argv[2],"-help") == 0) || (argc > 2 && strcmp(argv[2], "-h") == 0)) {
		statsInfo();
		return 1;
	}

	int i = 2;
	unsigned int j = 0;
	while (i < argc) {
		if (strcmp(argv[i], "-doubleton") == 0) {
			if ((i = popInput(argc, i+1, argv, pops[0], 2, popid)-2) < 0) {
				std::cerr << "-doubleton needs space-delimited list of '<first allele pop ID,second allele pop ID>'\n";
				return -1;
			}
			writeorder.push_back(0);
		} else if (strcmp(argv[i], "-pi") == 0) {
			if ((i = popInput(argc, i+1, argv, pops[1], 1, popid)-2) < 0) {
				std::cerr << "-pi needs a space-delimited list of single population IDs\n";
				return -1;
			}
			writeorder.push_back(1);
		} else if (strcmp(argv[i], "-thetaw") == 0) {
			if ((i = popInput(argc, i+1, argv, pops[2], 1, popid)-2) < 0) {
				std::cerr << "-thetaw needs a space-delimited list of single population IDs\n";
				return -1;
			}
			writeorder.push_back(2);
		} else if (strcmp(argv[i], "-tajimaD") == 0) {
			if ((i = popInput(argc, i+1, argv, pops[3], 1, popid)-2) < 0) {
				std::cerr << "-tajimaD needs a space-delimited list of single population IDs\n";
				return -1;
			}
			writeorder.push_back(3);
		} else if (strcmp(argv[i], "-sfs") == 0) {
			if((i = popInput(argc, i+1, argv, pops[4], 2, popid, 0)-2) < 0) {
				std::cerr << "-sfs needs space-delimited list of '<sfs category>,<population ID>'\n";
				return -1;
			}
			writeorder.push_back(4);
		} else if (strcmp(argv[i], "-segsites") == 0) {
			if ((i = popInput(argc, i+1, argv, pops[5], 1, popid)-2) < 0) {
				std::cerr << "-segsites needs a space-delimited list of single population IDs\n";
				return -1;
			}
			writeorder.push_back(5);
		} else if (strcmp(argv[i], "-fst") == 0) {
			if ((i = popInput(argc, i+1, argv, pops[6], 2, popid)-2) < 0) {
				std::cerr << "-fst needs a space-delimited list of '<pop1 ID>,<pop2 ID>'\n";
				return -1;
			}
			for (j=0; j<pops[6].size(); ++j) {
				if (pops[6][j][0] == pops[6][j][1]) {
					std::cerr << "Population IDs for pairwise Fst calculations must be different\n";
					return -1;
				}
			}
			writeorder.push_back(6);
		} else if (strcmp(argv[i], "-arpfile") == 0) {
			arpname = argv[i+1];
		} else if (strcmp(argv[i], "-out") == 0) {
			outfstream.open(argv[i+1]);
			if (!outfstream) {
				std::cerr << "Problem opening output file " << argv[i+1] << "\n";
				return -1;
			}
		} else if (strcmp(argv[i], "-sfsprob") == 0) {
			sfsprob = 1;
			--i;
		} else if (strcmp(argv[i], "-fold") == 0) {
			fold = 1;
			--i;
		} else if (strcmp(argv[i], "-seqlen") == 0) {
			seqlen = atoi(argv[i+1]);
			if (seqlen <= 0) {
				std::cerr << "-seqlen needs a sequence length greater than zero\n";
				return -1;
			}
		} else if (strcmp(argv[i], "-noheader") == 0) {
			writeheader = 0;
			--i;
		} else if (strcmp(argv[i], "-verbose") == 0) {
			verbose = 1;
			--i;
		} else {
			std::cerr << "Unknown argument " << argv[i] << "\n";
			return -1;
		}

		i += 2;
	}

	if (outfstream.is_open()) {
		os.rdbuf(outfstream.rdbuf());
	}

	return 0;
}

std::string checkPopIDs (std::vector<std::string> &id, int npops) {

	std::vector<std::string>::const_iterator it;
	for (it = id.begin(); it != id.end(); ++it) {
		int i = 0;
		std::string pos = std::to_string(i);
		while (pos != *it) {
			++i;
			if (i >= npops) {
				return *it;
			}
			pos = std::to_string(i);
		}
	}

	return "";
}

int mafidx (int* countarr, int n) {
	int maxidx = 1;
	int mafidx = 0;

	if (n > 1) {
		if (countarr[mafidx] > countarr[maxidx]) {
			maxidx = 0;
			mafidx = 1;
		}
	}

	if (n > 2) {
		for (int i=2; i<n; ++i) {
			if (countarr[i] > countarr[mafidx]) mafidx = i;
			if (countarr[mafidx] > countarr[maxidx]) {
				mafidx = maxidx;
				maxidx = i;
			}
		}
	}

	return mafidx;
}

int calcSFS (std::vector< std::vector<AlleleFreq> > &af, std::vector< std::vector<int> > &sfs, unsigned int** double_counts, int nsnp) {
	// consider total population minor allele derived

	unsigned int totalpop_idx = af.size()-1;
	unsigned int pop, i;
	int n [totalpop_idx];
	for (unsigned int k=0; k<totalpop_idx; ++k) {
		n[k] = sfs[k].size()-1;
	}
	int doublepops [2];

	for (int snp = 0; snp < nsnp; ++snp) {
		int didx = 0;

		for (pop = 0; pop < totalpop_idx; ++pop) {
			if (af[pop][snp].allele == af[totalpop_idx][snp].allele) {
				// population minor allele is the same as total population minor, i.e. derived
				i = af[pop][snp].freq;
			} else {
				// population minor allele different from total population minor, i.e. ancestral
				if (af[pop][snp].freq > 0) {
					i = n[pop] - af[pop][snp].freq;
				} else {
					// fix allele identity so that minor allele is total population minor allele
					af[pop][snp].allele = af[totalpop_idx][snp].allele;
					i = 0;
				}
			}

			++sfs[pop][i];

			// record populations constituting a doubleton site
			if (af[totalpop_idx][snp].freq == 2) {
				if (didx > 2) {
					std::cerr << "More derived alleles found among subpopulations than in the total population\n";
					return -1;
				}
				for (unsigned int p=0; p<i; ++p) {
					doublepops[didx] = pop;
					++didx;
				}
			}

		}

		// iterate doubleton counts
		if (af[totalpop_idx][snp].freq == 2) {
			++double_counts[doublepops[0]][doublepops[1]];
		}
	}

	return 0;
}

int probDoubleShare (unsigned int** double_counts, int npops, const std::vector<int*> &pops, std::vector<double>* prob) {

	prob->resize(pops.size());

	int maxidx = npops-1;
	double total_counts = 0;

	for (int i=0; i<npops; ++i) {
		for (int j=0; j<npops; ++j) {
			total_counts += double_counts[i][j];
		}
	}

	for (unsigned int i=0; i<pops.size(); ++i) {
		int p1 = pops[i][0];
		int p2 = pops[i][1];

		if (p1 > maxidx || p2 > maxidx) {
			std::cerr << "Max population ID is " << maxidx << " for calculating doubleton allele sharing\n";
			return -1;
		}

		if (p1 == p2) {
			(*prob)[i] = double_counts[p1][p2] / total_counts;
		} else {
			(*prob)[i] = (double_counts[p1][p2] + double_counts[p2][p1]) / total_counts;
		}
	}

	return 0;
}

void calcSegSites (const std::vector< std::vector<int> > &sfs, std::vector<int> &s) {
	// this will calculate the number of SNPs in each population
	s.resize(sfs.size(), 0);

	for (unsigned int i=0; i<sfs.size(); ++i) {
		for (unsigned int j=1; j < sfs[i].size()-1; ++j) {
			s[i] += sfs[i][j];
		}
	}
}

unsigned int binomCoeff (int n, int k) {
	double res = 1;

	// C(n,k) = C(n,n-k)
	if (k > n-k)
		k = n-k;

	// calculate [n*(n-1) *...* (n-k+1)] / [k * (k-1) *...* 1]
	for (int i=0; i<k; ++i) {
		res *= (n-i);
		res /= (i+1);
	}

	return res;
}

void calcWatterson (const std::vector<int> &s, const std::vector<int> &n, std::vector<double> &watterson) {
	// s = number segregating sites
	// n = haploid sample size

	watterson.resize(s.size());

	for (unsigned int i=0; i<s.size(); ++i) {
		double a1 = 0;
		for (int j=1; j<n[i]; ++j) {
			a1 += 1.0/j;
		}

		watterson[i] = s[i]/a1;
	}
}

void calcPi (const std::vector< std::vector<int> > &sfs, std::vector<double> &pi) {

	pi.resize(sfs.size());

	for (unsigned int i=0; i<sfs.size(); ++i) {
		int n = sfs[i].size()-1;
		double ndiff = 0;

		for (int j=1; j<n; ++j) {
			ndiff += j*(n-j)*sfs[i][j];
		}

		pi[i] = ndiff/binomCoeff(n, 2);
	}
}

int calcTajimaD(const std::vector< std::vector<int> > &sfs, const std::vector<int> &s, const std::vector<int> &n, const std::vector<int*> &pops,
		std::vector<double> &watterson, std::vector<double> &pi, std::vector<double> &tajimad) {

	tajimad.resize(pops.size());

	calcWatterson(s, n, watterson);
	calcPi(sfs, pi);

	for (unsigned int i=0; i<pops.size(); ++i) {
		int pop = pops[i][0];
		if (s[pop] < 1) {
			std::cerr << "Population " << pop << " has no variable sites - can't calculate Tajima's D\n";
			return -1;
		}
		double nseq = n[pop];
		double a1 = 0;
		double a2 = 0;

		for (int j=1; j<nseq; ++j) {
			a1 += 1.0/j;
			a2 += 1.0/(j*j);
		}

		double b1 = (nseq+1)/(3*(nseq-1));
		double b2 = (2*(nseq*nseq + nseq + 3))/(9*nseq*(nseq-1));
		double c1 = b1 - 1/a1;
		double c2 = b2 - (nseq + 2)/(a1*nseq) + a2/(a1*a1);
		double e1 = c1/a1;
		double e2 = c2/(a1*a1 + a2);
		double dvar = sqrt(e1*s[pop] + e2*s[pop]*(s[pop]-1));

		tajimad[i] = (pi[pop] - watterson[pop])/dvar;

		/*
		std::cerr << "\n-----" << pop << "-----\n";
		std::cerr << "n: " << nseq << "\n";
		std::cerr << "a1: " << a1 << "\n";
		std::cerr << "a2: " << a2 << "\n";
		std::cerr << "b1: " << b1 << "\n";
		std::cerr << "b2: " << b2 << "\n";
		std::cerr << "c1: " << c1 << "\n";
		std::cerr << "c2: " << c2 << "\n";
		std::cerr << "e1: " << e1 << "\n";
		std::cerr << "e2: " << e2 << "\n";
		std::cerr << "dvar: " << dvar << "\n";
		std::cerr << "pi: " << pi[pop] << "\n";
		std::cerr << "watterson: " << watterson[pop] << "\n";
		std::cerr << "D: " << tajimad[i] << "\n";
		std::cerr << "------------\n";
		*/
	}

	return 0;
}

int sfsBins(const std::vector< std::vector<int> > &sfs, const std::vector<int*> &pops, int fold, const std::vector<int> &s, int doprobs, std::vector<double> &bins) {

	// consider total population minor allele derived

	bins.resize(pops.size());

	for (unsigned int i=0; i<pops.size(); ++i) {
		int pop = pops[i][1];
		int nderived = pops[i][0];
		int n = sfs[pop].size()-1;

		if (nderived > n) {
			std::cerr << "SFS category " << nderived << " exceeds population " << pop << " sample size\n";
			return -1;
		}

		bins[i] = sfs[pop][nderived];
		if (fold) {
			bins[i] += sfs[pop][n-nderived];
		}
		if (doprobs) {
			bins[i] /= s[pop];
		}
	}

	return 0;
}

void globalFst (const std::vector< std::vector<AlleleFreq> > &af, const std::vector<int> &n, const std::vector<int*> &pops, std::vector<double> &fst, unsigned int start = 0, unsigned int end = 0) {
	// calculate Reynolds, Weir, Cockerham 1983 Fst estimator for 2 populations assuming biallelic SNP
	// implemented using formula from Fumagalli et al. 2013

	// start: index of starting SNP for region to do calculation over
	// end: index of end SNP for region to do calculation over
	int totalidx = af.size()-1;
	int npair = pops.size();
	int pair, pop;

	double* varcomp [npair];
	for (pair=0; pair<npair; ++pair) {
		varcomp[pair] = new double [2];
		varcomp[pair][0] = 0; // a component
		varcomp[pair][1] = 0; // a+b component
	}

	// calculate terms that don't depend on allele frequency
	double* nhap [npair];
	double* ndip [npair];
	double npoolhap [npair];
	double npooldip [npair];
	double d1[npair];
	double d2[npair];

	for (pair=0; pair<npair; ++pair) {
		nhap[pair] = new double [2];
		ndip[pair] = new double [2];
		npoolhap[pair] = 0;
		npooldip[pair] = 0;
		for (pop=0; pop<2; ++pop) {
			nhap[pair][pop] = n[pops[pair][pop]];
			ndip[pair][pop] = 0.5 * nhap[pair][pop];
			npoolhap[pair] += nhap[pair][pop];
			npooldip[pair] += ndip[pair][pop];
		}
		d1[pair] = 2.0 * ((2.0*ndip[pair][0]*ndip[pair][1])/npooldip[pair]);
		d2[pair] = npooldip[pair] - 1.0;
	}

	double f [3] = {0}; // [0] = pop1 freq, [1] = pop2 freq, [2] = pooled pop freq

	// iterate over SNPs and calculate variance components

	end = end == 0 ? af[pops[0][0]].size():end+1;

	for (unsigned int site=start; site<end; ++site) {

		for (pair=0; pair<npair; ++pair) {

			f[2] = 0;

			for (pop = 0; pop<2; ++pop) {
				if (af[pops[pair][pop]][site].allele == af[totalidx][site].allele) {
					f[pop] = af[pops[pair][pop]][site].freq;
					f[2] += f[pop];
					f[pop] /= nhap[pair][pop];
				} else {
					f[pop] = nhap[pair][pop] - af[pops[pair][pop]][site].freq;
					f[2] += f[pop];
					f[pop] /= nhap[pair][pop];
				}
			}

			f[2] /= npoolhap[pair];

			double fdiff0 = f[0]-f[2];
			double fdiff1 = f[1]-f[2];
			double c0 = 4.0 * ndip[pair][0] * fdiff0 * fdiff0;
			double c1 = 4.0 * ndip[pair][1] * fdiff1 * fdiff1;
			double alpha0 = 2.0 * f[0] * (1.0 - f[0]);
			double alpha1 = 2.0 * f[1] * (1.0 - f[1]);
			double bs = (ndip[pair][0]*alpha0 + ndip[pair][1]*alpha1) / d2[pair]; // between pop variance
			double as = (c0 + c1 - bs)/d1[pair]; // total pop variance

			varcomp[pair][0] += as;
			varcomp[pair][1] += as + bs;
		}

	}

	// calculat fst for region
	fst.resize(npair);

	for (pair=0; pair<npair; ++pair) {
		fst[pair] = varcomp[pair][1] > 0 ? varcomp[pair][0]/varcomp[pair][1] : 0.0;
		if (fst[pair] < 0)
			fst[pair] = 0.0;
	}

	// free memory
	for (pair=0; pair<npair; ++pair) {
		delete [] varcomp[pair];
		delete [] nhap[pair];
		delete [] ndip[pair];
	}

}


void writeStats (std::ostream &outstream, const std::vector< std::vector<int*> > &pops, const std::vector<int> &printorder, int writeheader, const std::vector<int> &segsites,
		const std::vector<double> &pvar, const std::vector<double> &sfsbin, const std::vector<double> &doubleprob, const std::vector<double> &thetaw,
		const std::vector<double> &pi, const std::vector<double> &tajimad, const std::vector<double> &fst) {

	unsigned int j, k;
	std::string whitespace (" \t");
	std::size_t lastchar;

	// print header info

	if (writeheader) {
		std::stringstream hs;
		for (k=0; k<printorder.size(); ++k) {

			int i = printorder[k];

			switch (i) {
				case 0:
					if (!doubleprob.empty()) {
						for (j=0; j<pops[i].size(); ++j) {
							hs << "double" << pops[i][j][0] << "_" << pops[i][j][1] << "\t";
						}
					}
					break;
				case 1:
					if (!pi.empty()) {
						for (j=0; j<pops[i].size(); ++j) {
							hs << "pi" << pops[i][j][0] << "\t";
						}
					}
					break;
				case 2:
					if (!thetaw.empty()) {
						for (j=0; j<pops[i].size(); ++j) {
							hs << "thetaW" << pops[i][j][0] << "\t";
						}
					}
					break;
				case 3:
					if (!tajimad.empty()) {
						for (j=0; j<pops[i].size(); ++j) {
							hs << "tajimaD" << pops[i][j][0] << "\t";
						}
					}
					break;
				case 4:
					if (!sfsbin.empty()) {
						for (j=0; j<pops[i].size(); ++j) {
							hs << "sfs" << pops[i][j][0] << "_" << pops[i][j][1] << "\t";
						}
					}
					break;
				case 5:
					if (!segsites.empty()) {
						for (j=0; j<pops[i].size(); ++j) {
							hs << "S" << pops[i][j][0] << "\t";
						}
					}
					break;
				case 6:
					if (!fst.empty()) {
						for (j=0; j<pops[i].size(); ++j) {
							hs << "fst" << pops[i][j][0] << "_" << pops[i][j][1] << "\t";
						}
					}
			}
		}
		std::string header = hs.str();
		if ((lastchar = header.find_last_not_of(whitespace)) != std::string::npos) {
			header.erase(lastchar+1);
		}
		outstream << header << "\n";

	}

	// print stats

	std::stringstream ss;

	for (k=0; k<printorder.size(); ++k) {

		int i = printorder[k];

		switch (i) {
			case 0:
				if (!doubleprob.empty()) {
					for (j=0; j<doubleprob.size(); ++j) {
						ss << doubleprob[j] << "\t";
					}
				}
				break;
			case 1:
				if (!pi.empty()) {
					for (j=0; j<pi.size(); ++j) {
						ss << pi[j] << "\t";
					}
				}
				break;
			case 2:
				if (!thetaw.empty()) {
					for (j=0; j<thetaw.size(); ++j) {
						ss << thetaw[j] << "\t";
					}
				}
				break;
			case 3:
				if (!tajimad.empty()) {
					for (j=0; j<tajimad.size(); ++j) {
						ss << tajimad[j] << "\t";
					}
				}
				break;
			case 4:
				if (!sfsbin.empty()) {
					for (j=0; j<sfsbin.size(); ++j) {
						ss << sfsbin[j] << "\t";
					}
				}
				break;
			case 5:
				if (!segsites.empty()) {
					for (j=0; j<segsites.size(); ++j) {
						if (pvar.size() == segsites.size()) {
							ss << pvar[j] << "\t";
						} else {
							ss << segsites[j] << "\t";
						}
					}
				}
				break;
			case 6:
				if (!fst.empty()) {
					for (j=0; j<fst.size(); ++j) {
						ss << fst[j] << "\t";
					}
				}
				break;
		}
	}

	std::string statstr = ss.str();
	if ((lastchar = statstr.find_last_not_of(whitespace)) != std::string::npos) {
		statstr.erase(lastchar+1);
	}
	outstream << statstr << "\n";
}

int doStats (int argc, char** argv) {
	int rv = 0;
	int verbose = 0; // controls amount of output to screen
	int writeheader = 1;
	int sfsprob = 0; // 1 = turn SFS categories into probabilities, 0 = leave SFS categories as counts
	int fold = 0; // 1 = folded sfs, 0 = unfolded
	unsigned int seqlen = 0; // sequence length for calculating per site statistics
	std::string arpname;
	std::ostream outstream(std::cout.rdbuf());
	std::ofstream outfstream;

	// container of what populations to calculate certain stats for
	// 0 = probability shared doubleton allele
	// 1 = pi (theta_tajima)
	// 2 = theta_watterson
	// 3 = Tajima's D
	// 4 = SFS, format: [SFS CATEGORY, POPULATION]
	// 5 = S, number of segregating sites within each population
	// 6 = Fst (Reynolds et al. 1983)

	const int nstats = 7; // increase this when adding more stats
	std::vector< std::vector<int*> >pops;
	pops.resize(nstats);

	// parse args
	std::vector<int> printorder;
	printorder.reserve(nstats);
	std::vector<std::string> popid;
	popid.reserve(50);

	if ((rv = statArgs(argc, argv, arpname, outfstream, outstream, pops, popid, printorder, fold, sfsprob, seqlen, writeheader, verbose))) {
		return rv;
	}

	// parse arp dna sequence input
	PopSeq seq(arpname.c_str());
	if (seq.nseq() < 1) {
		std::cerr << "No sequences read from " << arpname << "\n";
		return -1;
	}

	// make sure user-supplied pops IDs are legit

	std::string badpop = checkPopIDs(popid, seq.popn());
	if (!badpop.empty()) {
		std::cerr << "Invalid population ID " << badpop << "\n";
		return -1;
	}

	// calculate stats

	// vector to hold allele frequencies
	std::vector< std::vector<AlleleFreq> > af; // stores allele freq for every population + total population
	af.resize(seq.popn()+1);
	for (unsigned int i=0; i<af.size(); ++i) {
		af[i].resize(seq.snpn());
	}

	// matrix to hold doubleton counts
	unsigned int** double_counts = new unsigned int* [seq.popn()];
	for (int i=0; i<seq.popn(); ++i) {
		double_counts[i] = new unsigned int [seq.popn()];
		for (int j=0; j<seq.popn(); ++j) {
			double_counts[i][j] = 0;
		}
	}

	// SFS vector
	std::vector< std::vector<int> > sfs;
	sfs.resize(seq.popn());
	for (unsigned int i=0; i<sfs.size(); ++i) {
		sfs[i].resize(seq.sampsz[i]+1, 0); // sfs has categories 0,1,..,2n
	}

	// calculate minor allele frequencies
	int total_counts [5] = {0}; // [A,C,G,T,N]
	int pop_counts [5] = {0}; // [A,C,G,T,N]
	int allele_idx;

	for (unsigned int snp = 0; snp < seq.seq[0].length(); ++snp) {
		for (int j = 0; j<5; ++j) total_counts[j] = 0;
		int seqstr = 0;
		for (int popid = 0; popid < seq.popn(); ++popid) {
			for (int j = 0; j<5; ++j) pop_counts[j] = 0;
			int n = 0;
			while (n < seq.sampsz[popid]) {
				switch (seq.seq[seqstr][snp]) {
					case 'A' :
						allele_idx = 0;
						break;
					case 'C' :
						allele_idx = 1;
						break;
					case 'G' :
						allele_idx = 2;
						break;
					case 'T' :
						allele_idx = 3;
						break;
					case '?' :
						allele_idx = 4;
						break;
					default :
						std::cerr << "Unknown allele " << seq.seq[seqstr][snp] << "\n";
						return -1;
				}
				++total_counts[allele_idx];
				++pop_counts[allele_idx];
				++seqstr;
				++n;
			}
			// record population allele frequency
			af[popid][snp].allele = mafidx(pop_counts, 5);
			af[popid][snp].freq = pop_counts[af[popid][snp].allele];
			// iterate population SFS
		}
		// record total population allele frequency
		int tidx = seq.popn();
		af[tidx][snp].allele = mafidx(total_counts, 5);
		af[tidx][snp].freq = total_counts[af[tidx][snp].allele];
	}

	// calculate SFS and doubleton counts
	if (calcSFS(af, sfs, double_counts, seq.snpn())) {
		return -1;
	}

	if (verbose) {
		// print unfolded SFS to screen
		for (unsigned int i = 0; i < sfs.size(); ++i) {
			std::cerr << "population " << i << " SFS\n";
			std::cerr << sfs[i][0];
			for (unsigned int j = 1; j < sfs[i].size(); ++j) {
				std::cerr << " " << sfs[i][j];
			}
			std::cerr << "\n";
		}
	}

	// calculate number of segregating sites
	std::vector<int> segsites;
	std::vector<double> pvar;
	calcSegSites(sfs, segsites);
	if (seqlen && !pops[5].empty()) {
		pvar.reserve(segsites.size());
		for (unsigned int i=0; i<segsites.size(); ++i)
		pvar.push_back(static_cast<double>(segsites[i])/seqlen);
	}

	// calculate doubleton allele sharing probabilities
	std::vector<double> doubleprob;

	if (!pops[0].empty()) {
		if(probDoubleShare(double_counts, seq.popn(), pops[0], &doubleprob)) {
			return -1;
		}
	}

	// calculate thetas and Tajima's D

	std::vector<double> theta_w;
	std::vector<double> pi;
	std::vector<double> tajimaD;

	if (!pops[3].empty()) {
		// do Tajima's D first because this will calculate thetas
		if (calcTajimaD(sfs, segsites, seq.sampsz, pops[3], theta_w, pi, tajimaD)) {
			return -1;
		}
	}

	if (!pops[1].empty()) {
		if (pi.empty()) {
			calcPi(sfs, pi);
		}
		if (seqlen) {
			for (unsigned int i=0; i<pi.size(); ++i) {
				pi[i] /= seqlen;
			}
		}
	}

	if (!pops[2].empty()) {
		if (theta_w.empty()) {
			calcWatterson(segsites, seq.sampsz, theta_w);
		}
		if (seqlen) {
			for (unsigned int i=0; i<theta_w.size(); ++i) {
				theta_w[i] /= seqlen;
			}
		}
	}

	// calculate SFS category

	std::vector<double> sfsbin;

	if (!pops[4].empty()) {
		if(sfsBins(sfs, pops[4], fold, segsites, sfsprob, sfsbin)) {
			return -1;
		}
	}

	// calculate Fst
	std::vector<double> fst;

	if (!pops[6].empty()) {
		globalFst (af, seq.sampsz, pops[6], fst);
	}

	// output stats
	writeStats(outstream, pops, printorder, writeheader, segsites, pvar, sfsbin, doubleprob, theta_w, pi, tajimaD, fst);

	if (outfstream.is_open()) {
		outfstream.close();
	}

	// free memory
	for (int i=0; i<seq.popn(); ++i) {
			delete [] double_counts[i];
	}
	delete [] double_counts;


	for (int i=0; i<nstats; ++i) {
		for (unsigned int j=0; j < pops[i].size(); ++j) {
			delete [] pops[i][j];
		}
	}

	return rv;
}

void mainInfo () {
	int w = 12;

	std::cerr << "\nABC++utils is a program that integrates into ABC demographic inference pipelines\n";
	std::cerr << "\nABC++utils [FUNCTION] [FUNCTION ARGUMENTS]";
	std::cerr << "\n\nFUNCTIONS:\n\n"
	<< std::setw(w) << std::left << "doStats" << "Calculate summary statistics from Fastsimcoal2 DNA sequence\n\n";
}

int parseArg (int argc, char** argv) {
	int rv = 0;

	if (argc < 2 || strcmp(argv[1], "help") == 0 || strcmp(argv[1], "-help") == 0) {
		mainInfo();
		rv = 0;
	} else if (strcmp(argv[1], "doStats") == 0) {
		rv = 1;
	} else {
		std::cerr << "Unknown function " << argv[1] << "\n";
		rv = -1;
	}

	return rv;
}

int main (int argc, char** argv) {
	int rv = 0;
	int subroutine = parseArg(argc, argv);

	switch (subroutine) {
		case -1 :
			break;
		case 0 :
			break;
		case 1 :
			if ((rv = doStats(argc, argv)) >= 0) {
				rv = 0;
			}
			break;
	}

	return rv;
}

