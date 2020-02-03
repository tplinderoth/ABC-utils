/*
 * ABC++utils.h
 *
 *  Created on: Jan 18, 2020
 *      Author: tyler
 */

#ifndef ABC__UTILS_H_
#define ABC__UTILS_H_

#include <vector>
#include <cstdlib>
#include <string>

struct AlleleFreq {
	AlleleFreq ();
	int allele; // [0=A, 1=C, 2=G, 3=T, 4=?]
	unsigned int freq;
};

class PopSeq {
public:
	PopSeq (const char* filename = NULL);
	int parseArp (const char* filename);
	unsigned int nseq ();
	int popn ();
	unsigned int snpn ();
	std::vector<int> popmap; // identifies what population sequences belong to
	std::vector<std::string> seq; // holds the sequences
	std::vector<int> sampsz; // holds the haploid sample size for each population
private:
	int npops; // number of populations
	unsigned int nsnp; // number of polymorphic sites
	std::vector<char> ref; // holds reference (major) allele
	std::vector<char> alt; // holds alternate (minor) allele
};

#endif /* ABC__UTILS_H_ */
