/*
 * ChromListMaker.cpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Hani Zakaira Girgis
 */

#include "ChromListMaker.h"

namespace nonltr {

ChromListMaker::ChromListMaker(string seqFileIn) {
	seqFile = seqFileIn;
	chromList = new vector<Chromosome *>();
}

ChromListMaker::~ChromListMaker() {
	Util::deleteInVector(chromList);
	delete chromList;
}

const vector<Chromosome *> * ChromListMaker::makeChromList() {
	ifstream in(seqFile.c_str());
	bool isFirst = true;
	Chromosome * chrom;

	while (in.good()) {
		string line;
		getline(in, line);
		if (line[0] == '>') {
			if (!isFirst) {
				chrom->finalize();
				chromList->push_back(chrom);
			} else {
				isFirst = false;
			}

			chrom = new Chromosome();
			chrom->setHeader(line);
		} else {
			chrom->appendToSequence(line);
		}
	}
	chrom->finalize();
	chromList->push_back(chrom);
	in.close();

	return chromList;
}

const vector<Chromosome *> * ChromListMaker::makeChromOneDigitList() {
	ifstream in(seqFile.c_str());
	bool isFirst = true;
	ChromosomeOneDigit * chrom;

	while (in.good()) {
		string line;
		getline(in, line);
		if (line[0] == '>') {
			if (!isFirst) {
				chrom->finalize();
				chromList->push_back(chrom);
			} else {
				isFirst = false;
			}

			chrom = new ChromosomeOneDigit();
			chrom->setHeader(line);
		} else {
			chrom->appendToSequence(line);
		}
	}

	chrom->finalize();
	chromList->push_back(chrom);
	in.close();

	return chromList;
}

}
/* namespace nonltr */
