/* -*- C++ -*-
 *
 * needleman_wunsch.h
 *
 * Author: Benjamin T James
 */

#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include <string>

class needleman_wunsch {
public:
	needleman_wunsch(const std::string& s1, const std::string& s2, int match_,  int mismatch_, int sigma_, int epsilon_);
	~needleman_wunsch() {
		delete[] score;
		delete[] direction;
		delete[] horiz_gap_len;
		delete[] vert_gap_len;
	}
	double identity(std::pair<std::string, std::string> p) const;
	std::pair<std::string, std::string>
	align();
private:
	int gap(int gap_len) const;
	int match_score(char a, char b) const;
	inline int at(int a, int b) const { return a * l2 + b; };
	void fill(int,int);
	std::pair<std::string, std::string> backtrack();
	int match, mismatch;
	int sigma, epsilon;
	std::string s1, s2;
	int l1, l2;

	int *score;
	uint8_t *direction;
	int *horiz_gap_len;
	int *vert_gap_len;
};


#endif
