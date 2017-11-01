/* -*- C++ -*-
 *
 * Runner.h
 *
 * Author: Benjamin T James
 */
#ifndef RUNNER_H
#define RUNNER_H

#include <iostream>
#include <map>
#include "Point.h"
using namespace std;

class Runner {
public:
	Runner(int argc, char** argv);
	~Runner() {};
	int run() const;
private:
	template<class T> int do_run() const;
	template<class T> void print_output(const map<Point<T>*, vector<Point<T>*>*> &m) const;
	int k = -1;
        int bandwidth;
	double similarity = 90;
	long longest_length;
	int iterations = 15;
	int delta = 5;
	std::vector<std::string> files;
	string output = "output.clstr";
	void get_opts(int argc, char** argv);
	pair<int,uint64_t> find_k() const;
	void compute_bandwidth();
};
#endif
