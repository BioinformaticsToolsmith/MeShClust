/* -*- C++ -*-
 *
 * MeanShiftContext.h
 *
 * Author: Benjamin T James
 */
#ifndef MEANSHIFTCONTEXT_H
#define MEANSHIFTCONTEXT_H

#include <map>
#include <vector>
#include "MeanShift.h"
#include "DivMatrix.h"
/*
 * Context for the mean shift strategy
 */
template<class T>
class MeanShiftContext {
public:
/*
 * Constructor for this class.
 *
 * This takes two pointers to a vector of points (centers and data points),
 * the bandwidth parameter b,
 * the number of iterations per center n_iterations,
 * and the mean shift strategy
 */
        MeanShiftContext(std::vector<Point<T>*> *c, const std::vector<Point<T>*> *p, T band_width, int n_iterations, MeanShift<T> *ms, DivMatrix<T> &mat, std::map<Point<T>*,std::vector<Point<T>*>*> &m, int begin_indx=0, int end_indx=0) : centers(c), points(p), bandwidth(band_width), num_iterations(n_iterations), ms(ms), dm(mat), partition(m), begin_idx(begin_indx), end_idx((end_indx==0)?(p->size()-1):end_indx) {
	};
	MeanShiftContext(std::vector<Point<T>*> *c, const std::vector<Point<T>*> *p, T band_width, MeanShift<T> *ms, DivMatrix<T> &mat, std::map<Point<T>*,std::vector<Point<T>*>*> &m, int begin_indx=0, int end_indx=0) : centers(c), points(p), bandwidth(band_width), num_iterations(5), ms(ms), dm(mat), partition(m), begin_idx(begin_indx), end_idx(end_indx==0?p->size()-1:end_indx){
	};
	~MeanShiftContext() { };
	void mean_shift(int dlta=5, bool prune=true);
	Point<T>* find_nearest_neighbor(Point<double>& center, int before, int after);
        int cluster(std::map<Point<T>*,std::vector<Point<T>*>*> &m) const {
		m = partition;
		return 0;
	};
	void clear();
	void print() const;
	std::vector<Point<T>*> get_centers() const { return *centers; };
	void printOutput(string filename) const;
private:
	void compute_bandwidths(bool display=true);
	void compute_cluster(const vector<Point<T>*> &cntrs, bool do_compute=true);
	void compute_cluster_bounded(const vector<Point<T>*> &cntrs, bool do_compute=true);
	Point<T>* quadratic_cluster(Point<T>& p, const vector<Point<T>*> &cntrs) const;
	void mean_shift_iter(Point<T> *p, Point<double> &top, Point<double> &temp, int dta);
	void check_center(int index, const vector<Point<T>*> &cntr, int delta);
	void clean_up();
	std::map<Point<T>*,std::vector<Point<T>*>*> partition;
	std::map<Point<T>*, T> bandwidths;
	std::vector<Point<T>*> *centers;
	const std::vector<Point<T>*> *points;
	T bandwidth;
	int num_iterations;
	const int begin_idx, end_idx;
	T cluster_factor = 1.00;
	MeanShift<T> *ms;
	DivMatrix<T> &dm;
	const int atom_size = 5;
};

#ifdef HEADER_HACK
#ifndef MEANSHIFTCONTEXT_C
#define MEANSHIFTCONTEXT_C
#include "MeanShiftContext.cpp"
#endif
#endif

#endif
