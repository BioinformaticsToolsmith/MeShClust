/* -*- C++ -*-
 *
 * CenterCreator.h
 *
 * Author: Benjamin T James
 */
#ifndef CENTERCREATOR_H
#define CENTERCREATOR_H

#include <vector>
#include <random>
#include <map>
#include "Point.h"
#include "DivMatrix.h"
#ifdef OPENMP
#include <omp.h>
#endif

template<class T>
class CenterCreator {
public:
	CenterCreator(DivMatrix<T> &mat, T bndwidth, int iter) : dm(mat), bandwidth(bndwidth), iterations(iter) {
		gen = new std::mt19937(rd());
	};
	~CenterCreator() {
		delete gen;
	};
	std::pair<
	std::vector<Point<T>*>,
	std::map<Point<T>*, std::vector<Point<T>*>*>
	> get_centers(
		const std::vector<Point<T>*> &points, std::function<void(vector<Point<T>*> &pts)> func,
		int partitions=10);
	std::vector<Point<T>*> get_gibbs_centers(
		const std::vector<Point<T>*> &points,
		int take_at_time=500);
	std::vector<Point<T>*> get_heuristic_centers(
		const std::vector<Point<T>*> &points,
		int take_at_time=20);
private:
	T distance(Point<T>* a, Point<T>* b);
	T avg_distance(const std::vector<Point<T>*> &vec);
	T avg_distance(Point<T> *a, const std::vector<Point<T>*> &vec);
	pair<vector<vector<Point<T>*>>,vector<vector<Point<T>*>>>
		split_vector(const std::vector<Point<T>*> &pts, int parts) const;
	pair<int,bool> get_rand_index(
		const bool *const marked,
		std::uniform_int_distribution<> &idis,
		int max_to_consider=30) const;
	std::random_device rd;
	std::mt19937 *gen;
	const T bandwidth;
	const int iterations;
	DivMatrix<T> &dm;
	#ifdef OPENMP
	omp_lock_t distance_lock;
	#endif
};

#endif
