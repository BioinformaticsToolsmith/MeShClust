/* -*- C++ -*-
 *
 * ConcreteMeanShift.h
 *
 * Author: Benjamin T James
 */
#ifndef CONCRETEMEANSHIFT_H
#define CONCRETEMEANSHIFT_H

#include "MeanShift.h"
#include "DivMatrix.h"
#include <vector>
/*
 * Concrete strategy for the mean shift algorithm
 */
template<class T>
class ConcreteMeanShift : public MeanShift<T> {
public:
	ConcreteMeanShift(DivMatrix<T> &mat) : dm(mat) {};
	~ConcreteMeanShift() { }
	double kernel(Point<double>& x, Point<T>& y, T bandwidth);
	T weight(Point<T>& x);
	T variance(Point<T> &m, const std::vector<Point<T> *> &pts) const;
private:
	T get_rate(Point<T> &x, Point<T> &y) const;
	DivMatrix<T> &dm;
};

#ifdef HEADER_HACK
#ifndef CONCRETEMEANSHIFT_C
#define CONCRETEMEANSHIFT_C
#include "ConcreteMeanShift.cpp"
#endif
#endif

#endif
