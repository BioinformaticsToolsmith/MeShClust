/* -*- C++ -*-
 *
 * MeanShift.h
 *
 * Author: Benjamin T James
 */
#ifndef MEANSHIFT_H
#define MEANSHIFT_H

#include "Point.h"

/*
 * Strategy interface for the mean shift algorithm
 */
template<class T>
class MeanShift {
public:
	virtual ~MeanShift() {}
//	virtual T kernel(T x) = 0;
	virtual double kernel(Point<double>& x, Point<T>& y, T bandwidth) = 0;
	virtual T weight(Point<T>& x) = 0;
	virtual T variance(Point<T>& p, const std::vector<Point<T>*> &pts) const = 0;
};
#endif
