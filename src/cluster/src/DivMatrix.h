/* -*- C++ -*-
 *
 * DivMatrix.h
 *
 * Author: Benjamin T James
 */

#ifndef DIVMATRIX_H
#define DIVMATRIX_H
#include <vector>
#include <iostream>
#include <functional>
#include <tuple>
#include "Point.h"
#include "Mat.h"
using namespace std;
template<class T>
class DivMatrix {
public:
	DivMatrix(vector<Point<T>*> &v) : vec(v) {
/*		m = new Mat<tuple<T,T,T>>([&](int i, int j) {
				return make_tuple(v[i]->distance(*v[j]),
						  v[i]->prob_under(*v[j]),
						  v[j]->prob_under(*v[i]));
			},
			v.size());
*/
	}
	~DivMatrix() { // delete m;
	};
	void fill();
	void reset(vector<Point<T>*> &v) {
		vec = v;
		// delete m;
		// m = new Mat<tuple<T,T,T>>([&](int i, int j) {
		// 		return make_tuple(vec[i]->distance(*vec[j]),
		// 				  vec[i]->prob_under(*vec[j]),
		// 				  vec[j]->prob_under(*vec[i]));
		// 	},
		// 	v.size());
		last_used_valid = false;
	}
	bool exists(Point<T> &i1, Point<T> &i2) const;
	T get_div(Point<T> &i1, Point<T> &i2);
	double get_prob(Point<T> &i1, Point<T> &i2);
	double get_prob_d(Point<double> &i1, Point<T> &i2);

private:
	Mat<tuple<T,T,T>> *m;

	bool last_used_valid = false;
	pair<pair<int, int>, tuple<T,T,T>> last_used;

	tuple<T,T,T> lookup(int i, int j);

	vector<Point<T>*> &vec;
};
#endif
