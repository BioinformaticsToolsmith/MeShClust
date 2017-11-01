#ifndef FEATURES_H
#define FEATURES_H

#include "SingleFeature.h"

template<class T>
class Feature {
public:
	Feature(std::function<double(vector<double>)> combination, std::vector<SingleFeature<T> > sf)
		: features(sf), combo(combination) {}
	double operator()(Point<T>*, Point<T>*) const;

	static double manhattan(Point<T>& p, Point<T>& q);
	static double length_difference(Point<T>& p, Point<T>& q);
	static double n2rrc(Point<T>& p, Point<T>& q, const vector<int>&, const vector<int> &);
	static double rree_k_r(Point<T>& p, Point<T>& q);
	static double intersection(Point<T>& p, Point<T>& q);
	static double jenson_shannon(Point<T>& p, Point<T>& q);
	static double pearson(Point<T>& p, Point<T>& q);
	static double simratio(Point<T>& a, Point<T>& b);
	static double squaredchord(Point<T>& a, Point<T>& b);
private:
	vector<SingleFeature<T> > features;
	std::function<double(vector<double>)> combo;
};
#endif
