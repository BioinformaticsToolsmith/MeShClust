/* -*- C++ -*-
 *
 * ConcreteMeanShift.cpp
 *
 * Author: Benjamin T James
 */
#ifndef HEADER_HACK
#include "ConcreteMeanShift.h"
#endif

#include <math.h>
#include <iostream>

/*
 * Gaussian kernel
 */
template<class T>
double ConcreteMeanShift<T>::kernel(Point<double>& x, Point<T>& y, T bandwidth)
{
// This one works
// return dm.get_prob_d(x, y);


// if (rate < 0.69) {
	// 	return 0;
	// }


//	const auto dist = dm.get_div(x, y);
	//double dist = p.distance_d(c);
	//auto r = -1 * dist * dist; /// (bandwidth * bandwidth);
	//return exp(r);
	return dm.get_prob_d(x, y);
}

template<class T>
T ConcreteMeanShift<T>::weight(Point<T>& x)
{
	return 1; // we don't want any weight factor
}

template<class T>
T ConcreteMeanShift<T>::variance(Point<T> &center, const std::vector<Point<T> *> &pts) const
{

 	int count = 0;
 	T mean = 0, sim_mean = 0, prob_mean = 0;
	T var = 0, sim_var = 0, prob_var = 0;
	for (int i = 0; i < pts.size(); i++) {
		T dist = dm.get_div(*pts[i], center);
		T prob_dist = dm.get_prob(*pts[i], center);
		mean += dist;
		prob_mean += prob_dist;
	}
	mean /= pts.size();
	sim_mean /= pts.size();
	prob_mean /= pts.size();
	for (int i = 0; i < pts.size(); i++) {
		T dist = dm.get_div(*pts[i], center);
		T prob_dist = dm.get_prob(*pts[i], center);

		var += (dist - mean) * (dist - mean);
//		sim_var += (sim_dist - sim_mean) * (sim_dist - sim_mean);
		prob_var += (prob_dist - prob_mean) * (prob_dist - prob_mean);
	}
	var /= (pts.size() * pts.size());
	sim_var /= pts.size();
	prob_var /= pts.size();
	// cout << "Mean: " << mean << "\t" << sim_mean << "\t" << prob_mean << endl;
	// cout << "SD: " << var << "\t" << sim_var << "\t" << prob_var << endl;
	// cout << "Size: " << pts.size() << endl;
	return var;
// //	cout << "Cluster size: " << pts.size() << endl;
// 	T avg = sum / count;
// //	cout << "Mean: " << avg << endl;
// 	sum = 0;
// 	for (int i = 0; i < pts.size(); i++) {
// 		for (int j = i + 1; j < pts.size(); j++) {
// 			auto dist = (*pts[i] - *pts[j]) - avg;
// 			sum += dist * dist;
// 		}
// 	}
// //	auto var = sum / count;
// 	// cout << "Standard deviation: " << sqrt(var) << endl;
// 	// return var;

//  // This produced 100%
// 	// T var = 0;
// 	// auto n = pts.size();
// 	// for (auto const p : pts) {
// 	// 	T dist = (*p - m);
// 	// 	T prob = p->prob_under(m);
// 	// 	var += dist * dist * prob;
// 	// }
//         // auto result = var / (n * n);
// //	cout << "Variance: " << result << endl;
// //	return result;

}


template<class T>
T ConcreteMeanShift<T>::get_rate(Point<T>& x, Point<T>& y) const
{

//	T sim = x.alt_distance(y);
	return  x.prob_under(y);
//	return (sim < 0.35) ? 0 : (sim * x.prob_under(y));
//	return sim * x.prob_under(y);
//	T rate = 1.0;
//	if (use_alt) {
//		rate *= x.alt_distance(y);
		//	cout << "Using truncated kernel with similarity " << rate;
		// if (rate < 0.25) {
		//  	rate = 0;
		//  	return 0;
		// }

//		auto prob = x.prob_under(y);
		//cout << "Divergence: " << (x - y) << "\tIntersection: " << rate << "\tProbability: " << prob << endl;
//		rate *= prob; // rate == similarity
		//cout << " Probability: " << prob;
		//cout << "Intersection: " << rate << endl;
		// if (rate == 0.0) {
		// 	return 0;
		// }
//	}
//	return rate;
//		return 1.0;
}
#ifndef HEADER_HACK
template class ConcreteMeanShift<double>;
template class ConcreteMeanShift<int>;
template class ConcreteMeanShift<uint64_t>;
template class ConcreteMeanShift<uint32_t>;
template class ConcreteMeanShift<uint16_t>;
template class ConcreteMeanShift<uint8_t>;
#endif
