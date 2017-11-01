/* -*- C++ -*-
 *
 * DivMatrix.cpp
 *
 * Author: Benjamin T James
 */
#include "DivMatrix.h"
#include <assert.h>
using namespace utility;

template<class T>
void DivMatrix<T>::fill()
{
//	m->fill();
}


template<class T>
tuple<T,T,T> DivMatrix<T>::lookup(int i, int j) // check
{
	auto pr = i > j ? make_pair(i, j) : make_pair(j, i);
	if (last_used_valid && last_used.first == pr) {
		return last_used.second;
	}
	auto t = (*m)[pr];
	last_used = make_pair(pr, t);
	last_used_valid = true;
	return t;
}

template<class T>
bool DivMatrix<T>::exists(Point<T> &i1, Point<T> &i2) const
{
	return true;
	//return m->exists(i1.get_id(), i2.get_id());
}

template<class T>
T DivMatrix<T>::get_div(Point<T> &i1, Point<T> &i2)
{


	return vec[i1.get_id()]->distance(*vec[i2.get_id()]);

	#ifdef DEBUG
	auto a1 = vec[i1.get_id()];
	auto a2 = vec[i2.get_id()];
	assert(a1->alt_distance(i1) == 0);
	assert(a2->alt_distance(i2) == 0);
	#endif
        T div = get<0>(lookup(i1.get_id(), i2.get_id()));
	#ifdef DEBUG
	// T d2 = i1 - i2;
	// if (div != d2) {
	// 	cout << "Not equal: " << div << " " << d2 << endl;
	// }
	#endif
	return div;
}

template<class T>
double DivMatrix<T>::get_prob_d(Point<double> &i1, Point<T> &i2)
{
	int a = i1.get_id();
	int b = i2.get_id();
	return vec[a]->prob_under(*vec[b]);
}
template<class T>
double DivMatrix<T>::get_prob(Point<T> &i1, Point<T> &i2)
{
	int a = i1.get_id();
	int b = i2.get_id();
	return vec[a]->prob_under(*vec[b]);


	auto t = lookup(a, b);
	T ret;
	if (a <= b) {
		ret = std::get<1>(t);
	} else {
		ret = std::get<2>(t);
	}
	return ret;
	//T prob = i1.prob_under(i2);
	// if (prob != get<2>(t) && prob != get<1>(t) ) {
	// 	cout << get<2>(t) << " and " << get<1>(t) << " do NOT equal " << prob << endl;
	// 	auto a1 = vec[i1.get_id()];
	// 	auto a2 = vec[i2.get_id()];
	// 	cout << "Dist: " << a1->alt_distance(i1) << " and " << a2->alt_distance(i2) << endl;
	// } else {
	// 	cout << "Works" << endl;
	// }
//	return prob;
}
template class DivMatrix<int>;
template class DivMatrix<double>;
template class DivMatrix<uint64_t>;
template class DivMatrix<uint32_t>;
template class DivMatrix<uint16_t>;
template class DivMatrix<uint8_t>;
