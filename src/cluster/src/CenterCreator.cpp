/* -*- C++ -*-
 *
 * CenterCreator.cpp
 *
 * Author: Benjamin T James
 */
#include "CenterCreator.h"

#include "MeanShiftContext.h"
#include "ConcreteMeanShift.h"
#ifdef OPENMP
#include <omp.h>
#endif

template<class T>
std::vector<Point<T>*> CenterCreator<T>::get_gibbs_centers(const std::vector<Point<T>*> &points, int take_at_time)
{
	auto n = points.size();
	std::uniform_real_distribution<> rdis(0, 1);
	std::uniform_int_distribution<> idis(0, n-1);

	vector<Point<T>*> total_centers;
#ifdef OPENMP
	omp_lock_t writelock;
        omp_init_lock(&writelock);
	omp_init_lock(&distance_lock);
	#pragma omp parallel for
#endif
	for (int q = 0; q < 8; q++) {
		int num_take_at_time = q >= 4 ? take_at_time / 2 : take_at_time;
		std::vector<Point<T>*> centers;
		bool *marked = new bool[n]();
		int irand = get_rand_index(marked, idis).first;
		centers.push_back(points[irand]->clone());
		T old_avg_dist = 0;
		bool do_continue = true;
		while (do_continue) {
			int to_add = -1;
			std::map<int, T> scores;
			double total = 0;
			for (int i = 0; i < num_take_at_time; i++) {
				int r = get_rand_index(marked, idis).first;
				T avg = avg_distance(points[r], centers);
				scores[r] = avg;
				total += avg;
			}
			double rnum = rdis(*gen);
			int last = -1;
			double prob = 0;
			for (auto const& kv : scores) {
				prob += kv.second / total;
				to_add = kv.first;
				if (prob >= rnum) {
					break;
				}
			}
			marked[to_add] = true;
			centers.push_back(points[to_add]->clone());
			//cout << points[to_add]->get_header() << endl;

			T avg_dist = avg_distance(centers);
			do_continue = (avg_dist > old_avg_dist);
			old_avg_dist = avg_dist;
		}
		auto p = centers[centers.size()-1];
		centers.pop_back();
		delete p;
		delete[] marked;

#ifdef OPENMP
		omp_set_lock(&writelock);
#endif
		cout << "Adding " << centers.size() << " centers" << endl;
		for (auto c : centers) {
			total_centers.push_back(c);
		}
#ifdef OPENMP
		omp_unset_lock(&writelock);
#endif

	}
        #ifdef OPENMP
	omp_destroy_lock(&writelock);
	omp_destroy_lock(&distance_lock);
	#endif
	//cout << "Found " << centers.size() << " centers" << endl;
	return total_centers;
}
template<class T>
std::vector<Point<T>*> CenterCreator<T>::get_heuristic_centers(const std::vector<Point<T> *> &points, int take_at_time)
{
	const auto n = points.size();
	vector<Point<T>*> centers;
	std::uniform_int_distribution<> idis(0, n-1);
	bool *marked = new bool[n]();

	auto rand_set = [&]() {
		vector<int> v;
		for (int j = 0; j < take_at_time; j++) {
			int idx = get_rand_index(marked, idis).first;
			v.push_back(idx);
		}
		return v;
	};
	auto reduce = [&](const vector<int> &vec) -> vector<int> {
		vector<int> ret;
		bool *flags = new bool[vec.size()]();
		for (auto p : vec) {
			for (auto b : vec) {
				if (b == p || flags[b]) {
					continue;
				}
				T dist = distance(points[p], points[b]);
				if (dist < bandwidth) {
					flags[p] = true;
					break;
				}
			}
			if (!flags[p]) {
				ret.push_back(p);
			}
		}
		delete[] flags;
		return ret;
	};
	#ifdef OPENMP
	omp_lock_t writelock;
        omp_init_lock(&writelock);
	omp_init_lock(&distance_lock);
	#pragma omp parallel for
	#endif
	for (int i = 0; i < points.size() / take_at_time; i++) {
		auto set = reduce(rand_set());
		#ifdef OPENMP
		omp_set_lock(&writelock);
		#endif
		cout << "Pushing " << set.size() << " centers [" << i * take_at_time / points.size() * 100.0 << "%]" << endl;
		for (auto s : set) {
			marked[s] = true;
			centers.push_back(points[s]->clone());
		}
		#ifdef OPENMP
		omp_unset_lock(&writelock);
		#endif

	}
	#ifdef OPENMP
	omp_destroy_lock(&writelock);
	omp_destroy_lock(&distance_lock);
	#endif

	delete[] marked;
	return centers;
}

template<class T>
pair<vector<vector<Point<T>*>>,vector<vector<Point<T>*>>> CenterCreator<T>::split_vector(const std::vector<Point<T>*> &pts, int vec_size) const
{
	vector<vector<Point<T>*>> vec;
	vector<vector<Point<T>*>> pt_vec;
	int i = 0;
	vector<Point<T>*> temp;
	vector<Point<T>*> pt_temp;
	for (auto p : pts) {
		temp.push_back(p->clone());
		pt_temp.push_back(p);
		if (i % vec_size == vec_size - 1) {
			vec.push_back(temp);
			pt_vec.push_back(pt_temp);
			temp.clear();
			pt_temp.clear();
		}
		i++;
	}
	if (temp.size() != 0) {
		vec.push_back(temp);
		pt_vec.push_back(pt_temp);
	}
	cout << "Split into " << vec.size() << " sub-MS" << endl;
	for (auto p : pt_vec) {
		cout << "\t:" << p.size() << endl;
	}
	return make_pair(vec, pt_vec);
}
template<class T>
std::pair<
	std::vector<Point<T>*>,
	std::map<Point<T>*, std::vector<Point<T>*>*>
	> CenterCreator<T>::get_centers(const std::vector<Point<T> *> &points, std::function<void(vector<Point<T>*> &pts)> func, const int partitions)
{
	cout << "Size: " << points.size() << endl;
	for (auto p : points) {
		cout << "PRE " << p->get_header() << endl;
	}
	#ifdef OPENMP
	omp_lock_t writelock;
        omp_init_lock(&writelock);
	#endif

	auto n = points.size();
	std::vector<Point<T>*> centers;
        auto vec = split_vector(points, partitions);
	std::map<Point<T>*, std::vector<Point<T>*>*> partition;
	#ifdef OPENMP
        #pragma omp parallel for shared(dm)
	#endif
	for (int i = 0; i < vec.first.size(); i++) {
	           std::map<Point<T>*, std::vector<Point<T>*>*> p;

		   int begin = vec.second[i][0]->get_id();
		   int end = vec.second[i].back()->get_id();
		   MeanShiftContext<T> context(&vec.first[i], &points, bandwidth, iterations, new ConcreteMeanShift<T>(dm), dm, p, begin, end); //TODO race condition because of DivMatrix
		   context.mean_shift(3, false);
		auto ms_centers = context.get_centers();

		context.cluster(p);

                #ifdef OPENMP
		omp_set_lock(&writelock);
		#endif
		//centers.insert(std::end(centers), std::begin(ms_centers), std::end(ms_centers));
		// for (auto c : ms_centers) {
	        //        centers.push_back(c);
                // }
		for (auto const& kv : p) {
	                partition[kv.first] = kv.second;
			centers.push_back(kv.first);
			cout << "Center: " << kv.first->get_header() << " ";
			for (auto pt : *kv.second) {
	cout << pt->get_header() << " ";
}
			cout << endl;
                }

		context.printOutput("file_" + to_string(i) + ".clstr");
		#ifdef OPENMP
		omp_unset_lock(&writelock);
		#endif
	}
	#ifdef OPENMP
	omp_destroy_lock(&writelock);
	#endif
	return make_pair(centers, partition);
}

template<class T>
pair<int,bool> CenterCreator<T>::get_rand_index(const bool *const marked, std::uniform_int_distribution<> &idis, int max_to_consider) const
{
        bool bad = true;
        int r;
        for (int i = 0; i < max_to_consider && bad; i++) {
                bad = false;
                r = idis(*gen);
		if (marked[r]) {
			bad = true;
		}
        }
        return make_pair(r, bad);
}


template<class T>
T CenterCreator<T>::avg_distance(const vector<Point<T>*> &vec)
{
	T sum = 0;
	for (auto p : vec) {
		sum += avg_distance(p, vec);
	}
	return sum / vec.size();
}
template<class T>
T CenterCreator<T>::avg_distance(Point<T> *a, const std::vector<Point<T>*> &vec)
{
	T sum = 0;
	for (auto p : vec) {
		sum += distance(a, p);
	}
	return sum / vec.size();
}


template<class T>
T CenterCreator<T>::distance(Point<T>* a, Point<T>* b)
{
	auto p = a > b ? make_pair(a, b) : make_pair(b, a);
	return dm.get_div(*a, *b);

// 	if (!dm.exists(*a, *b)) {
// #ifdef OPENMP
// 		omp_set_lock(&distance_lock);
// #endif

// #ifdef OPENMP
// 		omp_unset_lock(&distance_lock);
// #endif

// 	}

}

template class CenterCreator<int>;
template class CenterCreator<double>;
template class CenterCreator<uint64_t>;
template class CenterCreator<uint32_t>;
template class CenterCreator<uint16_t>;
template class CenterCreator<uint8_t>;
