/* -*- C++ -*-
 *
 * MeanShiftContext.cpp
 *
 * Author: Benjamin T James
 */
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "MeanShiftContext.h"
#include "ConcreteMeanShift.h"
#include <cfenv>
#include <regex>
#include <algorithm>

/*
 * Checks to see if any other centers are
 * within the distance of the bandwidth parameter
 * If a center is, it is staged for deletion later.
 */
template<class T>
void MeanShiftContext<T>::check_center(int index, const vector<Point<T>*> &cntr, int delta)
{
	Point<T> *cur = cntr[index];
	int start = std::max(0, index - delta);
	int end = std::min((int)cntr.size()-1, index + delta);
//	for (Point<T> *p : cntr) {
	for (int i = start; i <= end; i++) {
//		check_center(index, i, cntr);
		Point<T> *p = cntr[i];
		if (p == cur || p->is_to_delete()) {
			continue;
		}
//		cout << "Magnitude: " << p->magnitude() << endl;
		T dist = dm.get_div(*cur, *p);
//		T dist = cur->distance(*p);
		//	cout << "Trying to merge: " << dist << endl;

		if (dist < cluster_factor * bandwidth) {
//    			cout << "Distance to center: " << dist << endl;
//			cout << "Bandwidth: " << bandwidth << endl;
			try {
//				auto similarity = cur->alt_distance(*p);
				double prob = cur->prob_under(*p);
//				cout << "Similarity: " << similarity << endl;
//				cout << "Prob: " << prob << endl;
				cur->set_to_delete(true);
			} catch (std::exception& e) {
				cout << "Center 1" << endl;
				cur->display();
				cout << endl << "Center 2" << endl;
				p->display();
				cout << endl;
			}
			break;
		}
	}
}

template<class T>
void MeanShiftContext<T>::compute_bandwidths(bool display)
{
	if (display) {
		cout << "Computing bandwidths... ";
		cout.flush();
	}
	bandwidths.clear();
	bool redo = false;
	for (auto const& kv : partition) {
		bandwidths[kv.first] = ms->variance(*kv.first, *kv.second);
		// T best = 0;
		// bool is_set = false;
		// T total = 0;
		// for (auto p : *kv.second) {
		// 	T tmp = kv.first->distance(*p);
		// 	total += tmp;
		// 	if (!is_set || tmp > best) {
		// 		is_set = true;
		// 		best = tmp;
		// 	}
		// }
		// bandwidths[kv.first] = best;
		// auto avg = total / kv.second->size();
		// cout << "Best: " << best << "\tAverage: " << avg << endl;
	}
	if (display) {
		cout << "Done" << endl;
	}
}

template<class T>
int get_template(Point<T>& p) {
	string header = p.get_header();
	int result = -1;
	try {
		std::regex re("template_(\\d+)");
		std::smatch match;
		if (std::regex_search(header, match, re) && match.size() > 1) {
			string res = match.str(1);
			result = std::stoi(res, NULL, 10);
		}
	} catch (std::regex_error& e) {
		std::cerr << "regex error: " << e.what() << std::endl;
	}
	return result;

}
template<class T>
Point<T>* MeanShiftContext<T>::quadratic_cluster(Point<T>& p, const vector<Point<T>*> &cntrs) const
{
	T best_dist;
	Point<T> *c = NULL;
	for (int j = 0; j < cntrs.size(); j++) {
		Point<T> *cur = cntrs[j];
		T dist = dm.get_div(p, *cur);
		if (c == NULL || dist < best_dist) {
			best_dist = dist;
			c = cur;
		}
	}
	return c;
}

template<class T>
void MeanShiftContext<T>::compute_cluster_bounded(const vector<Point<T>*> &cntrs, bool display)
{
	cout << "Centers: " << cntrs.size() << endl;
	std::map<Point<T>*,std::vector<Point<T>*>*> m;
	for (Point<T> *c : cntrs) {
		m[c] = new std::vector<Point<T>*>();

	}

	int before_idx = 0;
	int after_idx = 0;
	int window = 3;
	for (int i = begin_idx; i <= end_idx; i++) {
		int last = before_idx;
		// makes `last` the largest index before index I.
		while (before_idx < cntrs.size() && cntrs[before_idx]->get_id() <= i) {
			last = before_idx;
			before_idx++;
		}
		before_idx = last;
		const int n = (int)cntrs.size() - 1;
		after_idx = n;//min(last+10, n); // start high, move down
		//after_idx should be >= I but closest to I
		while (after_idx >= 0 && cntrs[after_idx]->get_id() >= i) {
			last = after_idx;
			after_idx--;
		}
		after_idx = last;


		Point<T> *point = points->at(i);
//		cout << "Index: " << point->get_id() << " Before: " << before_idx << " [" << cntrs[before_idx]->get_id() << "] After: " << after_idx << " [" << cntrs[after_idx]->get_id() << "] ";

		assert(before_idx <= after_idx);
		assert(after_idx < cntrs.size());
		int pt_template = get_template(*point);
		Point<T> *cnter = NULL;
		if (before_idx == after_idx) {
			cnter = cntrs[before_idx];
		} else {
			vector<Point<T>*> to_try;
			for (int j = max(0,before_idx-window); j <= min(n,after_idx+window); j++) {
				to_try.push_back(cntrs[j]);
			}
		        cnter = quadratic_cluster(*point, to_try);
		}
		double relax_pct = 0.9;
		if(display){
			double relaxed_bandwidth = ((1-relax_pct) * point->size()) + (relax_pct * bandwidth);
			T dist = cnter->distance(*point);
			if (dist > relaxed_bandwidth) {
//			cout << "Point " << point->get_id() << " matches " << endl;
				point->set_data_str("beyond");
				cnter = quadratic_cluster(*point, *centers);
			}
		}
//		cout << "Relaxed bandwidth: " << relaxed_bandwidth <<  " from " << bandwidth << "\t";
		m[cnter]->push_back(point);
		// if (get_template(*cnter) == pt_template) {
		// 	cout << "CORRECT" << endl;
		// } else {
		// 	cout << "INCORRECT" << endl;
		// }

	}
	for (auto const& kv : partition) {
		delete kv.second;
	}
	partition = m;
	if(display){
		print();
	}
}
/*
 * Clusters the points according to the closest centers
 * using the distance to the closest center
 */
template<class T>
void MeanShiftContext<T>::compute_cluster(const vector<Point<T>*> &cntrs, bool do_compute)
{
	return compute_cluster_bounded(cntrs, do_compute);
	bool display = true;
	if (display) {
		cout << "Getting clusters... ";
		cout.flush();
	}
	if (do_compute) {
		std::map<Point<T>*,std::vector<Point<T>*>*> m;
		for (Point<T> *c : cntrs) {
			m[c] = new std::vector<Point<T>*>();
		}
		if (cntrs.size() == 0) {
			cerr << "Empty partition" << endl;
			partition = m;
			return;
		}
		for (int i = begin_idx; i <= end_idx; i++) {
			Point<T> *point = points->at(i);
			cout << "Point " << point->get_id() << " of " << point->get_header();// << endl;
			Point<T> *c = cntrs[0];
			if (c == NULL) {
				perror("Null center");
				return;
			}
			T c_prob = dm.get_prob(*point, *c);

			int best_index = 0;
			for (int j = 0; j < cntrs.size(); j++) {
				Point<T> *cur = cntrs[j];
				T prob = dm.get_prob(*point, *cur);
				if (prob > c_prob) {
					c_prob = prob;
					c = cur;
					best_index = j;
				}
			}
			cout << " Best center: " << best_index << " [" << c->get_id() << "] from " << c->get_header() << endl;
			m[c]->push_back(point);
		}
		for (auto const& kv : partition) {
			delete kv.second;
		}
		partition = m;
	}
	if (display) {
		cout << "Found" << endl;
	}

	compute_bandwidths(display);
}

/*
 * frees the memory allocated by the driver class for the centers and points
 */
template<class T>
void MeanShiftContext<T>::clear()
{
	// for (Point<T> *p : *centers) {
	// 	delete p;
	// }
	// for (Point<T> *p : *points) {
	// 	delete p;
	// }
	// centers->clear();
}


template<class T>
void print_output(const string& output, const map<Point<T>*,vector<Point<T>*>*> & partition)
{
	cout << "Printing output" << endl;
	std::ofstream ofs;
	ofs.open(output, std::ofstream::out);
	int counter = 0;
	for (auto const& kv : partition) {
		if (kv.second->size() == 0) {
			continue;
		}
		ofs << ">Cluster " << counter << endl;
		int pt = 0;
		for (auto p : *kv.second) {
			string s = p->get_header();
			ofs << pt << "\t" << p->get_length() << "nt, " << s << "... " << endl;
//			string fa = am.get(p->get_id());
//			ofs << writefa(fa) << endl;
			pt++;
		}
		counter++;
	}
	ofs.close();
}
template<class T>
void MeanShiftContext<T>::printOutput(string filename) const
{
	print_output(filename, partition);
}

template<class T>
void MeanShiftContext<T>::print() const
{
	map<int,int> table;
	vector<int> ids;
	for (int i = 0; i < centers->size(); i++) {
		ids.push_back(centers->at(i)->get_id());
	}
	for (auto const& kv : partition) {
//		kv.first is center
//              *kv.second is points
		for (auto p : *kv.second) {
			table[p->get_id()] = kv.first->get_id();
		}
	}
	for (auto const& kv : table) {
		string hdr = points->at(kv.first)->get_header();
		cout << kv.first << "\t" << hdr << "\t" << kv.second << "\t";
//		if (kv.first != 0) {
			T dist = points->at(kv.first)->distance(*points->at(kv.second));
			cout << dist << "\t";
//		}
		if (kv.first == kv.second) {
			cout << "* ";
			for (int i = 0; i < ids.size(); i++) {
				if (kv.first == ids[i]) {
					if (i > 0) {
						cout << ids[i - 1] << "=";
						cout << points->at(ids[i-1])->distance(*points->at(kv.first));
					}
					if (i + 1 < ids.size()) {
						cout << " " << ids[i + 1] << "=";
						cout << points->at(ids[i+1])->distance(*points->at(kv.first));
					}
				}
			}

		}
		cout << points->at(kv.first)->get_data_str();
		cout << endl;
	}
}


/*
 * Persistent memory store of dataPts
 */
template<class T>
Point<T>* MeanShiftContext<T>::find_nearest_neighbor(Point<double> &center, int before, int after)
{
	// cout << "Correcting centers...";// Quadratic bottleneck
	// cout.flush();
	// for (auto *c : *centers) {
	// 	auto clstr = partition[c]; // center, although moves, still points to cluster;
	Point<T> *best = points->at(before); // pointer still is center
	T min = best->distance_d(center);
	for (int i = before + 1; i <= after; i++) {
		Point<T>* p = points->at(i);
		T dist = p->distance_d(center);
		if (dist < min) {
			best = p;
			min = dist;
		}
	}
	return best;
	// }
	// cout << " Done" << endl;

	// cout << "Correcting centers...";// Quadratic bottleneck
	// cout.flush();
	// for (auto *c : *centers) {
	// 	auto clstr = partition[c]; // center, although moves, still points to cluster;
	// 	Point<T> *best = points->at(0); // pointer still is center
	// 	T min = best->alt_distance(*c);
	// 	for (auto *p : *points) {
	// 		T dist = p->alt_distance(*c);
	// 		if (dist < min) {
	// 			best = p;
	// 			min = dist;
	// 		}
	// 	}
	// 	c->set(*best);
	// }
	// cout << " Done" << endl;
}

/*
 * Mean shift algorithm
 *
 * Iterates through each center (independently)
 * and uses gradient ascent to reach maxima.
 *
 * N.B.: Points 'top' and 'temp' are created on the heap
 * here to not stall performace inside the inner loops
 * by allocating and freeing.
 *
 * Parallelization is very possible for every center
 * if each iteration of the for loop is put on a separate
 * thread. This requires that a mutex be put on the deletion
 * flag of the point so the 'deleted' centers aren't included.
 */
template<class T>
void MeanShiftContext<T>::mean_shift(int delta, bool prune)
{
	cout << "Beginning Mean Shift" << endl;
	//cout << "Range: [" << begin_idx << "," << end_idx << "]" << endl;
	// // First thing we do is merge
	// for (auto *c : *centers) {
	// 	check_center(c, *centers);
	// }
	// for (int i = begin_idx; i <= end_idx; i++) {
	// 	cout << "BEGIN  " << points->at(i)->get_id() << "\t" << points->at(i)->get_header() << endl;
	// }
	Point<double> *top = (*centers)[0]->create_double();
	Point<double> *temp = top->create();
	for (int i = 0; i < num_iterations; i++) {
		cout << "iter " << i + 1 << "/" << num_iterations << endl;
		// compute_cluster(*centers, i != 0);
		// string filename = "output_" + std::to_string(i+1) + ".clstr";
		// print_output(filename, partition);
		cout << "Centers: " << centers->size() << endl;
		for (auto *c : *centers) {
			mean_shift_iter(c, *top, *temp, delta);
		}
//		find_nearest_neighbor();
		for (int a = 0; a < centers->size(); a++) {
			check_center(a, *centers, delta + 1);
		}
		clean_up();
		std::sort(centers->begin(), centers->end(), [](const Point<T>* a, const Point<T> *b) { return a->get_id() < b->get_id(); });
		// if (i == num_iterations - 1) {
		// 	for (auto c : *centers) {
		// 		cout << "CENTER: " << c->get_header() << endl;
		// 	}
		// }
	}
	//compute_cluster(*centers, true);
	delete top;
	delete temp;
	clean_up();


	bool found_atoms = true;
//	std::map<Point<T>*,std::vector<Point<T>*>*> partition;
//	cluster(partition);
	compute_cluster(*centers, false);

	while (prune && found_atoms && centers->size() > 0) {
		cout << "Finding atoms..."<< endl;
		cout.flush();
		found_atoms = false;
		for (auto const& kv : partition) {
			if (kv.second->size() <= atom_size) {
				found_atoms = true;
				int i = 0;
				for (auto v : *centers) {
					if (v == kv.first) {
						v->set_to_delete(true);
						break;
					}
				}
			}
		}
		clean_up();
		compute_cluster(*centers, false);
//		cluster(partition);
	}
	compute_cluster(*centers, true);

	// for (Point<T> *center : *centers) {
	// 	auto pts_to_consider = partition[center];
	// 	vector<string> indices;
	// 	for (auto p : *pts_to_consider) {
	// 		indices.push_back(p->get_header());
	// 	}
	// 	for (auto idx : indices) {
	// 		cout << idx << " ";
	// 	}
	// 	cout << endl << endl;
	// }
}

/*
 * Inner loop of mean_shift() which calculates
 * the numerator and denominator sums and sets
 * the new center to be the quotient
 */
template<class T>
void MeanShiftContext<T>::mean_shift_iter(Point<T> *center, Point<double> &top, Point<double> &temp, int delta)
{
	double bottom = 0;
	top.zero();
	std::feclearexcept(FE_OVERFLOW);
	std::feclearexcept(FE_UNDERFLOW);

	int before_idx = center->get_id() - delta;
	int after_idx = center->get_id() + delta;
	if (before_idx < 0) {
		before_idx = begin_idx;
	}
	if (after_idx >= points->size()) {
		after_idx = end_idx;
	}
//	cout << "Range: [" << before_idx << "," << after_idx << "]" << endl;
	for (int i = before_idx; i <= after_idx; i++) {
		points->at(i)->set_arg_to_this_d(temp);
//		temp.set_template<typeid(T).name()>(*points->at(i));
		// When used with bounds it has the same effect as the truncated flat kernel
		double multiplier = ms->kernel(temp, *center, bandwidth);// * ms->weight(*temp); // e^-x^2 is only != 0 when x = 0
		//temp->subOne();
		temp *= multiplier;
//		cout << "Point in ms:" << endl;
//		temp->display();
		top += temp;
		bottom += multiplier;
//		cout << "Multiplier: " << multiplier << " Bottom: " << bottom << endl;
		if ((bool)std::fetestexcept(FE_UNDERFLOW)) {
			cout << "Underflow!" << endl;
			throw std::exception();
		}
	}
//	cout << "Denom: " << bottom << endl;
	if (bottom) {
		top /= bottom;
	}
//	cout << "Real magnitude: " << top->getRealMagnitude() << endl;
//	T dist = *center - top;
//	cout << "Center distance: " << dist << endl;

	//center->addOne();

	center->set(*find_nearest_neighbor(top, before_idx, after_idx));

}

/*
 * Removes and frees the centers previously lazy-deleted
 */
template<class T>
void MeanShiftContext<T>::clean_up()
{
	vector<Point<T>*> *to_keep = new vector<Point<T>*>();
	for (auto &c : *centers) {
		if (!c->is_to_delete()) {
			to_keep->push_back(c);
		} else {
			delete c;
			//c = NULL;
		}
	}
	centers->clear();
//	delete centers;
	//TODO: memory leak
	centers = to_keep;
}

#ifndef HEADER_HACK
template class MeanShiftContext<double>;
template class MeanShiftContext<int>;
template class MeanShiftContext<uint64_t>;
template class MeanShiftContext<uint32_t>;
template class MeanShiftContext<uint16_t>;
template class MeanShiftContext<uint8_t>;
#endif
