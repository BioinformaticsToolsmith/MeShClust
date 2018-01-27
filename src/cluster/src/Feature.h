#ifndef FEATURES_H
#define FEATURES_H

#include "SingleFeature.h"
#include <cmath>
#include "LogTable.h"
#include <map>

#define FEAT_ALIGN         (1 << 0)
#define FEAT_LD            (1 << 1)
#define FEAT_MANHATTAN     (1 << 2)
#define FEAT_SQCHORD       (1 << 3)
#define FEAT_INTERSECTION  (1 << 4)
#define FEAT_PEARSON       (1 << 5)
#define FEAT_SIMRATIO      (1 << 6)
#define FEAT_N2RRC         (1 << 7)
#define FEAT_JENSONSHANNON (1 << 8)
#define FEAT_RREE_K_R      (1 << 9)
#define FEAT_KULCZYNSKI2   (1 << 10)

#define COMBO_SQUARED 1
#define COMBO_SELF    2


/*
 * Usage:
 *   add_feature(FEAT_LD | FEAT_INTERSECTION, COMBO_SELF);
 *   add_feature(FEAT_LD | FEAT_JENSONSHANNON, COMBO_SELF);
 *
 *   normalize(some_pairs_to_normalize)
 *   normalize(more_pairs_to_normalize)
 *   finalize()
 *
 *   add_feature(....);
 *
 *   normalize(some_pairs_to_normalize)
 *   normalize(more_pairs_to_normalize)
 *   finalize()
 *
 *   compute(p,q)
 *   for (size_t i = 0; i < feature.size(); i++) {
 *       cout << feature[i] << endl;
 *   }
 */
#include "LogTable.h"
template<class T>
class Feature {
public:
	Feature(const uint64_t N, const double* tbl_, double coeff_) : tbl(tbl_) {
		flags = 0;
		coeff = coeff_;
//		tbl = new LogTable(1000000, 2);
	}
	void add_feature(uint16_t f_flags, int combo=COMBO_SELF);

	void finalize();
	void remove_feature() {
		auto indices_to_rm = combos.back().second;
		combos.pop_back();
		throw "not implemented";

	}
	void normalize(const vector<pair<Point<T>*,Point<T>*> > &pairs);
	vector<double> compute(Point<T>& p, Point<T>& q) {
		vector<double> cache = compute_all_raw(p, q);
		normalize_cache(cache);
		return cache;
	};
	double operator()(int col, const vector<double>& cache) const {
		auto pr = combos.at(col);
		int combo = pr.first;
		auto indices = pr.second;
		if (combo == COMBO_SELF) {
			double prod = 1;
			for (auto idx : indices) {
				prod *= cache[idx];
			}
			return prod;
		} else if (combo == COMBO_SQUARED) {
			double prod = 1;
			for (auto idx : indices) {
				prod *= cache[idx] * cache[idx];
			}
			return prod;
		} else {
			throw "invalid combo";
		}
	}
	size_t size() const { return combos.size(); }
	void print_bounds() const {
		for (size_t i = 0; i < lookup.size(); i++) {
			cout << "bounds[" << i << "]: " << mins[i] << " to " << maxs[i] << endl;
		}
	}

	static double manhattan(Point<T>& p, Point<T>& q);
	static double length_difference(Point<T>& p, Point<T>& q);
	static double n2rrc(Point<T>& p, Point<T>& q, const vector<int>&, const vector<int> &);
	static double rree_k_r(Point<T>& p, Point<T>& q);
	static double intersection(Point<T>& p, Point<T>& q);
	double jenson_shannon(Point<T>& p, Point<T>& q) const;
	static double pearson(Point<T>& p, Point<T>& q);
	static double simratio(Point<T>& a, Point<T>& b);
	static double squaredchord(Point<T>& a, Point<T>& b);
	static double kulczynski2(Point<T>& a, Point<T>& b);
	static double align(Point<T>& a, Point<T>& b, std::map<std::pair<uintmax_t, uintmax_t>, double> &atable);
private:
	vector<double> compute_all_raw(Point<T>& p, Point<T>& q);
	void normalize_cache(vector<double>& cache) const;

	bool feat_is_sim(uint16_t single_flag) const;
        double raw(uint16_t single_flag, Point<T>& a, Point<T>& b);
	int index_of(uint16_t single_flag) const {
		for (size_t i = 0; i < lookup.size(); i++) {
			if (lookup[i] == single_flag) {
				return i;
			}
		}
		return -1;
	}

	uint16_t flags;
	std::vector<std::pair<int,
			      std::vector<int>
			      > > combos;
	std::vector<double> mins, maxs;
	std::vector<bool> is_sims, is_finalized;
	std::vector<uint16_t> lookup;
	const double* tbl;
	double coeff;


	std::map<std::pair<uintmax_t,uintmax_t>, double> atable;
};

// template<class T>
// class Feature {
// public:
// 	Feature(std::function<double(vector<double>)> combination, std::vector<SingleFeature<T> > sf)
// 		: features(sf), combo(combination) {}
// 	double operator()(Point<T>*, Point<T>*) const;


// 	static double manhattan(Point<T>& p, Point<T>& q);
// 	static double length_difference(Point<T>& p, Point<T>& q);
// 	static double n2rrc(Point<T>& p, Point<T>& q, const vector<int>&, const vector<int> &);
// 	static double rree_k_r(Point<T>& p, Point<T>& q);
// 	static double intersection(Point<T>& p, Point<T>& q);
// 	static double jenson_shannon(Point<T>& p, Point<T>& q);
// 	static double pearson(Point<T>& p, Point<T>& q);
// 	static double simratio(Point<T>& a, Point<T>& b);
// 	static double squaredchord(Point<T>& a, Point<T>& b);
// private:
// 	vector<SingleFeature<T> > features;
// 	std::function<double(vector<double>)> combo;
// };
#endif
