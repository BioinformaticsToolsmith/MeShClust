#include "Trainer.h"

#include <algorithm>
#include <set>
#include <map>
#include <cmath>
#include "../../utility/GlobAlignE.h"
#include "../../utility/AffineId.h"
#include "needleman_wunsch.h"
#include "GLM.h"
#include "Feature.h"
#include "Progress.h"
#include <random>

template<class T>
double Trainer<T>::align(Point<T> *a, Point<T>* b) const
{
	auto sa = a->get_data_str();
	auto sb = b->get_data_str();
	int la = sa.length();
	int lb = sb.length();

	// needleman_wunsch nw(sa, sb, 2, -3, 5, 2);
	// return nw.identity(nw.align());
	GlobAlignE galign(sa.c_str(), 0, la-1,
			  sb.c_str(), 0, lb-1,
			  1, -1, 2, 1);

	return galign.getIdentity();

}


template<class T>
std::tuple<Point<T>*,double,size_t,size_t> Trainer<T>::get_close(Point<T> *p, bvec_iterator<T> istart, bvec_iterator<T> iend, bool &is_min_r) const
{
	int ncols = weights.getNumRow();
#pragma omp declare reduction(pmax:std::tuple<Point<T>*,double,size_t,size_t>: \
			      omp_out = get<1>(omp_in) > get<1>(omp_out) ? omp_in : omp_out ) \
	initializer (omp_priv=std::make_tuple((Point<T>*)NULL,-1,0,0))

	std::tuple<Point<T>*,
		   double,
		   size_t,
		   size_t> result = std::tuple<Point<T>*, double, size_t, size_t>(NULL,
				     -1,
				     0,
				     0);
	bool has_found = false;

	#ifdef DEBUG
	cout << "begin " << istart.r << " " << istart.c << " end " << iend.r << " " << iend.c << endl;
	for (auto data : *istart.col) {
		cout << "\t" << data.size() << endl;
	}
	#endif
// #pragma omp parallel for reduction(pmin:result), reduction(||:has_found)
// 	for (bvec_iterator<T> i = istart; i <= iend; i++) {
// 		if (i <= iend) {
// 		Point<T>* pt = (*i).first;
// 		double sum = weights.get(0, 0);
// 		double dist = 0;
// 		for (int col = 1; col < ncols; col++) {
// 			if (col == 1) {
// 				dist = ff.at(col-1)(pt, p);
// 				sum += weights.get(col, 0) * dist;
// 			} else {
// 				sum += weights.get(col, 0) * ff.at(col-1)(pt, p);
// 			}
// 		}
// 		double res = round(1.0 / (1 + exp(-sum)));

// // set second to true if result is not 1.0
// 		// which means it will be removed
// 		result = std::make_pair(pt, dist);
// 		has_found = (res != 1.0);
// 		(*i).second = (res != 1.0);
// 		}
// 	}
	bool is_min = true;
#pragma omp parallel for reduction(pmax:result), reduction(&&:is_min)
	for (bvec_iterator<T> i = istart; i <= iend; ++i) {
		Point<T>* pt = (*i).first;
		double sum = weights.get(0, 0);
		double dist = 0;
		auto cache = feat->compute(*pt, *p);
		for (int col = 1; col < ncols; col++) {
			if (col == 1) {
				dist = (*feat)(col-1, cache);
				sum += weights.get(col, 0) * dist;
			} else {
				sum += weights.get(col, 0) * (*feat)(col-1, cache);
			}
		}
		double res = round(1.0 / (1 + exp(-sum)));
		//cout << "res: " << res << " " << dist << endl;
// set second to true if result is not 1.0
		// which means it will be removed
		result = (dist > std::get<1>(result)) ? std::make_tuple(pt, dist, i.r, i.c) : result;
		is_min = is_min && (res != 1.0);
//		has_found = has_found || (res != 1.0);
		if (res == 1.0) {
			*i = std::make_pair(pt, true);
//			(*i).second = true;
		}
	}

//	is_min = !has_found;
	is_min_r = is_min;
//	return get<0>(result);
	return result;

}

template<class T>
long Trainer<T>::merge(vector<Center<T> > &centers, long current, long begin, long last) const
{
#pragma omp declare reduction(ldpmax:std::pair<long,double>:			\
			      omp_out = omp_in.second > omp_out.second ? omp_in : omp_out ) \
	initializer (omp_priv=std::make_pair(0, std::numeric_limits<double>::min()))
	std::pair<long,double> best = std::make_pair(0, std::numeric_limits<double>::min());
	Point<T>* p = centers[current].getCenter();
#pragma omp parallel for reduction(ldpmax:best)
	for (long i = begin; i <= last; i++) {
		double sum = weights.get(0, 0);
		double dist = 0;
		Point<T>* cen = centers[i].getCenter();
		auto cache = feat->compute(*cen, *p);
		for (int col = 1; col < weights.getNumRow(); col++) {
			double d = (*feat)(col-1, cache);
			if (col == 1) {
				dist = d;
			}
			sum += weights.get(col, 0) * d;
		}
		double res = round(1.0 / (1 + exp(-sum)));

		if (res == 1) {
			best = best.second > dist ? best : std::make_pair(i, dist);
		}
	}
	return best.first;
}

template<class T>
vector<pair<Point<T>*,Point<T>*> > resize_vec(vector<pair<pair<Point<T>*,Point<T>*>, double> > &vec, size_t new_size)
{
	cout << "Vector size: " << vec.size() << " min size: " << new_size << endl;
	vector<pair<Point<T>*, Point<T>*> > data;
	if (vec.size() <= new_size) {
		for (int i = 0; i < vec.size(); i++) {
			data.push_back(vec[i].first);
		}
		return data;
	}
	using k = pair<pair<Point<T>*,Point<T>*>, double>;
	std::sort(vec.begin(), vec.end(), [](const k& a, const k& b) {
			return a.second < b.second;
		});
	double interval = (double)vec.size() / (vec.size() - new_size);
	std::set<int> indices;
	int i = 0;
	for (double index = 0; round(index) < vec.size() && i < (vec.size() - new_size);
	     i++, index += interval) {
		int j = round(index);
		indices.insert(j);
	}

	std::cout << "index size: " << indices.size() << std::endl;

	// for (double index = 0; round(index) < vec.size() && indices.size() < new_size;
	//      index += interval) {
	// 	int j = round(index);
	// 	indices.insert(vec[j]);
	// }
	// vec.erase(vec.begin(), std::remove_if(vec.begin(), vec.end(), [&](const k& a) {
	// 			return indices.find(a) == indices.end();
	// 		}));
	for (auto iter = indices.rbegin(); iter != indices.rend(); iter++) {
		int idx = *iter;
		vec.erase(vec.begin() + idx);
	}
	if (vec.size() != new_size) {
		cerr << "sizes are not the same: " << vec.size() << " " << new_size <<  endl;
		throw "Resize did not work";
	}
	for (auto a : vec) {
		data.push_back(a.first);
	}
	return data;
}

struct rng {
	rng() {
		srand(0);
	}
	int operator()(int n) const {
		return rand() % n;
	}
};
template<class T>
	pair<vector<pair<Point<T>*,
			 Point<T>*
			 > >,
	     vector<pair<Point<T>*,
			 Point<T>*> > > Trainer<T>::get_labels(vector<pair<Point<T>*,Point<T>*> > &vec, double cutoff) const
{

	auto cmp = [](const pair<Point<T>*,Point<T>*> a, const pair<Point<T>*,Point<T>*> b) {
		return a.first->get_header().compare(b.first->get_header()) < 0
		||
		(a.first->get_header() == b.first->get_header() && a.second->get_header().compare(b.second->get_header()) < 0);
	};
	auto scmp = [](const pair<pair<Point<T>*,Point<T>*>,double> a, const pair<pair<Point<T>*,Point<T>*>, double> b) {
		return a.first.first->get_header().compare(b.first.first->get_header()) < 0
		||
		(a.first.first->get_header() == b.first.first->get_header() && a.first.second->get_header().compare(b.first.second->get_header()) < 0);
	};

	// todo: convert to std::map
	std::set<pair<pair<Point<T>*,Point<T>*>, double>, decltype(scmp)> buf_pos(scmp), buf_neg(scmp);
	std::vector<pair<pair<Point<T>*,Point<T>*>, double> > buf_vpos, buf_vneg;
//	std::sort(vec.begin(), vec.end(), cmp);
	// cout << "Before Pair: " << vec[0].first->get_header() << ", " << vec[0].second->get_header() << endl;
	// cout << "Before Pair: " << vec[vec.size()-1].first->get_header() << ", " << vec[vec.size()-1].second->get_header() << endl;

	rng gen;
	random_shuffle(vec.begin(), vec.end(), gen);
	// cout << "Pair: " << vec[0].first->get_header() << ", " << vec[0].second->get_header() << endl;
	// cout << "Pair: " << vec[vec.size()-1].first->get_header() << ", " << vec[vec.size()-1].second->get_header() << endl;
	vector<double> scores(vec.size());
	Progress p(vec.size(), "Alignment");
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < vec.size(); i++) {
		double algn = align(vec[i].first, vec[i].second);
		bool is_pos = algn >= cutoff;
#pragma omp critical
		{
			scores[i] = algn;
			p++;
			if (is_pos) {
				buf_pos.insert(make_pair(vec[i], algn));
				//cout << vec[i].first->get_header() << " " << vec[i].second->get_header() << " " << algn << endl;
			} else {
				buf_neg.insert(make_pair(vec[i], algn));
			}

#ifdef DEBUG
			cout << vec[i].first->get_header() << " WITH " << vec[i].second->get_header() << " " << algn << endl;
			#endif

		}
	}
	p.end();
	std::sort(scores.begin(), scores.end());
	std::cout << "positive=" << buf_pos.size() << " negative=" << buf_neg.size() << endl;
	if (buf_pos.empty() || buf_neg.empty()) {
		std::cout << "Identity value does not match sampled data: ";
		if (buf_pos.empty()) {
			std::cout << "Too many sequences below identity";
		} else {
			std::cout << "Too many sequences above identity";
		}
		std::cout << std::endl;
		exit(0);
	}
	size_t m_size = std::min(buf_pos.size(), buf_neg.size());

	std::cout << "resizing positive" << std::endl;
	for (auto p : buf_pos) {
		buf_vpos.push_back(p);
	}
	for (auto p : buf_neg) {
		buf_vneg.push_back(p);
	}
	auto bp = resize_vec(buf_vpos, m_size);
	std::cout << "resizing negative" << std::endl;
	auto bn = resize_vec(buf_vneg, m_size);
        auto ret = make_pair(bp, bn);
	std::cout << "positive=" << ret.first.size() << " negative=" << ret.second.size() << endl;
	return ret;

}
template<class T>
void Trainer<T>::filter(Point<T> *p, vector<pair<Point<T> *, bool> > &vec) const
{
	for (auto& pt : vec) {
		double sum = weights.get(0, 0);
		auto cache = feat->compute(*pt.first, *p);
		for (int col = 1; col < weights.getNumRow(); col++) {
			sum += weights.get(col, 0) * (*feat)(col-1, cache);
		}
		double res = round(1.0 / (1 + exp(-sum)));
		pt.second = (res != 1);
	}
	vec.erase(std::remove_if(vec.begin(), vec.end(), [](pair<Point<T>*, bool> p) {
				return p.second;
			}), vec.end());
}

template<class T>
Point<T>* Trainer<T>::closest(Point<double> *p, vector<pair<Point<T> *, bool> > &vec) const
{
	Point<T>* best_pt = NULL;
	double best_dist = 0;
	for (auto& pt : vec) {
		double sum = weights.get(0, 0);
		double dist = pt.first->distance_d(*p);
		if (best_pt == NULL || dist < best_dist) {
			best_dist = dist;
			best_pt = pt.first;
		}
	}
	return best_pt;
}

template<class T>
std::pair<matrix::Matrix,matrix::Matrix> Trainer<T>::generate_feat_mat(pair<vector<pair<Point<T> *, Point<T> *> >, vector<pair<Point<T> *, Point<T> *> > > &data, int ncols)
{
	int nrows = data.first.size() + data.second.size();
	matrix::Matrix feat_mat(nrows, ncols);
	matrix::Matrix labels(nrows, 1);
#pragma omp parallel for
	for (int i = 0; i < data.first.size(); i++) {
		auto kv = data.first[i];
		int row = i;
		auto cache = feat->compute(*kv.first, *kv.second);
		for (int col = 0; col < ncols; col++) {

			if (col == 0) {
				feat_mat.set(row, col, 1);
			} else {
//				double val = ff[col-1](kv.first, kv.second);
				////#pragma omp critical
				double val = (*feat)(col-1, cache);
				feat_mat.set(row, col, val);
			}

		}
		////#pragma omp critical
		labels.set(row, 0, 1);
	}
#pragma omp parallel for
	for (int i = 0; i < data.second.size(); i++) {
		auto kv = data.second[i];
		int row = data.first.size() + i;
		auto cache = feat->compute(*kv.first, *kv.second);
		for (int col = 0; col < ncols; col++) {

			if (col == 0) {
				feat_mat.set(row, col, 1);
			} else {
//				double val = ff[col-1](kv.first, kv.second);
				////#pragma omp critical
				double val = (*feat)(col-1, cache);
				feat_mat.set(row, col, val);
			}

		}
		////#pragma omp critical
		labels.set(row, 0, -1);
	}
	return std::make_pair(feat_mat, labels);
}
template<class T>
double Trainer<T>::train_n(pair<vector<pair<Point<T> *, Point<T> *> >, vector<pair<Point<T> *, Point<T> *> > > &data, int ncols)
{
	std::cout << "done" << endl;
	cout << "Training on " << ncols << " columns" << endl;
	int nrows = data.first.size() + data.second.size();

	matrix::Matrix feat_mat(nrows, ncols);
	matrix::Matrix labels(nrows, 1);
	double avg_label = 0;
#pragma omp parallel for
	for (int i = 0; i < data.first.size(); i++) {
		auto kv = data.first[i];
		int row = i;
		auto cache = feat->compute(*kv.first, *kv.second);
		for (int col = 0; col < ncols; col++) {

			if (col == 0) {
				feat_mat.set(row, col, 1);
			} else {
//				double val = ff[col-1](kv.first, kv.second);
				////#pragma omp critical
				double val = (*feat)(col-1, cache);
				feat_mat.set(row, col, val);
			}

		}
		////#pragma omp critical
		labels.set(row, 0, 1);
	}
#pragma omp parallel for
	for (int i = 0; i < data.second.size(); i++) {
		auto kv = data.second[i];
		int row = data.first.size() + i;
		auto cache = feat->compute(*kv.first, *kv.second);
		for (int col = 0; col < ncols; col++) {

			if (col == 0) {
				feat_mat.set(row, col, 1);
			} else {
//				double val = ff[col-1](kv.first, kv.second);
				////#pragma omp critical
				double val = (*feat)(col-1, cache);
				feat_mat.set(row, col, val);
			}

		}
		////#pragma omp critical
		labels.set(row, 0, -1);
	}
	for (int row = 0; row < nrows; row++) {
		for (int col = 0; col < ncols; col++) {
			double val = feat_mat.get(row, col);
			std::cout << val << "\t";
		}
		std::cout << endl;
	}
	glm.train(feat_mat, labels);
	weights = glm.get_weights();
	#ifdef DEBUG
	for (int i = 0; i < ncols; i++) {
		cout << "weight: " << weights.get(i, 0) << endl;

	}
	#endif
	matrix::Matrix p = glm.predict(feat_mat);
	for (int row = 0; row < nrows; row++) {
		if (p.get(row, 0) == 0) {
			p.set(row, 0, -1);
		}
	}
	auto tup = glm.accuracy(labels, p);
	return get<0>(tup);
}

template<class T>
void Trainer<T>::train(double acc_cutoff)
{
	pair<vector<pair<Point<T>*,
			 Point<T>*
			 > >,
	     vector<pair<Point<T>*,
			 Point<T>*> > > training, testing;
	if (k != 0) {
	std::cout << "Splitting data" << endl;
	auto _data = split();
// 	std::sort(_data.begin(), _data.end(), [](const pair<Point<T>*,Point<T>*> a, const pair<Point<T>*,Point<T>*> b) {
// 			return a.first->get_header().compare(b.first->get_header()) < 0
// ||
// 										      (a.first->get_header() == b.first->get_header() && a.second->get_header().compare(b.second->get_header()) > 0);
// 		});
	auto both = get_labels(_data, cutoff);

	for (int i = 0; i < both.first.size(); i++) {
		//training.first.push_back(both.first[i]);
		if (i % 2 == 0) {
			training.first.push_back(both.first[i]);
		} else {
			testing.first.push_back(both.first[i]);
		}
	}
	both.first.clear();
	for (int i = 0; i < both.second.size(); i++) {
		if (i % 2 == 0) {
			training.second.push_back(both.second[i]);
		} else {
			testing.second.push_back(both.second[i]);
		}
	}
	both.second.clear();
	if (testing.first.empty() || testing.second.empty()) {
		throw "not enough points to sample";
	}
	}
	vector<std::pair<uint16_t, uint16_t> > bit_feats;
	//bit_feats.push_back(std::make_pair(FEAT_ALIGN, COMBO_SELF));
	int feat_set = 1;
	if (k == 0) {
		feat->add_feature(FEAT_ALIGN, COMBO_SELF);
		feat->normalize(training.first);
		feat->finalize();
		weights = matrix::Matrix(2, 1);
		weights.set(0, 0, -1 * cutoff);
		weights.set(1, 0, 1);
		return;
	} else if (feat_set == 0) {
		bit_feats.push_back(std::make_pair(FEAT_LD | FEAT_INTERSECTION, COMBO_SELF));
		bit_feats.push_back(std::make_pair(FEAT_LD | FEAT_JENSONSHANNON, COMBO_SELF));
		bit_feats.push_back(std::make_pair(FEAT_SIMRATIO, COMBO_SELF));
		bit_feats.push_back(std::make_pair(FEAT_SQCHORD, COMBO_SELF));
	} else {
		bit_feats.push_back(std::make_pair(FEAT_INTERSECTION | FEAT_LD, COMBO_SELF));
		bit_feats.push_back(std::make_pair(FEAT_MANHATTAN | FEAT_LD, COMBO_SQUARED));
		bit_feats.push_back(std::make_pair(FEAT_PEARSON, COMBO_SELF));
		bit_feats.push_back(std::make_pair(FEAT_KULCZYNSKI2 | FEAT_LD, COMBO_SQUARED));
	}
// 	std::cout << "data split, normalizing..." << endl;
// #pragma omp parallel for schedule(dynamic)
// 	for (int i = 0; i < sf.size(); i++) {
// 		sf[i].normalize(data.first);
// 		sf[i].normalize(data.second);
// 		cout << i << ": " << sf[i].max << " " << sf[i].min << endl;
// 	}
	// ff.emplace_back([](vector<double> d) { return d[0] * d[1]; },
	//   		std::vector<SingleFeature<T> >({sf[2], sf[1]}));
	// ff.emplace_back([](vector<double> d) { return d[0] * d[0] * d[1] * d[1]; },
	// 		std::vector<SingleFeature<T> >({sf[0], sf[1]}));
	// ff.emplace_back([](vector<double> d) { return d[0]; },
	//  		std::vector<SingleFeature<T> >({sf[3]}));

	double prev_acc = -10000;
	vector<matrix::Matrix> matvec;
	vector<Feature<T> > features;
	const size_t min_no_features = std::max(1, (int)bit_feats.size()-1);
	for (size_t num_features = min_no_features; num_features <= bit_feats.size(); num_features++) {
		for (size_t j = feat->size(); j < num_features && j < bit_feats.size(); j++) {
			feat->add_feature(bit_feats[j].first, bit_feats[j].second);
		}
		feat->normalize(training.first);
		feat->normalize(training.second);
		feat->finalize();
		feat->print_bounds();
		auto mtraining = generate_feat_mat(training, num_features+1);
		auto mtesting = generate_feat_mat(testing, num_features+1);
		glm.train(mtraining.first, mtraining.second);
		weights = glm.get_weights();
		matrix::Matrix p = glm.predict(mtesting.first);
		for (int row = 0; row < testing.first.size() + testing.second.size(); row++) {
			if (p.get(row, 0) == 0) {
				p.set(row, 0, -1);
			}
		}
		double acc = get<0>(glm.accuracy(mtesting.second, p));
		matrix::Matrix q = glm.predict(mtraining.first);
		for (int row = 0; row < training.first.size() + training.second.size(); row++) {
			if (q.get(row, 0) == 0) {
				q.set(row, 0, -1);
			}
		}
		glm.accuracy(mtraining.second, q);
		if (acc - prev_acc <= 1 && acc >= 90.0) {
			weights = matvec.back();
			*feat = features.back();
			cout << "feat size is " << feat->size() << endl;
			break;
		}
		matvec.push_back(weights);
		features.push_back(*feat);
		prev_acc = acc;
		if (acc >= acc_cutoff) {
			cout << "breaking from acc cutoff" << endl;
			break;
		}
	}
	cout << "Final: feat size is " << feat->size() << endl;
	cout << "Using " << weights.getNumRow()-1 << " features " << __DATE__ << endl;
}

template<class T>
vector<pair<Point<T>*, Point<T>*> > Trainer<T>::split()
{
	// n_points total per side
	// max_pts_from_one on each side
	auto cmp = [](const pair<Point<T>*,Point<T>*> a, const pair<Point<T>*,Point<T>*> b) {
			return a.first->get_header().compare(b.first->get_header()) < 0
||
										      (a.first->get_header() == b.first->get_header() && a.second->get_header().compare(b.second->get_header()) < 0);
	};
        set<pair<Point<T>*, Point<T>*>, decltype(cmp)> pairs(cmp);
//	vector<pair<Point<T>*, Point<T>*> > pairs;
	const size_t total_num_pairs = n_points * 2;
	int aerr = 0;
	int bandwidth = (1.0 - cutoff) * 10000;
	vector<Point<T>*> indices;
	std::sort(points.begin(), points.end(), [](const Point<T>* a,
						   const Point<T>* b) -> bool {
			  return a->get_length() < b->get_length();
			  });
	Point<T> *begin_pt = points[points.size()/2];

	std::sort(points.begin(), points.end(), [&](const Point<T>* a,
							    const Point<T>* b) -> bool {
				  return a->distance(*begin_pt) < b->distance(*begin_pt);
			  });
	int num_iterations = ceil(((double)n_points) / max_pts_from_one) - 1;
	for (int i = 0; i <= num_iterations; i++) {
		int idx = i * (points.size()-1) / num_iterations;
		indices.push_back(points[idx]);
	}
	cout << "Point pairs: " << indices.size() << endl;
	size_t to_add_each = max_pts_from_one / 2;
	Progress prog(indices.size(), "Sorting data");
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < indices.size(); i++) {
		vector<Point<T>*> pts = points;
		Point<T>* p = indices[i];
		std::sort(pts.begin(), pts.end(), [&](const Point<T>* a,
						      const Point<T>* b) {
				  return a->distance(*p) < b->distance(*p);
			  });
		// do binary search with alignment
		size_t offset = pts.size() / 4;
		size_t pivot = offset;
		double closest_algn = 20000;
		size_t best_pivot = 2 * offset;
		for (pivot = 2 * offset; offset > 0; offset /= 2) {
			double algn = align(p, pts[pivot]);
			// cout << "Pivot: " << pivot << " point: " << pts[pivot]->get_header() << " sim: " << align(p, pts[pivot]) << endl;
			if (fabs(algn - cutoff) < closest_algn) {
				closest_algn = fabs(algn - cutoff);
				best_pivot = pivot;
			}
			if (algn < cutoff) {
				pivot -= offset;
			} else if (algn > cutoff) {
				pivot += offset;
			} else {
				break;
			}
		}
//		cout << "Pivot: " << pivot << " point: " << pts[pivot]->get_header() << " sim: " << align(p, pts[pivot]) << endl;
		// before: [0, pivot) size: to_add_each
		// after: [pivot, size) size: to_add_each
		double before_inc = (double)pivot / to_add_each;
		double after_inc = ((double)(pts.size() - pivot)) / to_add_each;
#pragma omp critical
		{
			prog++;
			if (before_inc < 1) {
				aerr = 1;
			} else if (after_inc < 1) {
				aerr = -1;
			}
		}
		double before_start = 0;
		double after_start = pivot;
		double top_start = 0;
		size_t size_before = pairs.size();
		vector<pair<Point<T>*,Point<T>*> > buf;
		// Adds points above cutoff by adding before_inc
		for (int i = 0; i < to_add_each; i++) {
			int idx = round(before_start);
			int dist = pts[idx]->distance(*p);
			//	cout << p->get_header() << " " << pts[idx]->get_header() << " " << dist << endl;
			auto pr = p->get_header().compare(pts[idx]->get_header()) < 0 ? make_pair(p, pts[idx]) : make_pair(pts[idx], p);
			buf.push_back(pr);
			before_start += before_inc;
		}
		// Adds points before cutoff by adding after_inc
		for (int i = 0; i < to_add_each && round(after_start) < pts.size(); i++) {
			int idx = round(after_start);
			int dist = pts[idx]->distance(*p);
			//		cout << p->get_header() << " " << pts[idx]->get_header() << " " << dist << endl;
			auto pr = p->get_header().compare(pts[idx]->get_header()) < 0 ? make_pair(p, pts[idx]) : make_pair(pts[idx], p);
			buf.push_back(pr);
			after_start += after_inc;
		}
#pragma omp critical
		{
			// Adds buffer to total pairs
		// 	for (auto p : buf) {
// 				pairs.push_back(p);
// 			}
			pairs.insert(std::begin(buf), std::end(buf));
		}
//			cout << "added " << pairs.size() - size_before << " pairs" << endl;
	}
	prog.end();
	if (aerr < 0) {
		cerr << "Warning: Alignment may be too small for sampling" << endl;
	} else if (aerr > 0) {
		cerr << "Warning: Alignment may be too large for sampling" << endl;
	}
	int i = 0;
	for (auto a : pairs) {
		cout << "Before Pair: " << a.first->get_header() << ", " << a.second->get_header() << endl;
		if (++i == 4) {
			break;
		}
	}
	return std::vector<std::pair<Point<T>*,Point<T>*> >(pairs.begin(), pairs.end());
}
template<class T>
std::pair<std::map<std::pair<Point<T>*, Point<T>*>, double>,
	  std::map<std::pair<Point<T>*, Point<T>*>, double> >
Trainer<T>::split_old() {
	using train_map = std::map<std::pair<Point<T>*, Point<T>*>, double>;
	std::pair<train_map, train_map> split;
	int bandwidth = (1.0 - cutoff) * 10000;
	size_t last_cutoff = points.size() / 2;
	while (split.first.size() < n_points) {
		Point<T> *p = points[last_cutoff];
		std::sort(points.begin(), points.end(), [&](const Point<T>* a,
							    const Point<T>* b) -> bool {
				  return a->distance(*p) < b->distance(*p);
			  });
		int b_cutoff = points.size() / 2;
		for (int offset = b_cutoff; offset >= 1; offset /= 2) {
			int dist = p->distance(*points[b_cutoff]);
			if (dist < bandwidth) {
				b_cutoff += offset;
			} else if (dist > bandwidth) {
				b_cutoff -= offset;
			} else {
				break;
			}
		}
		size_t cutoff_index = points.size();
		const size_t count = split.first.size();

		if (b_cutoff >= max_pts_from_one) {
			double ratio = (double)b_cutoff / max_pts_from_one;
			double sum = 0;
			for (size_t q = 0; q < max_pts_from_one; q++) {
				size_t i = round(sum);
				if (i >= points.size()) {
					cerr << "this shouldn't happen" << endl;
					throw "this shouldn't happen";
				}
				double alignment = align(p, points[i]);
				if (alignment < cutoff) {
					cutoff_index = i + 10;
					break;
				}
				if (split.first.size() < n_points) {
					split.first[make_pair(p, points[i])] = alignment;
				}
				sum += ratio;
			}
		} else {
			for (size_t i = 1; i < cutoff_index; i++) {
				double alignment = align(p, points[i]);
				if (alignment < cutoff) {
					cutoff_index = i + 10;
					break;
				}
				if (split.first.size() < n_points) {
					split.first[make_pair(p, points[i])] = alignment;
				}
			}
		}
		size_t similar_points_added = split.first.size() - count;
		size_t available_points = points.size() - cutoff_index;
		if (available_points == 0 || available_points <= similar_points_added) {
			cerr << "change cutoff value, points are too similar" << endl;
			throw "change cutoff value, points are too similar";
		}
		double ratio = (double)(available_points - 1.0) / (double)similar_points_added;
		double sum = 0;
		for (size_t q = 0; q < similar_points_added; q++) {
			size_t i = cutoff_index + round(sum);
			if (i >= points.size()) {
				break;
			}
			double alignment = align(p, points[i]);
			split.second[make_pair(p, points[i])] = alignment;
			sum += ratio;
		}
	        if (split.first.size() != split.second.size()) {
			cerr << "something happened";
			throw "something happened";
		}
		last_cutoff = cutoff_index;
	}
	for (auto p : points) {
		p->set_data_str("");
	}
	return split;
}


int gcd(int a, int b)
{
	if (b <= 0) {
		return a;
	}
	return gcd(b, a % b);
}
int gcd_vec(std::vector<int> v)
{
	int ret = v[0];
	for (size_t i = 1; i < v.size(); i++) {
		if (v[i] == 0) {
			continue;
		}
		ret = gcd(ret, v[i]);
	}
	return ret;
}

inline int sign(double x) {
	return (x > 0) - (x < 0);
}
void scale(double (&mat)[4][4], double &sigma, double& epsilon)
{
	double scale_factor = 100000;
	std::vector<int> signs, scaled;
	signs.push_back(sign(sigma));
	scaled.push_back(round(scale_factor * fabs(sigma)));
	signs.push_back(sign(epsilon));
	scaled.push_back(round(scale_factor * fabs(epsilon)));
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			signs.push_back(sign(mat[i][j]));
			scaled.push_back(round(scale_factor * fabs(mat[i][j])));
		}
	}
	double common_div = gcd_vec(scaled);
	sigma = signs[0] * scaled[0] / common_div;
	epsilon = signs[1] * scaled[1] / common_div;
	int count = 2;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mat[i][j] = signs[count] * scaled[count] / common_div;
			count++;
		}
	}
}

template<class T>
void Trainer<T>::init(double (&matrix)[4][4], double sig, double eps)
{
	scale(matrix, sig, eps);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mat[i][j] = (int)matrix[i][j];
		}
	}
	sigma = (int)sig;
	eps = (int)eps;
	// sf.emplace_back([](Point<T>* a, Point<T> *b) {
	// 	        return Feature<T>::manhattan(*a, *b);
	// 	}, false);
	// sf.emplace_back([](Point<T>* a, Point<T> *b) {
	// 		return Feature<T>::length_difference(*a, *b);
	// 	}, false);
	// sf.emplace_back([](Point<T>* a, Point<T> *b) {
	// 		return Feature<T>::rree_k_r(*a, *b);
	// 	}, false);
	// sf.emplace_back([](Point<T>* a, Point<T>* b) {
	// 		return Feature<T>::length_difference(*a, *b);
	// 	}, false);
	// sf.emplace_back([](Point<T>* a, Point<T>* b) {
	// 		return Feature<T>::intersection(*a, *b);
	// 	}, true);
	// sf.emplace_back([](Point<T>* a, Point<T>* b) {
	// 		return Feature<T>::jenson_shannon(*a, *b);
	// 	}, false);
	// sf.emplace_back([](Point<T>* a, Point<T>* b) {
	// 		return Feature<T>::simratio(*a, *b);
	// 	}, true);
	// sf.emplace_back([](Point<T>* a, Point<T>* b) {
	// 		return Feature<T>::squaredchord(*a, *b);
	// 	}, false);
	// sf.emplace_back([](Point<T>* a, Point<T>* b) {
	// 		return Feature<T>::manhattan(*a, *b);
	// 	}, false);
	// sf.emplace_back([](Point<T>* a, Point<T>* b) {
	// 		return Feature<T>::pearson(*a, *b);
	// 	}, true);

}
template class Trainer<uint8_t>;
template class Trainer<uint16_t>;
template class Trainer<uint32_t>;
template class Trainer<uint64_t>;
template class Trainer<int>;
template class Trainer<double>;
