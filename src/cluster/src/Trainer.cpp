#include "Trainer.h"

#include <algorithm>
#include <set>
#include <map>
#include <cmath>
#include "needleman_wunsch.h"
#include "GLM.h"
#include "Feature.h"
template<class T>
double Trainer<T>::align(Point<T> *a, Point<T>* b) const
{
	auto sa = a->get_data_str();
	auto sb = b->get_data_str();
	int la = sa.length();
	int lb = sb.length();
	needleman_wunsch nw(sa, sb, mat, sigma, epsilon);
	return nw.identity(nw.align());
}

template<class T>
Point<T>* Trainer<T>::merge(Point<T> *p, vector<pair<Point<T> *, double> > &vec) const
{
	size_t idx = 0;
	for (auto& pt : vec) {
		double sum = weights.get(0, 0);
		double dist;
		for (int col = 1; col < weights.getNumRow(); col++) {
			double d = ff[col-1](pt.first, p);
			if (col == 1) {
				dist = d;
			}
			sum += weights.get(col, 0) * d;
		}
		double res = round(1.0 / (1 + exp(-sum)));
		pt.second = (res == 1) ? dist : -1;
	}
	auto iter = std::max_element(vec.begin(), vec.end(), [](pair<Point<T>*, double> a, pair<Point<T>*, double> b) {
			return a.second < b.second;
		});
	if (iter != vec.end() && iter->second != -1) {
		return iter->first;
	} else {
		return NULL;
	}
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
template<class T>
	pair<vector<pair<Point<T>*,
			 Point<T>*
			 > >,
	     vector<pair<Point<T>*,
			 Point<T>*> > > Trainer<T>::get_labels(vector<pair<Point<T>*,Point<T>*> > &vec, double cutoff) const
{
	std::vector<pair<pair<Point<T>*,Point<T>*>, double> > buf_pos, buf_neg;
	random_shuffle(vec.begin(), vec.end());
	vector<double> scores(vec.size());
#pragma omp parallel for
	for (int i = 0; i < vec.size(); i++) {
		double algn = align(vec[i].first, vec[i].second);
		bool is_pos = algn >= cutoff;
#pragma omp critical
		{
			scores[i] = algn;
			if (is_pos) {
				buf_pos.push_back(make_pair(vec[i], algn));
			} else {
				buf_neg.push_back(make_pair(vec[i], algn));
			}
		}
	}
	std::sort(scores.begin(), scores.end());
	for (double algn : scores) {
		cout << "alignment: " << algn << endl;
	}
	std::cout << "positive=" << buf_pos.size() << " negative=" << buf_neg.size() << endl;

	size_t m_size = std::min(buf_pos.size(), buf_neg.size());

	std::cout << "resizing positive" << std::endl;
	auto bp = resize_vec(buf_pos, m_size);
	std::cout << "resizing negative" << std::endl;
	auto bn = resize_vec(buf_neg, m_size);
        auto ret = make_pair(bp, bn);
	std::cout << "positive=" << ret.first.size() << " negative=" << ret.second.size() << endl;
	return ret;

}
template<class T>
void Trainer<T>::filter(Point<T> *p, vector<pair<Point<T> *, bool> > &vec) const
{
	for (auto& pt : vec) {
		double sum = weights.get(0, 0);
		for (int col = 1; col < weights.getNumRow(); col++) {
			sum += weights.get(col, 0) * ff[col-1](pt.first, p);
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
vector<pair<int, double> > Trainer<T>::get_close(Point<T> *p, const vector<pair<Point<T> *, int> > &vec, bool &is_min) const
{

	int nrows = vec.size();
	int ncols = weights.getNumRow();
	vector<pair<int, double> > close(nrows);
	vector<double> min_distances(nrows);
#pragma omp parallel for
	for (int row = 0; row < nrows; row++) {
		double sum = weights.get(0, 0);
		double dist = 0;
		for (int col = 1; col < ncols; col++) {
			if (col == 1) {
				dist = ff.at(col-1)(vec[row].first, p);
				sum += weights.get(col, 0) * dist;
			} else {
				sum += weights.get(col, 0) * ff.at(col-1)(vec[row].first, p);
			}
		}
		double res = round(1.0 / (1 + exp(-sum)));
		int idx = res == 1 ? vec[row].second : -1;
		close[row] = make_pair(idx, dist);
		min_distances[row] = dist;
	}
	is_min = false;
	close.erase(std::remove_if(close.begin(), close.end(), [](pair<int,double> p) {
				return p.first == -1;
			}), close.end());
	if (close.empty() && !vec.empty()) {
		auto iter = std::max_element(min_distances.begin(), min_distances.end());
		int min_manhattan_idx = std::distance(min_distances.begin(), iter);
		//	std::cout << "using min: [" << min_manhattan_idx << "] " << min_manhattan << endl;
		close.push_back(make_pair(vec[min_manhattan_idx].second, 0));
		is_min = true;
	}

	std::sort(close.begin(), close.end(), [](pair<int,double> a, pair<int,double> b) {
				return a.second > b.second;
		});
//	std::cout << "size(): " << close.size() << " " << is_min << std::endl;
	return close;
// 	// code: labels = features * weights
// 	matrix::Matrix feat_mat(nrows, ncols);
// 	auto begin = clock();
// #pragma omp parallel for collapse(2)
// 	for (int row = 0; row < nrows; row++) {
// 		for (int col = 0; col < ncols; col++) {
// 			if (col == 0) {
// 				feat_mat.set(row, col, 1);
// 			} else {
// 				double val = ff[col-1](vec[row].first, p);
// 				feat_mat.set(row, col, val);
// 			}
// 		}
// 	}
// 	auto end = clock();
// 	double diff = ((double)end - begin) / CLOCKS_PER_SEC;
// //	std::cout << "Loop diff: " << diff << endl;
// 	matrix::Matrix pred = glm.predict(feat_mat);

// //	std::cout << "close to " << p->get_header() << endl;
// 	for (int row = 0; row < nrows; row++) {
// 		if (pred.get(row, 0) == 1) {
// 			//std::cout << vec[row].first->get_header() << endl;
// 			close.push_back(vec[row].second);
// 		}
// 	}
// 	return close;
}
template<class T>
double Trainer<T>::train_n(pair<vector<pair<Point<T> *, Point<T> *> >, vector<pair<Point<T> *, Point<T> *> > > &data, int ncols)
{
	std::cout << "done" << endl;
	int nrows = data.first.size() + data.second.size();

	matrix::Matrix feat_mat(nrows, ncols);
	matrix::Matrix labels(nrows, 1);
	double avg_label = 0;
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < data.first.size(); i++) {
		auto kv = data.first[i];
		int row = i;
		for (int col = 0; col < ncols; col++) {
			if (col == 0) {
				feat_mat.set(row, col, 1);
			} else {
				double val = ff[col-1](kv.first, kv.second);
				feat_mat.set(row, col, val);
			}
		}
		labels.set(row, 0, 1);
	}
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < data.second.size(); i++) {
		auto kv = data.second[i];
		int row = data.first.size() + i;
		for (int col = 0; col < ncols; col++) {
			if (col == 0) {
				feat_mat.set(row, col, 1);
			} else {
				double val = ff[col-1](kv.first, kv.second);
				feat_mat.set(row, col, val);
			}
		}
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
	for (int i = 0; i < ncols; i++) {
		cout << "weight: " << weights.get(i, 0) << endl;

	}
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
	std::cout << "Splitting data" << endl;
	auto _data = split();
	auto data = get_labels(_data, cutoff);
	if (data.first.empty() || data.second.empty()) {
		throw "not enough points to sample";
	}
	std::cout << "data split, normalizing..." << endl;
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < sf.size(); i++) {
		sf[i].normalize(data.first);
		sf[i].normalize(data.second);
#pragma omp critical
		cout << i << ": " << sf[i].max << " " << sf[i].min << endl;
	}
	// ff.emplace_back([](vector<double> d) { return d[0] * d[1]; },
	//   		std::vector<SingleFeature<T> >({sf[2], sf[1]}));
	// ff.emplace_back([](vector<double> d) { return d[0] * d[0] * d[1] * d[1]; },
	// 		std::vector<SingleFeature<T> >({sf[0], sf[1]}));
	// ff.emplace_back([](vector<double> d) { return d[0]; },
	//  		std::vector<SingleFeature<T> >({sf[3]}));
	auto prod2 = [](vector<double> d) { return d[0] * d[1]; };
	auto pself = [](vector<double> d) { return d[0]; };
	auto prod2sq = [](vector<double> d) { return d[0] * d[1] * d[0] * d[1]; };
	ff.emplace_back(prod2, std::vector<SingleFeature<T> >({sf[0], sf[1]}));
	ff.emplace_back(prod2, std::vector<SingleFeature<T> >({sf[0], sf[2]}));
	ff.emplace_back(pself, std::vector<SingleFeature<T> >({sf[3]}));
	ff.emplace_back(pself, std::vector<SingleFeature<T> >({sf[4]}));

	double prev_acc = 0;
	vector<matrix::Matrix> matvec;
	for (size_t i = 0; i < ff.size(); i++) {
		double acc = train_n(data, i + 2);

		if (acc - prev_acc <= 1 && acc >= 90.0) {
			weights = matvec.back();
			break;
		}
		matvec.push_back(weights);
		prev_acc = acc;
		if (acc >= acc_cutoff) {
			break;
		}
	}
	cout << "Using " << weights.getNumRow()-1 << " features " << endl;
}

template<class T>
vector<pair<Point<T>*, Point<T>*> > Trainer<T>::split()
{
	// n_points total per side
	// max_pts_from_one on each side
	vector<pair<Point<T>*, Point<T>*> > pairs;
	const size_t total_num_pairs = n_points * 2;

	int bandwidth = (1.0 - cutoff) * 10000;
	cout << "bandwidth: " << bandwidth << endl;
	vector<Point<T>*> indices;
	Point<T> *begin_pt = points[0];
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
		if (before_inc < 1) {
			cerr << "Too large similarity" << endl;
		} else if (after_inc < 1) {
			cerr << "Too small similarity" << endl;
		}
		double before_start = 0;
		double after_start = pivot;
		double top_start = 0;
		size_t size_before = pairs.size();
		vector<pair<Point<T>*,Point<T>*> > buf;
		for (int i = 0; i < to_add_each; i++) {
			int idx = round(before_start);
			int dist = pts[idx]->distance(*p);
			//	cout << p->get_header() << " " << pts[idx]->get_header() << " " << dist << endl;
			buf.push_back(make_pair(p, pts[idx]));
			before_start += before_inc;
		}
		for (int i = 0; i < to_add_each && round(after_start) < pts.size(); i++) {
			int idx = round(after_start);
			int dist = pts[idx]->distance(*p);
			//		cout << p->get_header() << " " << pts[idx]->get_header() << " " << dist << endl;
			buf.push_back(make_pair(p, pts[idx]));
			after_start += after_inc;
		}
#pragma omp critical
		pairs.insert(std::end(pairs), std::begin(buf), std::end(buf));
//			cout << "added " << pairs.size() - size_before << " pairs" << endl;
	}
	return pairs;
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
	sf.emplace_back([](Point<T>* a, Point<T>* b) {
			return Feature<T>::length_difference(*a, *b);
		}, false);
	sf.emplace_back([](Point<T>* a, Point<T>* b) {
			return Feature<T>::intersection(*a, *b);
		}, true);
	sf.emplace_back([](Point<T>* a, Point<T>* b) {
			return Feature<T>::jenson_shannon(*a, *b);
		}, false);
	sf.emplace_back([](Point<T>* a, Point<T>* b) {
			return Feature<T>::simratio(*a, *b);
		}, true);
	sf.emplace_back([](Point<T>* a, Point<T>* b) {
			return Feature<T>::squaredchord(*a, *b);
		}, false);
	// sf.emplace_back([](Point<T>* a, Point<T>* b) {
	// 		return Feature<T>::manhattan(*a, *b);
	// 	}, false);
	// sf.emplace_back([](Point<T>* a, Point<T>* b) {
	// 		return Feature<T>::pearson(*a, *b);
	// 	}, true);
	return;
	int four_k = pow(4, k);
	vector<int> reverse;
	vector<int> reverse_complement;
	auto freverse = [](int idx, int k) {
		int sum = 0;
		for (int i = 0; i < k; i++) {
			int rem = idx % 4;
			idx /= 4;
			sum = 4 * sum + rem;

		}
		return sum;
	};
	auto freverse_complement = [](int idx, int k) {
		std::vector<int> v;
		for (int i = 0; i < k; i++) {
			v.push_back(3 - idx % 4);
			idx /= 4;
		}
		int sum = 0;
		for (auto val : v) {
			sum = 4 * sum + val;
		}
		return sum;
	};
	auto print_index = [](int idx, int k) {
		for (int i = 0; i < k; i++) {
			int num = idx % 4;
			idx /= 4;
			cout << num;
		}
		cout << endl;
	};
	for (int i = 0; i < four_k; i++) {
		reverse.push_back(freverse(i, k));
		reverse_complement.push_back(freverse_complement(i, k));
	}
	cout << "After computing: size: " << reverse.size() << " " << reverse_complement.size() << endl;
	sf.emplace_back([](Point<T>* a, Point<T> *b, const vector<int>& c, const vector<int> &d) {
			return Feature<T>::n2rrc(*a, *b, c, d);
				}, reverse, reverse_complement, true);

}
template class Trainer<uint8_t>;
template class Trainer<uint16_t>;
template class Trainer<uint32_t>;
template class Trainer<uint64_t>;
template class Trainer<int>;
template class Trainer<double>;
