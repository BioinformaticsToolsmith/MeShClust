#include "Feature.h"
#include "DivergencePoint.h"
#include <cmath>
#include <limits>
#include "../../utility/GlobAlignE.h"

template<class T>
void Feature<T>::add_feature(uint16_t f_flags, int combo)
{
	cout << "Adding combo " << f_flags << endl;
	if (combo != COMBO_SQUARED && combo != COMBO_SELF) {
		throw "invalid combo";
	}
	vector<int> indices;
	for (uint16_t f = 1; f <= f_flags; f = (f << 1)) {
		// it is in the new parameter but not currently in store
		if ((f_flags & f) != 0) {
			if ((flags & f) == 0) {
				lookup.push_back(f);
				cout << "new single feature " << f << endl;
				mins.push_back(std::numeric_limits<double>::max());
				maxs.push_back(std::numeric_limits<double>::min());
				is_sims.push_back(feat_is_sim(f));
				is_finalized.push_back(false);
				flags |= f;
			}
			indices.push_back(index_of(f));
		}
	}
	combos.push_back(std::make_pair(combo, indices));
}


template<class T>
void Feature<T>::finalize()
{
	for (size_t i = 0; i < is_finalized.size(); i++) {
		is_finalized[i] = true;
	}
}
template<class T>
void Feature<T>::normalize_cache(vector<double> &cache) const
{
	for (size_t i = 0; i < lookup.size(); i++) {
		double val = (cache[i] - mins[i]) / (maxs[i] - mins[i]);
		if (is_sims[i]) {
			cache[i] = val;
		} else {
			cache[i] = 1 - val;
		}
	}
}
template<class T>
vector<double> Feature<T>::compute_all_raw(Point<T> &p, Point<T> &q)
{
	vector<double> cache(lookup.size());
	uint16_t done = 0;
	// if (flags & (FEAT_INTERSECTION | FEAT_JENSONSHANNON | FEAT_SIMRATIO)) {
	// 	auto three = three_features(p, q);
	// 	size_t inter_i = index_of(FEAT_INTERSECTION);
	// 	size_t js_i = index_of(FEAT_JENSONSHANNON);
	// 	size_t sim_i = index_of(FEAT_SIMRATIO);
	// 	cache[inter_i] = get<0>(three);
	// 	cache[js_i] = get<1>(three);
	// 	cache[sim_i] = get<2>(three);
	// 	done |= (FEAT_INTERSECTION | FEAT_JENSONSHANNON | FEAT_SIMRATIO);
	// } else if (flags & (FEAT_INTERSECTION | FEAT_JENSONSHANNON)) {
	// 	auto two = two_features(p, q);
	// 	size_t inter_i = index_of(FEAT_INTERSECTION);
	// 	size_t js_i = index_of(FEAT_JENSONSHANNON);
	// 	cache[inter_i] = two.first;
	// 	cache[js_i] = two.second;
	// 	done |= (FEAT_INTERSECTION | FEAT_JENSONSHANNON);
	// }

#pragma omp parallel for
	for (size_t i = 0; i < lookup.size(); i++) {
		if ((lookup[i] & done) == 0) {
			auto rres = raw(lookup[i], p, q);
			cache[i] = rres;
		}
	}
	return cache;
}

template<class T>
void Feature<T>::normalize(const vector<pair<Point<T>*,Point<T>*> > &pairs)
{

	for (size_t i = 0; i < lookup.size(); i++) {
		double small = mins[i], big = maxs[i];
		if (lookup[i] == FEAT_ALIGN) {
			mins[i] = 0;
			maxs[i] = 1;
			continue;
		}
		if (is_finalized[i]) {
			continue;
		}
#pragma omp parallel for reduction(min:small), reduction(max:big)
		for (size_t j = 0; j < pairs.size(); j++) {
			double val = raw(lookup[i], *pairs[j].first, *pairs[j].second);
			if (val < small) {
				small = val;
			}
			if (val > big) {
				big = val;
			}
		}

		mins[i] = small;
		maxs[i] = big;
	}
};


template<class T>
double Feature<T>::raw(uint16_t single_flag, Point<T>& a, Point<T>& b)
{
	double val = 0;
	switch (single_flag) {
	case FEAT_ALIGN:
		val = align(a, b, atable);
		break;
	case FEAT_LD:
		val = length_difference(a, b);
		break;
	case FEAT_MANHATTAN:
		val = manhattan(a, b);
		break;
	case FEAT_SQCHORD:
		val = squaredchord(a, b);
		break;
	case FEAT_INTERSECTION:
		val = intersection(a, b);
		break;
	case FEAT_PEARSON:
		val = pearson(a, b);
		break;
	case FEAT_SIMRATIO:
		val = simratio(a, b);
		break;
	case FEAT_N2RRC:
		cerr << "n2rrc not implemented" << endl;
		throw single_flag;
	case FEAT_JENSONSHANNON:
		val = jenson_shannon(a, b);
		break;
	case FEAT_RREE_K_R:
		val = rree_k_r(a, b);
		break;
	case FEAT_KULCZYNSKI2:
		val = kulczynski2(a, b);
		break;
	default:
		cerr << "bad feature flag " << single_flag << endl;
		throw single_flag;
	}
	return val;
}
template<class T>
bool Feature<T>::feat_is_sim(uint16_t single_flag) const
{
	bool is_sim = true;
	switch (single_flag) {
	case FEAT_ALIGN:
		is_sim = true;
		break;
	case FEAT_LD:
		is_sim = false;
		break;
	case FEAT_MANHATTAN:
		is_sim = false;
		break;
	case FEAT_SQCHORD:
		is_sim = false;
		break;
	case FEAT_INTERSECTION:
		is_sim = true;
		break;
	case FEAT_PEARSON:
		is_sim = false;
		break;
	case FEAT_SIMRATIO:
		is_sim = true;
		break;
	case FEAT_N2RRC:
		cerr << "n2rrc not implemented" << endl;
		throw single_flag;
	case FEAT_JENSONSHANNON:
		is_sim = false;
		break;
	case FEAT_RREE_K_R:
		is_sim = false;
		break;
	case FEAT_KULCZYNSKI2:
		is_sim = true;
		break;
	default:
		cerr << "bad feature flag " << single_flag << endl;
		throw single_flag;
	}
	return is_sim;
}

template<class T>
double Feature<T>::kulczynski2(Point<T> &a, Point<T> &b)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	const int N = p.points.size();
	uint64_t min_sum = 0;
	double ap = (double)p.getPseudoMagnitude() / N;
	double aq = (double)q.getPseudoMagnitude() / N;
	for (int i = 0; i < N; i++) {
		min_sum += std::min(p.points[i], q.points[i]);
	}
	double coeff = N * (ap + aq) / (2 * ap * aq);
	return coeff * min_sum;
}
template<class T>
double Feature<T>::align(Point<T> &a, Point<T> &b, std::map<std::pair<uintmax_t, uintmax_t>, double> &atbl)
{
	auto ai = a.get_id();
	auto bi = b.get_id();
	std::pair<uintmax_t, uintmax_t> pr = ai < bi ? std::make_pair(ai, bi) : std::make_pair(bi, ai);
	auto res = atbl.find(pr);
	if (res == atbl.end()) {
		auto sa = a.get_data_str();
		auto sb = b.get_data_str();
		int la = sa.length();
		int lb = sb.length();
		GlobAlignE galign(sa.c_str(), 0, la-1,
				  sb.c_str(), 0, lb-1,
				  1, -1, 2, 1);
		double val = galign.getIdentity();
#pragma omp critical
		atbl[pr] = val;
		return val;
	} else {
		return res->second;
	}
}

template<class T>
double Feature<T>::squaredchord(Point<T> &a, Point<T> &b)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	const int N = p.points.size();
	double sum = 0;
	for (int i = 0; i < N; i++) {
		sum += p.points[i] + q.points[i] - 2 * sqrt(p.points[i] * q.points[i]);
	}
	return sum;
}

template<class T>
double Feature<T>::intersection(Point<T> &a, Point<T> &b)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	const int N = p.points.size();
	uintmax_t dist = 0;
	uintmax_t mag = p.getPseudoMagnitude() + q.getPseudoMagnitude();
	#pragma omp simd
	for (int i = 0; i < N; i++) {
		dist += 2 * std::min(p.points[i], q.points[i]);
	}
	return (double)dist / (double)mag;
}

template<class T>
double Feature<T>::pearson(Point<T> &a, Point<T> &b)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	const int N = p.points.size();
	double dap = (double)p.getPseudoMagnitude() / N;
	double daq = (double)q.getPseudoMagnitude() / N;
	int ap = round(dap);
	int aq = round(daq);
	int dot = 0, np = 0, nq = 0;
	for (int i = 0; i < N; i++) {
	        auto dp = p.points[i] - ap;
	        auto dq = q.points[i] - aq;
		np += dp * dp;
		nq += dq * dq;
		dot += dp * dq;
	}
	return dot / sqrt(np * nq);
}

template<class T>
double Feature<T>::simratio(Point<T> &a, Point<T> &b)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	const int N = p.points.size();
	uintmax_t dot = 0, norm2 = 0;
	for (int i = 0; i < N; i++) {
		intmax_t diff = p.points[i] - q.points[i];
		dot += p.points[i] * q.points[i];
		norm2 += diff * diff;
	}
	return dot / (dot + sqrt(norm2));
}
template<class T>
double Feature<T>::manhattan(Point<T>& a, Point<T>& b)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	int N = p.points.size();
	int sum = 0;
	#pragma omp simd
	for (int i = 0; i < N; i++) {
		sum += p.points[i] > q.points[i] ? p.points[i] - q.points[i] : q.points[i] - p.points[i];
	}
//	std::cout << "manhattan: " << sum << std::endl;
	return sum;
}

template<class T>
double Feature<T>::length_difference(Point<T>& a, Point<T>& b)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	auto lp = p.get_length();
	auto lq = q.get_length();
	if (lp == 0 || lq == 0) {
		cerr << "lp: " << lp << " lq: " << lq << endl;
		throw 123;
	}
	auto ret = (lp > lq) ? (lp - lq) : (lq - lp);
//	std::cout << "length difference: " << ret << std::endl;
	return ret;
}


double neighbor(double *cp, double *cq, double ap, double aq, const int N)
{
	double sp = 0, sq = 0;
	#pragma omp simd
	for (int i = 0; i < N; i++) {
		double dp = cp[i] - ap;
		double dq = cq[i] - aq;
		sp += dp * dp;
		sq += dq * dq;
	}
	sp = sqrt(sp / N);
	sq = sqrt(sq / N);
	double psum = 0, qsum = 0;
	#pragma omp simd
	for (int i = 0; i < N; i++) {
		cp[i] = (cp[i] - ap) / sp;
		cq[i] = (cq[i] - aq) / sq;
		psum += cp[i] * cp[i];
		qsum += cq[i] * cq[i];
	}
	double total = 0;
	psum = sqrt(psum);
	qsum = sqrt(qsum);
	#pragma omp simd
	for (int i = 0; i < N; i++) {
		cp[i] /= psum;
		cq[i] /= qsum;
		total += cp[i] * cq[i];
	}
	return total;
}
template<class T>
double Feature<T>::n2rrc(Point<T>& a, Point<T>& b, const vector<int>& reverse, const vector<int>& reverse_complement)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	const int N = p.points.size();
	double *cp = new double[N];
	double *cq = new double[N];
	double ap = 0, aq = 0;
	for (int i = 0; i < N; i++) {
		int j = reverse.at(i);
		int h = reverse_complement.at(i);
		cp[i] = p.points[h] + p.points[i] + p.points[j];
		cq[i] = q.points[h] + q.points[i] + q.points[j];
		ap += cp[i];
		aq += cq[i];
	}
	ap /= N;
	aq /= N;
	double total = neighbor(cp, cq, ap, aq, N);
	delete[] cp;
	delete[] cq;
//	std::cout << "n2rrc: " << total << std::endl;
	return total;
}

/*
 * found at
 * http://www.machinedlearnings.com/2011/06/fast-approximate-logarithm-exponential.html
 */
inline float fastlog2(float x)
{
	union { float f; uint32_t i; } vx = { x };
	union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | (0x7e << 23) };
	float y = vx.i;
	y *= 1.0 / (1 << 23);
	return y - 124.22544637f - 1.498030302f * mx.f - 1.72587999f / (0.3520887068f + mx.f);
}
inline float fastlog4(float x)
{
	return fastlog2(x) / 2;
}
inline float fastlog(float x)
{
	return 0.69314718f * fastlog2(x);
}
template<class T>
double Feature<T>::jenson_shannon(Point<T> &a, Point<T> &b) const
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	uint64_t mp = p.getPseudoMagnitude();
	uint64_t mq = q.getPseudoMagnitude();
	double sum = 0;
	const int N = p.points.size();
        #pragma omp simd reduction(+:sum)
	for (int i = 0; i < N; i++) {
		double pp = (double)p.points[i] / mp;
		double pq = (double)q.points[i] / mq;
		double avg = 0.5 * (pp + pq);
		#ifdef USELOG
		double lp = // tbl[(int)(coeff * pp / avg)];
			log(pp / avg);
		double lq = // tbl[(int)(coeff * pq / avg)];
			log(pq / avg);
		#else
		double lp = tbl[(int)(coeff * pp / avg)];
		double lq = tbl[(int)(coeff * pq / avg)];
		#endif
	        sum += pp * lp + pq * lq;
	}
	return sum / 2;
}

template<class T>
double Feature<T>::rree_k_r(Point<T>& a, Point<T>& b)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	const int N = p.points.size();
	double op = 0, oq = 0;
	for (size_t offset = 0; offset < N; offset += 4) {
		int psum = 0, qsum = 0;
		#pragma omp simd
		for (int j = 0; j < 4; j++) {
			psum += p.points[offset+j];
			qsum += q.points[offset+j];
		}
		double ip = 0, iq = 0;
		for (int j = 0; j < 4; j++) {
			double cp = (double)p.points[offset + j] / psum;
			double cq = (double)q.points[offset + j] / qsum;
			double avg = 0.5 * (cp + cq);
			ip += cp * fastlog4(cp / avg);
			iq += cq * fastlog4(cq / avg);
		}
		op += ip;
		oq += iq;
	}
        double val = 0.5 * (op + oq);
//	std::cout << "RREE: " << val << std::endl;
	return val;
}

template<class T>
std::pair<double,double> two_features(Point<T>& a, Point<T>& b)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	const int N = p.points.size();
	const uint64_t mp = p.getPseudoMagnitude();
	const uint64_t mq = q.getPseudoMagnitude();
	const uint64_t mag = mp + mq;
	uint64_t dist = 0;
	double js_sum = 0;
#pragma omp parallel for reduction(+:dist), reduction(+:js_sum)
	for (int i = 0; i < N; i++) {
		auto pi = p.points[i];
		auto qi = q.points[i];
		dist += std::min(pi, qi);
		double pp = (double)pi / mp;
		double pq = (double)qi / mq;
		double avg = 0.5 * (pp + pq);
		pp = pp * fastlog(pp / avg);
		pq = pq * fastlog(pq / avg);
		js_sum += pp + pq;
	}

	return std::make_pair((double)dist / mag,
			      js_sum / 2.0);
}

template<class T>
std::tuple<double,double,double> three_features(Point<T>& a, Point<T>& b)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	const int N = p.points.size();
	const uint64_t mp = p.getPseudoMagnitude();
	const uint64_t mq = q.getPseudoMagnitude();
	const uint64_t mag = mp + mq;
	uint64_t dist = 0;
	double js_sum = 0;
	int64_t dot = 0, norm2 = 0;
	for (int i = 0; i < N; i++) {
		auto pi = p.points[i];
		auto qi = q.points[i];
		dist += std::min(pi, qi);
		dot += pi * qi;
	        int64_t diff = pi - qi;
		norm2 += diff * diff;
		double pp = (double)pi / mp;
		double pq = (double)qi / mq;
		double avg = 0.5 * (pp + pq);
		pp = pp * fastlog(pp / avg);
		pq = pq * fastlog(pq / avg);
		js_sum += pp + pq;
	}
	dist *= 2;
	return std::make_tuple((double)dist / mag,
			       js_sum / 2.0,
			       dot / (dot + sqrt(norm2)));
}
template<class T>
std::tuple<double,double,double,double> all_features(Point<T>& a, Point<T>& b)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	const int N = p.points.size();
	const uint64_t mp = p.getPseudoMagnitude();
	const uint64_t mq = q.getPseudoMagnitude();
	const uint64_t mag = mp + mq;
	double sum = 0, dot = 0, norm2 = 0;
	uintmax_t dist = 0;
	double sqchord = 0;
	for (int i = 0; i < N; i++) {
		auto pi = p.points[i];
		auto qi = q.points[i];
		double diff = pi - qi;
		double prod = pi * qi;
		dist += std::min(pi, qi);
		dot += prod;
		sqchord += sqrt(prod);
		norm2 += diff * diff;
		double pp = pi / mp;
		double pq = qi / mp;
		double avg = 0.5 * (pp + pq);
		double lp = fastlog(pp / avg);
		double lq = fastlog(pq / avg);
		sum += pp * lp + pq * lq;
	}
	dist *= 2;

	sqchord *= -2;
	sqchord += mag;

	return make_tuple((double)dist / (double)mag,
			  sum / 2,
			  dot / (dot + sqrt(norm2)),
			  sqchord
		);
}
template class Feature<uint8_t>;
template class Feature<uint16_t>;
template class Feature<uint32_t>;
template class Feature<uint64_t>;
template class Feature<int>;
template class Feature<double>;
