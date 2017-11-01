#include "Feature.h"
#include "DivergencePoint.h"
#include <cmath>

template<class T>
double Feature<T>::operator()(Point<T> *a, Point<T> *b) const
{
	vector<double> values;
//	std::cout << "call:" << endl;
	for (auto f : features) {
		double val = f(a, b);
		//std::cout << "value: " << val << endl;
		values.push_back(val);
	}
        double ret = combo(values);
//	std::cout << "Time: " << diff << std::endl;
	return ret;
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
	uintmax_t mag = 0;
	for (int i = 0; i < N; i++) {
		dist += 2 * std::min(p.points[i], q.points[i]);
		mag += p.points[i] + q.points[i];
	}
	return (double)dist / (double)mag;
}

template<class T>
double Feature<T>::pearson(Point<T> &a, Point<T> &b)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	const int N = p.points.size();
	double ap = p.getPseudoMagnitude() / N;
	double aq = q.getPseudoMagnitude() / N;
	double dot = 0, np = 0, nq = 0;
	for (int i = 0; i < N; i++) {
		double dp = p.points[i] - ap;
		double dq = q.points[i] - aq;
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
	double dot = 0, norm2 = 0;
	for (int i = 0; i < N; i++) {
		double diff = p.points[i] - q.points[i];
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
		sum += p.points[i] > q.points[i] ?
			p.points[i] - q.points[i] :
			q.points[i] - p.points[i];
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
double Feature<T>::jenson_shannon(Point<T> &a, Point<T> &b)
{
	const DivergencePoint<T>& p = dynamic_cast<const DivergencePoint<T>&>(a);
	const DivergencePoint<T>& q = dynamic_cast<const DivergencePoint<T>&>(b);
	double mp = 0, mq = 0;
	double sum = 0;
	const int N = p.points.size();
	for (int i = 0; i < N; i++) {
		mp += p.points[i];
		mq += q.points[i];
	}
	for (int i = 0; i < N; i++) {
		double pp = p.points[i] / mp;
		double pq = q.points[i] / mq;
		double avg = 0.5 * (pp + pq);
		double lp = fastlog(pp / avg);
		double lq = fastlog(pq / avg);
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
template class Feature<uint8_t>;
template class Feature<uint16_t>;
template class Feature<uint32_t>;
template class Feature<uint64_t>;
template class Feature<int>;
template class Feature<double>;
