/* -*- C++ -*-
 *
 * Mat.h
 *
 * Author: Benjamin T James
 */
#ifndef MAT_H
#define MAT_H
#include <iostream>
#include <functional>
using namespace std;
template<class T>
class Mat {
public:
	Mat(function<T(int,int)> func, const long size) : n(size), table_size(size*(size+1)/2), compute(func) {
		if (size <= 0) {
			throw "Invalid size";
		}
		table = new T[table_size];
		set = new bool[table_size]();
	};
	~Mat() {
		delete[] table;
		delete[] set;
	};
	void fill() {
		unsigned long long count = 0;
		#ifdef OPENMP
                #pragma omp parallel for collapse(2) shared(set)
		#endif
		for (long i = 0; i < n; i++) {
			for (long j = 0; j < n; j++) {
				const auto idx = addr(i, j);
				if (!set[idx]) {
					auto res = compute(i, j);
					table[idx] = res;
					set[idx] = true;
					count++;
				}
				if (count % 10000 == 0) {
					cout << count << " / " << table_size << endl;
				}
			}
		}

	};
	T& operator[](pair<int,int> index) {
		const unsigned long idx = addr(index.first, index.second);
		if (!set[idx]) {
			table[idx] = compute(index.first, index.second);
			set[idx] = true;
		}
		return table[idx];
	};
	bool exists(int i, int j) const {
		return set[addr(i, j)];
	}
private:
	T* table;
	bool* set;
	const unsigned long table_size;
	const unsigned long n;
	function<T(int,int)> compute;

	unsigned long addr(unsigned long i, unsigned long j) const {
		if (i <= j) {
			return i * n - (i - 1) * i / 2 + j - i;
		} else {
			return j * n - (j - 1) * j / 2 + i - j;
		}
	};
};
#endif
