/* -*- C++ -*- */
#ifndef TRAINER_H
#define TRAINER_H

#include "Point.h"
#include "GLM.h"
#include "Feature.h"
#include "bvec.h"
#include "Center.h"
#include "LogTable.h"
#include <set>
template<class T>
class Trainer {
public:
	Trainer(std::vector<Point<T>*> v, size_t num_points, int largest_count, double cutoff_, size_t max_pts_from_one_, double (&matrix)[4][4], double sig, double eps, int ksize) : points(v), n_points(num_points), max_pts_from_one(max_pts_from_one_), cutoff(cutoff_), k(ksize) {
		init(matrix, sig, eps);
		uintmax_t size = 1000 * 1000 * 10;
		log_table = new double[size];
		log_coeff = size / 2;
		double lsize = log(size);
		log_table[0] = 0;
		for (uintmax_t i = 1; i < size; i++) {
			log_table[i] = log(2 * i) - lsize;
		}
		feat = new Feature<T>(largest_count, log_table, log_coeff);
	};
	~Trainer() { delete feat_mat; delete feat; delete[] log_table;}
	std::pair<std::map<std::pair<Point<T>*, Point<T>*>, double>,
		  std::map<std::pair<Point<T>*, Point<T>*>, double> > split_old();
        vector<std::pair<Point<T>*,Point<T>*> > split();
	double train_n(pair<vector<pair<Point<T>*,
			 Point<T>*
			 > >,
	     vector<pair<Point<T>*,
		       Point<T>*> > > &data, int ncols);
	void train(double acc_cutoff=97.5);
	std::tuple<Point<T>*,double,size_t,size_t> get_close(Point<T>*, bvec_iterator<T> istart, bvec_iterator<T> iend,  bool& is_min) const;
//	vector<pair<int, double> > get_close(Point<T>*, const vector<pair<Point<T>*,int> > &,  bool& is_min) const;
	void filter(Point<T>*, vector<pair<Point<T>*,bool> >&) const;
	Point<T>* closest(Point<double>*, vector<pair<Point<T>*,bool> >&) const;
	long merge(vector<Center<T> > &centers, long current, long begin, long end) const;
//	Point<T>* merge(Point<T>*, vector<pair<Point<T>*,double> >&) const;
	double raw_classify(Point<T>*,Point<T>*) const;
private:
	matrix::GLM glm;
	matrix::Matrix weights;
	double align(Point<T>* a, Point<T>* b) const;
	std::pair<matrix::Matrix,matrix::Matrix> generate_feat_mat(pair<vector<pair<Point<T>*,
			 Point<T>*
			 > >,
	     vector<pair<Point<T>*,
					 Point<T>*> > > &data, int ncols);
	void init(double (&matrix)[4][4], double sig, double eps);


	pair<vector<pair<pair<Point<T>*,Point<T>*>, double> >,
     vector<pair<pair<Point<T>*,Point<T>*>, double > > > get_labels(vector<std::pair<Point<T>*,Point<T>*> >&, double cutoff) const;
	Feature<T> *feat;
	double *log_table;
	int mat[4][4];
	int sigma, epsilon;
	std::vector<Point<T>*> points;
	matrix::Matrix *feat_mat = NULL;
	size_t n_points, max_pts_from_one;
	double cutoff, log_coeff;
	int k;
	LogTable *tbl = NULL;
};
#endif
