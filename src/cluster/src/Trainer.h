/* -*- C++ -*- */
#ifndef TRAINER_H
#define TRAINER_H

#include "Point.h"
#include "GLM.h"
#include "Feature.h"
template<class T>
class Trainer {
public:
	Trainer(std::vector<Point<T>*> v, size_t num_points, double cutoff_, size_t max_pts_from_one_, double (&matrix)[4][4], double sig, double eps, int ksize) : points(v), n_points(num_points), cutoff(cutoff_), max_pts_from_one(max_pts_from_one_), k(ksize) {
		init(matrix, sig, eps);
	};
	~Trainer() { delete feat_mat; }
	std::pair<std::map<std::pair<Point<T>*, Point<T>*>, double>,
		  std::map<std::pair<Point<T>*, Point<T>*>, double> > split_old();
	vector<std::pair<Point<T>*,Point<T>*> > split();
	double train_n(pair<vector<pair<Point<T>*,
			 Point<T>*
			 > >,
	     vector<pair<Point<T>*,
		       Point<T>*> > > &data, int ncols);
	void train(double acc_cutoff=97.5);
	vector<pair<int, double> > get_close(Point<T>*, const vector<pair<Point<T>*,int> > &,  bool& is_min) const;
	void filter(Point<T>*, vector<pair<Point<T>*,bool> >&) const;
	Point<T>* closest(Point<double>*, vector<pair<Point<T>*,bool> >&) const;
	Point<T>* merge(Point<T>*, vector<pair<Point<T>*,double> >&) const;
private:
	matrix::GLM glm;
	matrix::Matrix weights;
	double align(Point<T>* a, Point<T>* b) const;
	void init(double (&matrix)[4][4], double sig, double eps);
	pair<vector<pair<Point<T>*,
			 Point<T>*
			 > >,
	     vector<pair<Point<T>*,
			 Point<T>*> > > get_labels(vector<std::pair<Point<T>*,Point<T>*> >&, double cutoff) const;
	vector<SingleFeature<T> > sf;
	vector<Feature<T> > ff;


	int mat[4][4];
	int sigma, epsilon;
	std::vector<Point<T>*> points;
	matrix::Matrix *feat_mat = NULL;
	size_t n_points, max_pts_from_one;
	double cutoff;
	int k;
};
#endif
