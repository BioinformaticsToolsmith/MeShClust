/* -*- C++ -*-
 *
 * NearestNeighbor.h
 *
 * Author: Benjamin T James
 */
#ifndef NEARESTNEIGHBOR_H
#define NEARESTNEIGHBOR_H
// #include <ANN/ANN.h>
// #include "Point.h"
// template<class T>
// class NearestNeighbor {
// public:
// 	NearestNeighbor(const vector<Point<T>*> &pts) : points(pts) {
// 		const int dim = pts[0]->get_data().size();
// 		const int maxPts = pts.size();
// 		dataPts = annAllocPts(maxPts, dim);
// 		queryPt = annAllocPt(dim);
// 		for (int nPts = 0; nPts < maxPts; nPts++) {
// 			auto vec = pts[nPts]->get_data();
// 			for (int i = 0; i < vec.size(); i++) {
// 				dataPts[nPts][i] = vec[i];
// 			}
// 		}
// 		kd_tree = new ANNkd_tree(dataPts, maxPts, dim);
// 		nnIdx = new ANNidx[1];
// 		dists = new ANNdist[1];
// 	};
// 	~NearestNeighbor() {
// 		delete[] nnIdx;
// 		delete[] dists;
// 		delete kd_tree;
// 		annClose();
// 	};
// 	void find_nearest_neighbor(Point<T> &center) const {
// 		auto vec = center.get_data();
// 		for (int i = 0; i < vec.size(); i++) {
// 			queryPt[i] = vec[i];
// 		}
// 		kd_tree->annkSearch(queryPt, 1, nnIdx, dists);
// 		ANNidx idx = nnIdx[0];
// 		center.set(*points[idx]);
// 	};
// private:
// 	ANNkd_tree *kd_tree = NULL;
// 	ANNpointArray dataPts;
// 	ANNpoint queryPt;
// 	ANNidxArray nnIdx;
// 	ANNdistArray dists;
// 	const vector<Point<T>*> &points;
// };
#endif
