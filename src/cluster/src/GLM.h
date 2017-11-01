/*
 * glm.h
 *
 * Created on: May 29, 2017
 * Author: Robert Geraghty, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 * Modified by Benjamin James
 */

#ifndef SRC_MATRIX_GLM_H_
#define SRC_MATRIX_GLM_H_

#include "Matrix.h"
#include <tuple>
namespace matrix {

class GLM {
private:
	Matrix weights;

public:
	void train(matrix::Matrix& features, matrix::Matrix& labels);
	Matrix predict(matrix::Matrix& features) const;
	std::tuple<double,double,double> accuracy(matrix::Matrix& oLabels, matrix::Matrix& pLabels);
	const Matrix& get_weights() const { return weights; };
};

}

#endif /* SRC_MATRIX_GLM_H_ */
