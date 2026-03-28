#ifndef _GAUSS_H_
#define _GAUSS_H_
#include <vector>
#include "Matrix.hpp"

std::vector<double> GaussSolve(Matrix A, std::vector<double> b);

std::vector<double> GaussSolvePivot(Matrix A, std::vector<double> b);

std::vector<double> BackSubstitution(const Matrix& U, const std::vector<double> &b);

std::vector<double> ForwardSubstitution(const Matrix &L, const std::vector<double> &b);

std::vector<double> LUsolve(const Matrix &A, const std::vector<double> &b);

void LUDecomposition(const Matrix &A, Matrix& L, Matrix &U);

double ResidualNorm(const Matrix &A, const std::vector<double> &x, const std::vector<double> &b);


#endif /* _GAUSS_H_ */
