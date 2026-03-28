#ifndef GENERATORS_HPP
#define GENERATORS_HPP

#include <random>
#include <vector>
#include "Matrix.hpp"

Matrix GenerateRandomMatrix(int n, double left, double right, unsigned seed);

std::vector<double> GenerateRandomVector(int n, double left, double right, unsigned seed);

Matrix GenerateHilbertMatrix(int n);

#endif /* GENERATORS_HPP */
