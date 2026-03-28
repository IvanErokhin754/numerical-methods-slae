#include "Matrix.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

const double eps = 1e-10;

double &Matrix::operator()(int i, int j) {
    if (i >= size_n || j >= size_m || i < 0 || j < 0)
        throw std::out_of_range("Index out of range");
    return A[i][j];
}

const double &Matrix::operator()(int i, int j) const {
    if (i >= size_n || j >= size_m || i < 0 || j < 0)
        throw std::out_of_range("Index out of range");
    return A[i][j];
}

Matrix &Matrix::operator=(const Matrix &other) {
    if (this == &other)
        return *this;
    
    size_n = other.size_n;
    size_m = other.size_m;
    A = other.A;

    return *this;
}

Matrix::Matrix() : size_n(0), size_m(0) {}
Matrix::Matrix(int n, int m) : size_n(n), size_m(m), A(n, std::vector<double>(m, 0.0)) {}
Matrix::Matrix(const Matrix &other) : size_n(other.size_n), size_m(other.size_m), A(other.A) {}

void Matrix::SetSizeN(int size) {
    if (size <= 0)
        throw std::invalid_argument("Matrix size must be positive");
    size_n = size;
}
void Matrix::SetSizeM(int size) {
    if (size <= 0)
        throw std::invalid_argument("Matrix size must be positive");
    size_m = size;
}

int Matrix::GetSizeN() const{
    return size_n;
}
int Matrix::GetSizeM() const{
    return size_m;
}

void Matrix::SetElement(int i, int j, double value) {
    if (i >= size_n || j >= size_m || i < 0 || j < 0) 
        throw std::out_of_range("Index out of range");
    A[i][j] = value;
}

double Matrix::GetElement(int i, int j) const {
    if (i >= size_n || j >= size_m || i < 0 || j < 0) 
        throw std::out_of_range("Index out of range");
    return A[i][j];
}

Matrix Matrix::Sum(const Matrix &other){
    if (size_n != other.size_n || size_m != other.size_m)
        throw std::invalid_argument("Matrix sizes do not match");

    Matrix Res(size_n, size_m);
    
    for (int i = 0; i < size_n; i++){
        for (int j = 0; j < size_m; j++){
            Res.A[i][j] = A[i][j] + other.A[i][j];
        }
    }
    return Res;
}
Matrix Matrix::Minus(const Matrix &other){
    if (size_n != other.size_n || size_m != other.size_m)
        throw std::invalid_argument("Matrix sizes do not match");

    Matrix Res(size_n, size_m);
    
    for (int i = 0; i < size_n; i++){
        for (int j = 0; j < size_m; j++){
            Res.A[i][j] = A[i][j] - other.A[i][j];
        }
    }
    return Res;
}
Matrix Matrix::Mult(const Matrix &other){
    if (size_m != other.GetSizeN())
        throw std::invalid_argument("Matrix sizes do not match");

    Matrix Res(size_n, other.GetSizeM());
    
    for (int i = 0; i < size_n; i++){
        for (int j = 0; j < other.GetSizeM(); j++){
            double sum = 0.0;
            for (int k = 0; k < size_m; k++){
                sum += A[i][k] * other.A[k][j];
            }
            Res.A[i][j] = sum;
        }
    }

    return Res;
}
void Matrix::AddLinearCombination(int target, int source, double alpha){
    if (target < 0 || target >= size_n || source < 0 || source >= size_n)
        throw std::out_of_range("Row index out of range");

    for (int j = 0; j < size_m; j++){
        A[target][j] += alpha * A[source][j];
    }
}

void Matrix::SwapRows(int row_1, int row_2){
    if (row_1 < 0 || row_1 >= size_n || row_2 < 0 || row_2 >= size_n)
        throw std::out_of_range("Row index out of range");

    if (row_1 == row_2) return;
    std::swap(A[row_1], A[row_2]);
}

void Matrix::Transp(){
    if (size_n != size_m) 
        throw std::invalid_argument("Matrix is not square");
    for (int i = 0; i < size_n; i++){
        for (int j = i + 1; j < size_m; j++){
            std::swap(A[i][j], A[j][i]);
        }
    }
}

double Norm(const std::vector<double> &V) {
    double res = 0.0;
    
    for (int i = 0; i < V.size(); i++) {
        res += V[i] * V[i];
    }

    return std::sqrt(res);
}

std::vector<double> MultiplyMatrixVector(const Matrix &A, const std::vector<double> &vector) {
    int n = A.GetSizeN();
    int m = A.GetSizeM();

    if ((int)vector.size() != m)
        throw std::invalid_argument("Matrix and vector sizes do not match");
    
    std::vector<double> result(n, 0.0);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += A.GetElement(i, j) * vector[j];
        }
    }

    return result;
}

double RelativeError(const std::vector<double> &x, const std::vector<double> &x_true) {
    if (x.size() != x_true.size())
        throw std::invalid_argument("Vectors sizes do not match");

    std::vector<double> relative_vector(x.size());
    double relative_error = 0.0;

    for (int i = 0; i < x.size(); i++){
        relative_vector[i] = x[i] - x_true[i];
    }
    double norm_x_true = Norm(x_true);
    if (norm_x_true == 0.0)
        throw std::runtime_error("Exact solution norm is zero");

    relative_error = Norm(relative_vector) / norm_x_true;
    return relative_error;
}