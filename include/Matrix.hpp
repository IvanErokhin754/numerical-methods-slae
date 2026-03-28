#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <vector>

class Matrix{
public:
    Matrix();
    Matrix(int n, int m);
    Matrix(const Matrix &other);
    Matrix &operator=(const Matrix &other);
    double &operator()(int i, int j);
    const double &operator()(int i, int j) const;
    void SetSizeN(int size);
    void SetSizeM(int size);
    int GetSizeN() const;
    int GetSizeM() const;
    void SetElement(int i, int j, double value);
    double GetElement(int i, int j) const;
    Matrix Sum(const Matrix &other);
    Matrix Minus(const Matrix &other);
    Matrix Mult(const Matrix &other);
    void AddLinearCombination(int target, int source, double alpha);
    void SwapRows(int row_1, int row_2);
    void Transp();

private:
    int size_n, size_m;
    std::vector<std::vector<double>> A;

};

std::vector<double> MultiplyMatrixVector(const Matrix &A, const std::vector<double> &vector);

double Norm(const std::vector<double> &V);

double RelativeError(const std::vector<double> &x, const std::vector<double> &x_true);

#endif /* MATRIX_HPP */

