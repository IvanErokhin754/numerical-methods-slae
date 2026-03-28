#include <vector>
#include <cmath>
#include <stdexcept>
#include "Gauss.hpp"

const double eps = 1e-12;

std::vector<double> ForwardSubstitution(const Matrix &L, const std::vector<double> &b)
{
    int n = L.GetSizeN();

    if (L.GetSizeN() != L.GetSizeM())
        throw std::invalid_argument("Matrix must be square!");
    if ((int)b.size() != n)
        throw std::invalid_argument("Wrong right-hand side size!");
    
    std::vector<double> y(n);

    for (int i = 0; i < n; i++){
        double sum = 0.0;

        for (int j = 0; j < i; j++){
            sum += L.GetElement(i, j) * y[j];
        }

        double diag = L.GetElement(i, i);

        if (std::fabs(diag) < eps)
            throw std::runtime_error("Zero diagonal during forward subtituion");

        y[i] = (b[i] - sum) / diag;
    }
    
    return y;
}

std::vector<double> BackSubstitution(const Matrix& U, const std::vector<double> &b)
{
    
    int n = U.GetSizeN();

    if (U.GetSizeN() != U.GetSizeM())
        throw std::invalid_argument("Matrix must be square");
    if ((int)b.size() != n)
        throw std::invalid_argument("Wrong right-hand side size");

    std::vector<double> x(n);



    for (int i = n - 1; i >= 0; i--){
        double sum = 0.0;

        for (int j = i + 1; j < n; j++){
            sum += U.GetElement(i, j) * x[j];
        }
        double diag = U.GetElement(i, i);

        if (std::fabs(diag) < eps)
            throw std::runtime_error("Zero diagonal during back substitution");

        x[i] = (b[i] - sum) / diag;
    }
    return x;
}

std::vector<double> GaussSolve(Matrix A, std::vector<double> b)
{
    int n = A.GetSizeN();

    if (A.GetSizeN() != A.GetSizeM())
        throw std::invalid_argument("Matrix must be square!");
    if ((int)b.size() != n)
        throw std::invalid_argument("Wrong right-hand side size!");

    for (int k = 0; k < n - 1; k++){
        double pivot = A.GetElement(k, k);

        if (std::fabs(pivot) < eps)
            throw std::runtime_error("Zero pivot in Gauss Method!");

        for (int i = k + 1; i < n; i++){
            double factor = A.GetElement(i, k) / pivot;

            A.AddLinearCombination(i, k, -factor);
            b[i] -= factor * b[k];
        }
    }

    return BackSubstitution(A, b);
}

std::vector<double> GaussSolvePivot(Matrix A, std::vector<double> b)
{
    int n = A.GetSizeN();

    if (A.GetSizeN() != A.GetSizeM())
        throw std::invalid_argument("Matrix must be square!");
    if ((int)b.size() != n)
        throw std::invalid_argument("Wrong right-hand side size!");

    for (int k = 0; k < n - 1; k++){
        int pivot_row = k;

        for (int i = k + 1; i < n; i++){
            if (std::fabs(A.GetElement(i, k)) > std::fabs(A.GetElement(pivot_row, k))){
                pivot_row = i;
            }
        }

        if (pivot_row != k){
            A.SwapRows(k, pivot_row);
            std::swap(b[k], b[pivot_row]);
        }
        
        double pivot = A.GetElement(k, k);

        if (std::fabs(pivot) < eps)
            throw std::runtime_error("Zero pivot in Gauss Method!");

        for (int i = k + 1; i < n; i++){
            double factor = A.GetElement(i, k) / pivot;

            A.AddLinearCombination(i, k, -factor);
            b[i] -= factor * b[k];
        }
    }
    
    return BackSubstitution(A, b);
}

void LUDecomposition(const Matrix &A, Matrix& L, Matrix &U)
{
    int n = A.GetSizeN();

    if (A.GetSizeN() != A.GetSizeM())
        throw std::invalid_argument("Matrix must be square!");
    
    L = Matrix(n, n);
    U = A;

    for (int i = 0; i < n; i++) {
        L(i, i) = 1.0;
    }

    for (int k = 0; k < n - 1; k++){
        double pivot = U.GetElement(k, k);

        if (std::fabs(pivot) < eps)
            throw std::runtime_error("Zero pivot in LU decomposition!");
        
        for (int i = k + 1; i < n; i++){
            double factor = U.GetElement(i, k) / pivot;

            L.SetElement(i, k, factor);
            U.AddLinearCombination(i, k, -factor);
        }
    }
}

std::vector<double> LUsolve(const Matrix &A, const std::vector<double> &b)
{
    Matrix L, U;
    LUDecomposition(A, L, U);

    std::vector<double> y = ForwardSubstitution(L, b);
    std::vector<double> x = BackSubstitution(U, y);

    return x;
}

double ResidualNorm(const Matrix &A, const std::vector<double> &x, const std::vector<double> &b) 
{
    std::vector<double> Ax = MultiplyMatrixVector(A, x);
    
    if (Ax.size() != b.size()) 
        throw std::invalid_argument("Vectors sizes do not match");  
    
    std::vector<double> r(Ax.size());

    for (int i = 0; i < Ax.size(); i++) {
        r[i] = Ax[i] - b[i];
    }

    return Norm(r);
}