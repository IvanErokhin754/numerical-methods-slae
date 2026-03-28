#include <gtest/gtest.h>
#include "Matrix.hpp"

std::vector<double> BackSubstitution(const Matrix &U, const std::vector<double> &b);
std::vector<double> ForwardSubstitution(const Matrix &L, const std::vector<double> &b);
std::vector<double> GaussSolve(Matrix A, std::vector<double> b);
std::vector<double> GaussSolvePivot(Matrix A, std::vector<double> b);
std::vector<double> LUsolve(const Matrix &A, const std::vector<double> &b);
void LUDecomposition(const Matrix &A, Matrix &L, Matrix &U);

TEST(BackSubstitution, Simple2x2)
{
    Matrix U;
    U.SetSizeM(2);
    U.SetSizeN(2);

    std::vector<std::vector<double>> data = {
        {2.0, 1.0},
        {0.0, 3.0}
    };
    U.SetMatr(data);
    
    std::vector<double> b = {5.0, 6.0};

    std::vector<double> x = BackSubstitution(U, b);

    EXPECT_NEAR(x[0], 1.5, 1e-9);
    EXPECT_NEAR(x[1], 2.0, 1e-9);
}

TEST(ForwardSubstitution, Simple2x2)
{
    Matrix L;
    L.SetSizeN(2);
    L.SetSizeM(2);

    std::vector<std::vector<double>> data = {
        {2.0, 0.0},
        {3.0, 1.0}
    };
    L.SetMatr(data);

    std::vector<double> b = {4.0, 5.0};

    std::vector<double> y = ForwardSubstitution(L, b);

    EXPECT_NEAR(y[0], 2.0, 1e-9);
    EXPECT_NEAR(y[1], -1.0, 1e-9);
}

TEST(GaussSolve, Simple2x2)
{
    Matrix A;
    A.SetSizeN(2);
    A.SetSizeM(2);

    std::vector<std::vector<double>> data = {
        {2.0, 1.0},
        {1.0, 3.0}
    };
    A.SetMatr(data);

    std::vector<double> b = {5.0, 5.0};

    std::vector<double> x = GaussSolve(A, b);

    EXPECT_NEAR(x[0], 2.0, 1e-9);
    EXPECT_NEAR(x[1], 1.0, 1e-9);
}

TEST(GaussSolvePivot, Simple2x2)
{
    Matrix A;
    A.SetSizeN(2);
    A.SetSizeM(2);

    std::vector<std::vector<double>> data = {
        {2.0, 1.0},
        {1.0, 3.0}
    };
    A.SetMatr(data);

    std::vector<double> b = {5.0, 5.0};

    std::vector<double> x = GaussSolvePivot(A, b);

    EXPECT_NEAR(x[0], 2.0, 1e-9);
    EXPECT_NEAR(x[1], 1.0, 1e-9);
}

TEST(GaussSolvePivot, ZeroPivotButSolvable)
{
    Matrix A;
    A.SetSizeN(2);
    A.SetSizeM(2);

    std::vector<std::vector<double>> data = {
        {0.0, 1.0},
        {2.0, 3.0}
    };
    A.SetMatr(data);

    std::vector<double> b = {1.0, 5.0};

    std::vector<double> x = GaussSolvePivot(A, b);

    EXPECT_NEAR(x[0], 1.0, 1e-9);
    EXPECT_NEAR(x[1], 1.0, 1e-9);
}

TEST(LUsolve, Simple2x2)
{   
    Matrix A;
    A.SetSizeN(2);
    A.SetSizeM(2);

    std::vector<std::vector<double>> data = {
        {2.0, 1.0},
        {1.0, 3.0}
    };
    A.SetMatr(data);

    std::vector<double> b = {5.0, 5.0};
    std::vector<double> x = LUsolve(A, b);

    EXPECT_NEAR(x[0], 2.0, 1e-9);
    EXPECT_NEAR(x[1], 1.0, 1e-9);
}

TEST(LUDecomposition, CheckProduct)
{
    Matrix A;
    A.SetSizeN(3);
    A.SetSizeM(3);

    std::vector<std::vector<double>> data = {
        {2.0, 1.0, 1.0},
        {4.0, 3.0, 3.0},
        {8.0, 7.0, 9.0}
    };
    A.SetMatr(data);
    
    Matrix L, U;
    LUDecomposition(A, L, U);

    Matrix R = L.Mult(U);

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            EXPECT_NEAR(R.GetElement(i, j), A.GetElement(i, j), 1e-9);
        }
    }
}