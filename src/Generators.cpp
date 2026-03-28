#include "Generators.hpp"

Matrix GenerateRandomMatrix(int n, double left, double right, unsigned seed)
{   
    if (n <= 0)
        throw std::invalid_argument("Matrix size must be positive");
    
    if (left > right)
        throw std::invalid_argument("Invalid range");
    
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> dist(left, right);
    Matrix A(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A.SetElement(i, j, dist(generator));
        }
    }
    
    return A;
}

std::vector<double> GenerateRandomVector(int n, double left, double right, unsigned seed)
{
    if (n <= 0)
        throw std::invalid_argument("Vector size must be positive");
    if (left > right)
        throw std::invalid_argument("Invalid range");
    
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> dist(left, right);

    std::vector<double> V(n);
    for (int i = 0; i < n; i++) {
        V[i] = dist(generator);
    }

    return V;
}

Matrix GenerateHilbertMatrix(int n)
{
    if (n <= 0)
        throw std::invalid_argument("Matrix size must be positive");
    
    Matrix H(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            H(i, j) = 1.0 / (i + j + 1);
        }
    }

    return H;
}

