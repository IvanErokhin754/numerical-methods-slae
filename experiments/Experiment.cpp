#include "Experiments.hpp"
#include <fstream>
#include <limits>

void SaveResultsToCSV(const std::string& filename,
                      const std::vector<ExperimentResult>& results)
{
    std::ofstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Failed to open output file");

    file << "n,k,method,time,residual,relative_error\n";

    for (const auto& r : results) {
        file << r.n << ","
             << r.k << ","
             << r.method << ","
             << r.time << ","
             << r.residual << ","
             << r.relative_error << "\n";
    }
}

void PrintResults(const std::vector<ExperimentResult> &results)
{
    for (const auto &r : results) {
        std::cout << "n = " << r.n
                  << ", k = " << r.k
                  << ", method = " << r.method
                  << ", time = " << r.time
                  << ", residual = " << r.residual
                  << ", relative_error = " << r.relative_error
                  << "\n";
    }
}

std::vector<ExperimentResult> ExperimentSingleSystem()
{
    std::vector<ExperimentResult> results;
    std::vector<int> sizes = {100, 200, 500, 1000};

    using Clock = std::chrono::high_resolution_clock;

    for (int n : sizes) {
        Matrix A = GenerateRandomMatrix(n, -1.0, 1.0, 751);
        std::vector<double> b = GenerateRandomVector(n, -1.0, 1.0, 752);

        auto start = Clock::now();
        std::vector<double> x1 = GaussSolve(A, b);
        auto end = Clock::now();
        double t1 = std::chrono::duration<double>(end - start).count();
        results.push_back({n, 0, "Gauss", t1, ResidualNorm(A, x1, b), 0.0});

        start = Clock::now();
        std::vector<double> x2 = GaussSolvePivot(A, b);
        end = Clock::now();
        double t2 = std::chrono::duration<double>(end - start).count();
        results.push_back({n, 0, "GaussPivot", t2, ResidualNorm(A, x2, b), 0.0});

        start = Clock::now();
        std::vector<double> x3 = LUsolve(A, b);
        end = Clock::now();
        double t3 = std::chrono::duration<double>(end - start).count();
        results.push_back({n, 0, "LU Total", t3, ResidualNorm(A, x3, b), 0.0});

        Matrix L, U;
        start = Clock::now();
        LUDecomposition(A, L, U);
        end = Clock::now();
        double tLUDecomp = std::chrono::duration<double>(end - start).count();
        results.push_back({n, 0, "LU Decomp", tLUDecomp, 0.0, 0.0});

        start = Clock::now();
        std::vector<double> y = ForwardSubstitution(L, b);
        std::vector<double> x4 = BackSubstitution(U, y);
        end = Clock::now();
        double tLUSolve = std::chrono::duration<double>(end - start).count();
        results.push_back({n, 0, "LU Solve", tLUSolve, ResidualNorm(A, x4, b), 0.0});
    }

    return results;
}

std::vector<ExperimentResult> ExperimentHilbert()
{
    std::vector<ExperimentResult> results;
    std::vector<int> sizes = {5, 10, 15};

    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (int n : sizes) {
        Matrix H = GenerateHilbertMatrix(n);
        std::vector<double> xTrue(n, 1.0);
        std::vector<double> b = MultiplyMatrixVector(H, xTrue);

        try {
            std::vector<double> x1 = GaussSolve(H, b);
            results.push_back({
                n, 0, "Gauss", 0.0,
                ResidualNorm(H, x1, b),
                RelativeError(x1, xTrue)
            });
        }
        catch (const std::exception &e) {
            std::cout << "Hilbert n = " << n
                      << ", method = Gauss failed: "
                      << e.what() << "\n";
            results.push_back({n, 0, "Gauss", 0.0, nan, nan});
        }

        try {
            std::vector<double> x2 = GaussSolvePivot(H, b);
            results.push_back({
                n, 0, "GaussPivot", 0.0,
                ResidualNorm(H, x2, b),
                RelativeError(x2, xTrue)
            });
        }
        catch (const std::exception &e) {
            std::cout << "Hilbert n = " << n
                      << ", method = GaussPivot failed: "
                      << e.what() << "\n";
            results.push_back({n, 0, "GaussPivot", 0.0, nan, nan});
        }

        try {
            std::vector<double> x3 = LUsolve(H, b);
            results.push_back({
                n, 0, "LU", 0.0,
                ResidualNorm(H, x3, b),
                RelativeError(x3, xTrue)
            });
        }
        catch (const std::exception &e) {
            std::cout << "Hilbert n = " << n
                      << ", method = LU failed: "
                      << e.what() << "\n";
            results.push_back({n, 0, "LU", 0.0, nan, nan});
        }
    }

    return results;
}

std::vector<ExperimentResult> ExperimentMultipleRightParts()
{
    std::vector<ExperimentResult> results;
    std::vector<int> ks = {1, 10, 100};
    const int n = 500;

    using Clock = std::chrono::high_resolution_clock;

    Matrix A = GenerateRandomMatrix(n, -1.0, 1.0, 751);

    for (int k : ks) {
        std::vector<std::vector<double>> rights;
        rights.reserve(k);

        for (int i = 0; i < k; i++) {
            rights.push_back(GenerateRandomVector(n, -1.0, 1.0, 100 + i));
        }

        auto start = Clock::now();
        for (int i = 0; i < k; i++) {
            std::vector<double> x = GaussSolvePivot(A, rights[i]);
        }
        auto end = Clock::now();
        double tGauss = std::chrono::duration<double>(end - start).count();
        results.push_back({n, k, "GaussPivot Total", tGauss, 0.0, 0.0});

        start = Clock::now();
        Matrix L, U;
        LUDecomposition(A, L, U);
        for (int i = 0; i < k; i++) {
            std::vector<double> y = ForwardSubstitution(L, rights[i]);
            std::vector<double> x = BackSubstitution(U, y);
        }
        end = Clock::now();
        double tLU = std::chrono::duration<double>(end - start).count();
        results.push_back({n, k, "LU Total", tLU, 0.0, 0.0});
    }

    return results;
}