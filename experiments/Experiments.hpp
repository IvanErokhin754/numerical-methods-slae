#ifndef EXPERIMENTS_HPP
#define EXPERIMENTS_HPP

#include "Gauss.hpp"
#include "Matrix.hpp"
#include "Generators.hpp"
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

struct ExperimentResult {
    int n;
    int k;
    std::string method;
    double time;
    double residual;
    double relative_error;
};

void SaveResultsToCSV(const std::string& filename, const std::vector<ExperimentResult>& results);

std::vector<ExperimentResult> ExperimentSingleSystem();

std::vector<ExperimentResult> ExperimentMultipleRightParts();

std::vector<ExperimentResult> ExperimentHilbert();

void PrintResults(const std::vector<ExperimentResult> &results);

#endif /* EXPERIMENTS_HPP */
