#include "../experiments/Experiments.hpp"

int main()
{
    auto r1 = ExperimentSingleSystem();
    PrintResults(r1);
    SaveResultsToCSV("single_system.csv", r1);

    auto r2 = ExperimentMultipleRightParts();
    PrintResults(r2);
    SaveResultsToCSV("multiple_rhs.csv", r2);

    auto r3 = ExperimentHilbert();
    PrintResults(r3);
    SaveResultsToCSV("hilbert.csv", r3);
    return 0;
}

