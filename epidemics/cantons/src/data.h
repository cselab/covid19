#pragma once

#include <map>
#include <string>
#include <vector>

/// A data value for the given region and day.
struct DataPoint {
    int day;
    int region;
    double value;
};

/// Model bulk parameters (not optimized for).
struct ModelData {
    // Note: py/model.py:ModelData.to_cpp depends on this structure.
    size_t numRegions;
    std::map<std::string, int> regionNameToIndex;
    std::vector<int> regionPopulation;
    std::vector<double> Mij;  // Row-major [to][from].
    std::vector<double> Cij;  // Row-major [to][from].
    std::vector<double> externalCases;  // Row-major [day][canton].

    double getExternalCasesAt(int day, int canton) const noexcept {
        int idx = day * (int)numRegions + canton;
        return idx < 0 || idx >= (int)externalCases.size() ? 0 : externalCases[idx];
    }
};

/// UQ-specific data
struct ReferenceData {
    // Note: py/model.py:ReferenceData.to_cpp depends on this structure.

    /// List of known number of cases.
    std::vector<DataPoint> cases;

    /// Represent all known data point values into a single vector.
    /// Used for setting up Korali.
    std::vector<double> getReferenceData() const;

    /*
    /// Extract number of infected people for days and regions for which we
    /// have measured data, in the same order as the values returned by
    /// `getReferenceData()`.
    std::vector<double> getReferenceEvaluations(
            const std::vector<State> &states) const;
    */
};

ModelData readModelData(const char *filename = "data/cpp_model_data.dat");
ReferenceData readReferenceData(const char *filename = "data/cpp_reference_data.dat");
