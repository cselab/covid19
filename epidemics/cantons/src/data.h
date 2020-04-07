#pragma once

#include "model.h"

#include <map>
#include <string>
#include <vector>

/*
 * Single known data point. Note that we do not know the number of infected
 * people for each day for each canton.
 */
struct DataPoint {
    int day;
    int country;
    double value;
};

/*
 * Collection of all canton data.
 */
struct CantonsData {
    size_t numCantons;
    std::map<std::string, int> nameToIndex;
    std::vector<int> population;
    std::vector<double> Mij;  // Row-major.
    std::vector<DataPoint> dataPoints;

    /// Represent all known data point values into a single vector.
    /// Used for setting up Korali.
    std::vector<double> getReferenceData() const;

    /// Extract number of infected people for days and cantons for which we
    /// have measured data, in the same order as the values returned by
    /// `getReferenceData()`.
    std::vector<double> getReferenceEvaluations(
            const std::vector<State> &states) const;
};

CantonsData readCantonsData(const char *filename = "data/cantons_data.dat");
