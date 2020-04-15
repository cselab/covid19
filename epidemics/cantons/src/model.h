#pragma once

#include <cassert>
#include <map>
#include <string>
#include <vector>

/*
 * Using custom types with boost::odeint is not that simple.
 * Instead, we use a single std::vector<double> to store the
 * whole state of the simulation, and the wrapper below to
 * access the data in a human-friendly way.
 */
using RawState = std::vector<double>;

/// Human-friendly wrapper around RawState.
struct State {
    static constexpr int kVarsPerRegion = 5;

    /// Create a state with all values set to 0.
    explicit State(size_t numRegions) :
        numRegions_{numRegions},
        v_(kVarsPerRegion * numRegions, 0.0)
    { }
    explicit State(RawState state) :
        numRegions_{state.size() / kVarsPerRegion},
        v_{std::move(state)}
    {
        assert(v_.size() % kVarsPerRegion == 0);
    }
    double &S (size_t i) { return v_[0 * numRegions_ + i]; }
    double &E (size_t i) { return v_[1 * numRegions_ + i]; }
    double &Ir(size_t i) { return v_[2 * numRegions_ + i]; }
    double &Iu(size_t i) { return v_[3 * numRegions_ + i]; }
    double &N (size_t i) { return v_[4 * numRegions_ + i]; }

    double S (size_t i) const { return v_[0 * numRegions_ + i]; }
    double E (size_t i) const { return v_[1 * numRegions_ + i]; }
    double Ir(size_t i) const { return v_[2 * numRegions_ + i]; }
    double Iu(size_t i) const { return v_[3 * numRegions_ + i]; }
    double N (size_t i) const { return v_[4 * numRegions_ + i]; }

    size_t numRegions() const noexcept { return numRegions_; }
    const RawState &raw() const & { return v_; }
    RawState raw() && { return std::move(v_); }

private:
    size_t numRegions_;
    RawState v_;
};


/// Lightweight parameters (optimized for).
struct Parameters {
    double beta;   // Transmission rate.
    double mu;     // Reduction factor for transmission rate of undocumented individuals.
    double alpha;  // Fraction of documented infections.
    double Z;      // Average latency period.
    double D;      // Average duration of infection.
    double theta;  // Corrective multiplicative factor for Mij.
};

/// Model bulk parameters (not optimized for).
struct ModelData {
    // Note: py/model.py:ModelData.to_cpp depends on this structure.
    size_t numRegions;
    std::map<std::string, int> regionNameToIndex;
    std::vector<int> regionPopulation;
    std::vector<double> Mij;  // Row-major.
    std::vector<double> externalCases;  // Row-major [day][canton].

    double getExternalCasesAt(int day, int canton) const noexcept {
        int idx = day * (int)numRegions + canton;
        return idx < 0 || idx >= (int)externalCases.size() ? 0 : externalCases[idx];
    }
};

/// A data value for the given region and day.
struct DataPoint {
    int day;
    int region;
    double value;
};

/// UQ-specific data
struct ReferenceData {
    // Note: py/model.py:ReferenceData.to_cpp depends on this structure.

    /// List of known number of cases.
    std::vector<DataPoint> cases;

    /// Represent all known data point values into a single vector.
    /// Used for setting up Korali.
    std::vector<double> getReferenceData() const;

    /// Extract number of infected people for days and regions for which we
    /// have measured data, in the same order as the values returned by
    /// `getReferenceData()`.
    std::vector<double> getReferenceEvaluations(
            const std::vector<State> &states) const;
};

ModelData readModelData(const char *filename = "data/cpp_model_data.dat");
ReferenceData readReferenceData(const char *filename = "data/cpp_reference_data.dat");
