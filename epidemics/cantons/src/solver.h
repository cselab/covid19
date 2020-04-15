#pragma once

#include "model.h"

class Solver {
public:
    Solver(ModelData modelData, bool verbose = false);

    double M(int from, int to) const {
        return modelData_.Mij[from * modelData_.numRegions + to];
    }

    std::vector<State> solve(const Parameters &parameters, State initialState, int days) const;
    std::vector<State> solve(const Parameters &parameters, RawState initialState, int days) const;

private:
    void deterministicRHS(int day, Parameters p, const RawState &x, RawState &dxdt) const;

    ModelData modelData_;
    bool verbose_;
};

using CheckSignalsFunc = void(*)();

// Initially nullptr. If set, it will be called from the solver on every time step.
extern CheckSignalsFunc check_signals_func;
