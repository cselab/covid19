/// Code common to country and cantons models.

#pragma once

#include "bindings.h"
#include <epidemics/integrator.hh>

namespace epidemics {

IntegratorSettings integratorSettingsFromKwargs(py::kwargs kwargs);

template <typename Solver, typename State, typename Parameters, typename PySolver>
void exportSolverSolve(
        PySolver &pySolver,
        const char *name)
{
    using namespace py::literals;
    pySolver.def(
            name,
            [](const Solver &solver,
               const Parameters &params,
               State state,
               const std::vector<double> &tEval,
               py::kwargs kwargs)
            {
                SignalRAII breakRAII;
                return solver.solve(params, std::move(state), tEval,
                                    integratorSettingsFromKwargs(kwargs));
            }, "parameters"_a, "initial_state"_a, "t_eval"_a);
}

}  // namespace epidemics
