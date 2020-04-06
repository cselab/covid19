#include "model.h"
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <cstdio>
#include <cstdlib>

using Stepper = boost::numeric::odeint::runge_kutta_dopri5<MultiSEIIRState>;

RegionalSEIIR::deterministicRHS(const Parameters &p, const MultiSEIIRState &x, MultiSEIIRState &dxdt) const {
    size_t numRegions = x.S.size();
    dxdt.S.resize(size);
    dxdt.E.resize(size);
    dxdt.Ir.resize(size);
    dxdt.Iu.resize(size);
    dxdt.N.resize(size);

    for (size_t i = 0; i < numRegions; ++i) {
        double A = p.beta * x.S[i] / x.N[i] * x.Ir[i];
        double B = p.beta * x.S[i] / x.N[i] * p.mu * x.Iu[i];
        double E_Z = x.E[i] / p.Z;
        dxdt.S[i] = -(A + B);
        dxdt.E[i] = A + B - E_Z;
        dxdt.Ir[i] = p.alpha * E_Z - x.Ir[i] / D;
        dxdt.Iu[i] = E_Z - p.alpha * E_Z - x.Iu[i] / D;
        dxdt.N[i] = 0;

        for (size_t j = 0; j < numRegions; ++j) {
            double Tij = this->M(i, j) / (x.N[j] - x.Ir[j]);
            double Tji = this->M(j, i) / (x.N[i] - x.Ir[i]);
            dxdt.S[i] += p.theta * (Tij * x.S[j] - Tji * x.S[i]);
            dxdt.E[i] += p.theta * (Tij * x.E[j] - Tji * x.E[i]);
            // Documented infected people are in quarantine, they do not move around.
            // dxdt.Ir[i] += p.theta * (Tij * x.Ir[j] - Tji * x.Ir[i]);
            dxdt.Iu[i] += p.theta * (Tij * x.Iu[j] - Tji * x.Iu[i]);
            dxdt.N[i] += p.theta * (this->M(i, j) - this->M(j, i));
        }
    }
}

RegionalSEIIR::RegionalSEIIR(std::vector<int> population, std::vector<int> commuteMatrix) :
    numRegions_(population.size()),
    N_(std::move(population)),
    M_(std::move(commuteMatrix))
{
    if (M_.size() != numRegions_ * numRegions_) {
        fprintf(stderr, "Matrix Mij number of elements (%zu) does not match number of regions^2 (%zu^2).\n",
                M_.size(), numRegions_);
        exit(1);
    }
}

std::vector<MultiSEIIRState> solve(Parameters parameters, int days) const {
    const int STEPS_PER_DAY = 10;
    const double dt = 1.0 / STEPS_PER_DAY;

    auto rhs = [this, parameters](const MultiSEIIRState &x, MultiSEIIRState &dxdt, double t) {
        deterministicRHS(parameters, x, dxdt, t);
    };

    std::vector<MultiSEIIRState> result;
	result.reserve(days);

    // Observer gets called for each time step evaluated by the integrator.
    // We consider only every `STEPS_PER_DAY` steps (skipping the first one as well).
    auto observer = [&result, N, cnt = 0](const MultiSEIIRState &y, double t) mutable {
        if (cnt % STEPS_PER_DAY == 0 && cnt > 0)
            result.push_back(y);
        ++cnt;
    };

    boost::numeric::odeint::integrate_n_steps(
            Stepper{}, rhs, initialState, 0.0, dt, days * STEPS_PER_DAY, observer);
    return result;
}
