#include "model.h"
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using Stepper = boost::numeric::odeint::runge_kutta_dopri5<MultiSEIIRState>;

/*
 * The view classes for reading and writing the state vector in a human-readable way.
 */
struct MultiSEIIRStateView {
    double &S (size_t i) const { return p_[0 * numRegions_ + i]; }
    double &E (size_t i) const { return p_[1 * numRegions_ + i]; }
    double &Ir(size_t i) const { return p_[2 * numRegions_ + i]; }
    double &Iu(size_t i) const { return p_[3 * numRegions_ + i]; }
    double &N (size_t i) const { return p_[4 * numRegions_ + i]; }

    MultiSEIIRStateView(size_t numRegions, double *p) :
        numRegions_{numRegions}, p_{p}
    { }

private:
    size_t numRegions_;
    double *p_;
};

struct MultiSEIIRStateConstView {
    double S (size_t i) const { return p_[0 * numRegions_ + i]; }
    double E (size_t i) const { return p_[1 * numRegions_ + i]; }
    double Ir(size_t i) const { return p_[2 * numRegions_ + i]; }
    double Iu(size_t i) const { return p_[3 * numRegions_ + i]; }
    double N (size_t i) const { return p_[4 * numRegions_ + i]; }

    MultiSEIIRStateConstView(size_t numRegions, const double *p) :
        numRegions_{numRegions}, p_{p}
    {}

private:
    size_t numRegions_;
    const double *p_;
};


static std::vector<double> transposeMatrix(const std::vector<double> &m, size_t N) {
    std::vector<double> out(N * N, 0.0);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            out[i * N + j] = m[j * N + i];
    return out;
}

MultiSEIIR::MultiSEIIR(std::vector<double> commuteMatrix) :
    numRegions_{static_cast<size_t>(std::sqrt(commuteMatrix.size()) + 0.1)},
    M_{std::move(commuteMatrix)},
    Mt_{transposeMatrix(M_, numRegions_)}
{
    if (M_.size() != numRegions_ * numRegions_) {
        fprintf(stderr, "The number of elements of Mij is not a perfect square: %zu\n", M_.size());
        exit(1);
    }
}

std::vector<MultiSEIIRState> MultiSEIIR::solve(Parameters parameters, MultiSEIIRState initialState, int days) const {
    const int STEPS_PER_DAY = 10;
    const double dt = 1.0 / STEPS_PER_DAY;

    auto rhs = [this, parameters](const MultiSEIIRState &x, MultiSEIIRState &dxdt, double /*t*/) {
        deterministicRHS(parameters, x, dxdt);
    };

    std::vector<MultiSEIIRState> result;
	result.reserve(days);

    // Observer gets called for each time step evaluated by the integrator.
    // We consider only every `STEPS_PER_DAY` steps (skipping the first one as well).
    auto observer = [&result, cnt = 0](const MultiSEIIRState &y, double t) mutable {
        if (cnt % STEPS_PER_DAY == 0 && cnt > 0)
            result.push_back(y);
        ++cnt;
    };

    evalCounter_ = 0;
    boost::numeric::odeint::integrate_n_steps(
            Stepper{}, rhs, initialState, 0.0, dt, days * STEPS_PER_DAY, observer);
    fprintf(stderr, "eval counter = %d\n", evalCounter_);
    return result;
}

void MultiSEIIR::deterministicRHS(
        Parameters p,
        const MultiSEIIRState &x_,
        MultiSEIIRState &dxdt_) const {
    ++evalCounter_;
    dxdt_.resize(x_.size());

    MultiSEIIRStateConstView x{numRegions_, x_.data()};
    MultiSEIIRStateView dxdt{numRegions_, dxdt_.data()};

    for (size_t i = 0; i < numRegions_; ++i) {
        double A = p.beta * x.S(i) / x.N(i) * x.Ir(i);
        double B = p.beta * x.S(i) / x.N(i) * p.mu * x.Iu(i);
        double E_Z = x.E(i) / p.Z;

        double dS = -(A + B);
        double dE = A + B - E_Z;
        double dIr = p.alpha * E_Z - x.Ir(i) / p.D;
        double dIu = E_Z - p.alpha * E_Z - x.Iu(i) / p.D;
        double dN = 0;

        double inv = 1 / (x.N(i) - x.Ir(i));
        for (size_t j = 0; j < numRegions_; ++j) {
            double Tij = this->M (i, j) / (x.N(j) - x.Ir(j));
            double Tji = this->Mt(i, j) * inv;
            dS += p.theta * (Tij * x.S(j) - Tji * x.S(i));
            dE += p.theta * (Tij * x.E(j) - Tji * x.E(i));
            // Documented infected people are in quarantine, they do not move around.
            // dIr += p.theta * (Tij * x.Ir(j) - Tji * x.Ir(i));
            dIu += p.theta * (Tij * x.Iu(j) - Tji * x.Iu(i));
            dN += p.theta * (this->M(i, j) - this->Mt(i, j));
        }

        dxdt.S(i) = dS;
        dxdt.E(i) = dE;
        dxdt.Ir(i) = dIr;
        dxdt.Iu(i) = dIu;
        dxdt.N(i) = dN;
    }
}
