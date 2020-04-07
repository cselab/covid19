#include <korali.hpp>
#include "model.h"
#include "data.h"
#include "utils.h"

#include <cstdlib>

// This does not really fit here...
std::vector<double> getStandardDeviationModel(
        const CantonsData &data,
        double sigma,
        const std::vector<State> &/*result*/) {
    return std::vector<double>(data.dataPoints.size(), sigma);
}

State createInitialState(const CantonsData &data) {
    State y0{data.numCantons};

    for (size_t i = 0; i < data.numCantons; ++i) {
        y0.S(i) = data.population[i];
        y0.E(i) = 0;
        y0.Ir(i) = 0;
        y0.Iu(i) = 0;
        y0.N(i) = data.population[i];
    }

    int ticino = data.nameToIndex.at("TI");
    y0.Ir(ticino) = 10;

    return y0;
}


auto makeKoraliWrapper(
        const CantonsData &data,
        const Solver &solver,
        const State &y0,
        int numDays) {
    return [&data, &solver, &y0, numDays](korali::Sample& k) {
        Parameters params{
            k["Parameters"][0],
            k["Parameters"][1],
            k["Parameters"][2],
            k["Parameters"][3],
            k["Parameters"][4],
            k["Parameters"][5],
        };
        printf("%g %g %g %g %g %g  %g\n",
                (double)k["Parameters"][0],
                (double)k["Parameters"][1],
                (double)k["Parameters"][2],
                (double)k["Parameters"][3],
                (double)k["Parameters"][4],
                (double)k["Parameters"][5],
                (double)k["Parameters"][6]);
        double sigma = k["Parameters"][6];
        auto result = solver.solve(params, y0, numDays);
        k["Reference Evaluations"] = data.getReferenceEvaluations(result);
        k["Standard Deviation Model"] = getStandardDeviationModel(data, sigma, result);
    };
}


int main() {
    korali::Engine k;
    korali::Experiment e;
    CantonsData data = readCantonsData();
    Solver solver{data.Mij};
    State y0 = createInitialState(data);
    int numDays = 50;

    e["Problem"]["Type"] = "Bayesian/Reference";
    e["Problem"]["Likelihood Model"] = "Additive Normal General";
    e["Problem"]["Reference Data"] = data.getReferenceData();
    e["Problem"]["Computational Model"] = makeKoraliWrapper(data, solver, y0, numDays);

    e["Solver"]["Type"] = "TMCMC";
    e["Solver"]["Population Size"] = 2000;
    // e["Solver"]["Termination Criteria"]["Max Generations"] = 30;

    e["Distributions"][0]["Name"] = "Prior for beta";
    e["Distributions"][0]["Type"] = "Univariate/Uniform";
    e["Distributions"][0]["Minimum"] = 0.8;
    e["Distributions"][0]["Maximum"] = 1.5;

    e["Distributions"][1]["Name"] = "Prior for mu";
    e["Distributions"][1]["Type"] = "Univariate/Uniform";
    e["Distributions"][1]["Minimum"] = 0.2;
    e["Distributions"][1]["Maximum"] = 1.0;

    e["Distributions"][2]["Name"] = "Prior for alpha";
    e["Distributions"][2]["Type"] = "Univariate/Uniform";
    e["Distributions"][2]["Minimum"] = 0.02;
    e["Distributions"][2]["Maximum"] = 1.0;

    e["Distributions"][3]["Name"] = "Prior for Z";
    e["Distributions"][3]["Type"] = "Univariate/Uniform";
    e["Distributions"][3]["Minimum"] = 2.0;  // Days.
    e["Distributions"][3]["Maximum"] = 5.0;  // Days.

    e["Distributions"][4]["Name"] = "Prior for D";
    e["Distributions"][4]["Type"] = "Univariate/Uniform";
    e["Distributions"][4]["Minimum"] = 2.0;  // Days.
    e["Distributions"][4]["Maximum"] = 15.0;  // Days.

    e["Distributions"][5]["Name"] = "Prior for theta";
    e["Distributions"][5]["Type"] = "Univariate/Uniform";
    e["Distributions"][5]["Minimum"] = 1.0;
    e["Distributions"][5]["Maximum"] = 1.75;

    e["Distributions"][6]["Name"] = "Prior for [Sigma]";
    e["Distributions"][6]["Type"] = "Univariate/Uniform";
    e["Distributions"][6]["Minimum"] = 0.0;
    e["Distributions"][6]["Maximum"] = +600.0;

    const char *names[] = {"beta", "mu", "alpha", "Z", "D", "theta", "[Sigma]"};
    for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); ++i) {
        e["Variables"][i]["Name"] = names[i];
        e["Variables"][i]["Prior Distribution"] = std::string("Prior for ") + names[i];
    }

    k["Conduit"]["Type"] = "Concurrent";
    k["Conduit"]["Concurrent Jobs"] = 8;

    k.run(e);

    return 0;
}
