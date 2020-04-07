#include <korali.hpp>
#include "model.h"

#include <cstdlib>
#include <map>
#include <string>
#include <vector>

void _die [[gnu::format(printf, 3, 4)]] (const char *filename, int line, const char *fmt, ...) {
    fprintf(stderr, "%s:%d ", filename, line);

    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);

    exit(1);
}

#define DIE(...) _die(__FILE__, __LINE__, __VA_ARGS__);

struct DataPoint {
    int day;
    int country;
    double value;
};

struct CantonsData {
    int numCantons;
    std::map<std::string, int> nameToIndex;
    std::vector<int> population;
    std::vector<double> Mij;  // Row-major.
    std::vector<DataPoint> dataPoints;

    std::vector<double> getReferenceData() const {
        size_t N = dataPoints.size();
        std::vector<double> out(N, 0.0);
        for (size_t i = 0; i < N; ++i)
            out[i] = dataPoints[i].value;
        return out;
    }
    std::vector<double> getReferenceEvaluations(
            const std::vector<std::vector<double>> &result) const {
        size_t N = dataPoints.size();
        std::vector<double> out(N, 0.0);
        for (size_t i = 0; i < N; ++i)
            out[i] = result[dataPoints[i].day][dataPoints[i].country];
        return out;
    }
};

// This does not really fit here...
std::vector<double> getStandardDeviationModel(
        const CantonsData &data,
        double sigma,
        const std::vector<std::vector<double>> &/*result*/) {
    size_t N = data.dataPoints.size();
    std::vector<double> out(N, sigma);
    return out;
}

CantonsData readCantonsData() {
    const char *filename = "data/cantons_data.dat";
    FILE *f = fopen(filename, "r");
    if (f == nullptr)
        DIE("Error opening file \"%s\". Aborting.\n", filename);

    int N;
    if (fscanf(f, "%d", &N) != 1)
        DIE("Reading N failed.");

    CantonsData out;
    out.numCantons = N;
    for (int i = 0; i < N; ++i) {
        char name[16];
        if (fscanf(f, "%s", name) != 1)
            DIE("Reading name of the canton #%d failed.\n", i);
        out.nameToIndex[name] = i;
    }

    out.population.resize(N);
    for (int &pop : out.population)
        if (fscanf(f, "%d", &pop) != 1)
            DIE("Reading population failed.\n");

    out.Mij.resize(N * N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (fscanf(f, "%lf", &out.Mij[i * N + j]) != 1)
                DIE("Reading Mij[%d][%d] failed.\n", i, j);
            if (i == j)
                out.Mij[i * N + j] = 0.0;
        }

    int M;
    if (fscanf(f, "%d", &M) != 1)
        DIE("Failed reading M.\n");
    out.dataPoints.resize(M);
    for (int i = 0; i < M; ++i) {
        DataPoint &dp = out.dataPoints[i];
        int tmp;
        if ((tmp = fscanf(f, "%d%d%lf", &dp.day, &dp.country, &dp.value)) != 3)
            DIE("Failed reading data point #%d/%d. tmp=%d\n", i, M, tmp);
    }

    fclose(f);

    return out;
}

std::vector<double> processResult(const CantonsData &data, const std::vector<std::vector<double>> &result) {
    std::vector<double> out;
    out.reserve(data.dataPoints.size());
    for (DataPoint dp : data.dataPoints)
        out.push_back(result[dp.day][dp.country]);
    return out;
}

State createInitialState(const CantonsData &data) {
    auto y0 = createEmptyState(data.numCantons);

    MultiRegionStateView view{y0};
    for (int i = 0; i < data.numCantons; ++i) {
        view.S(i) = data.population[i];
        view.E(i) = 0;
        view.Ir(i) = 0;
        view.Iu(i) = 0;
        view.N(i) = data.population[i];
    }

    int ticino = data.nameToIndex.at("TI");
    view.Ir(ticino) = 10;

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
    e["Solver"]["Population Size"] = 200;
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
