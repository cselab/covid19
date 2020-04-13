#include "model.h"
#include "utils.h"

std::vector<double> ValidationData::getReferenceData() const {
    size_t N = cases.size();
    std::vector<double> out(N, 0.0);
    for (size_t i = 0; i < N; ++i)
        out[i] = cases[i].value;
    return out;
}

std::vector<double> ValidationData::getReferenceEvaluations(
        const std::vector<State> &states) const {
    size_t N = cases.size();
    std::vector<double> out(N, 0.0);
    for (size_t i = 0; i < N; ++i) {
        out[i] = states[cases[i].day].Ir(cases[i].region);
    }
    return out;
}

ModelData readModelData(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (f == nullptr)
        DIE("Error opening file \"%s\". Did you forget to run ./py/data.py?\n", filename);

    int N;  // Number of regions.
    if (fscanf(f, "%d", &N) != 1)
        DIE("Reading number of regions failed.");

    ModelData out;
    out.numRegions = N;
    for (int i = 0; i < N; ++i) {
        char name[16];
        if (fscanf(f, "%s", name) != 1)
            DIE("Reading name of the region #%d failed.\n", i);
        out.regionNameToIndex[name] = i;
    }

    out.regionPopulation.resize(N);
    for (int &pop : out.regionPopulation)
        if (fscanf(f, "%d", &pop) != 1)
            DIE("Reading region population failed.\n");

    out.Mij.resize(N * N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (fscanf(f, "%lf", &out.Mij[i * N + j]) != 1)
                DIE("Reading Mij[%d][%d] failed.\n", i, j);
            if (i == j)
                out.Mij[i * N + j] = 0.0;
        }

    int numDays;
    if (fscanf(f, "%d", &numDays) != 1)
        DIE("Reading numDays for external cases failed.\n");
    out.externalCases.resize(N * numDays);
    for (int i = 0; i < numDays; ++i)
        for (int j = 0; j < N; ++j)
            if (fscanf(f, "%lf", &out.externalCases[i * N + j]) != 1)
                DIE("Reading externalCases[day=%d][canton=%d] failed.\n", i, j);

    fclose(f);
}

ValidationData readValidationData(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (f == nullptr)
        DIE("Error opening file \"%s\". Did you forget to run ./py/data.py?\n", filename);

    ValidationData out;
    int M;
    if (fscanf(f, "%d", &M) != 1)
        DIE("Failed reading M.\n");
    out.cases.resize(M);
    for (int i = 0; i < M; ++i) {
        DataPoint &dp = out.cases[i];
        if (fscanf(f, "%d%d%lf", &dp.day, &dp.region, &dp.value) != 3)
            DIE("Failed reading number of cases #%d/%d.\n", i, M);
    }

    fclose(f);

    return out;
}
