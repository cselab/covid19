#include "data.h"
#include "utils.h"

std::vector<double> CantonsData::getReferenceData() const {
    size_t N = dataPoints.size();
    std::vector<double> out(N, 0.0);
    for (size_t i = 0; i < N; ++i)
        out[i] = dataPoints[i].value;
    return out;
}

std::vector<double> CantonsData::getReferenceEvaluations(
        const std::vector<State> &states) const {
    size_t N = dataPoints.size();
    std::vector<double> out(N, 0.0);
    for (size_t i = 0; i < N; ++i) {
        out[i] = states[dataPoints[i].day].Ir(dataPoints[i].country);
    }
    return out;
}

void fixMij(int N, std::vector<double> &v) {
    // Clear the diagonal.
    for (int i = 0; i < N; ++i)
        v[i * N + i] = 0.0;

    // Symmetrize.
    // The data probably says how many people in the region i work in region j.
    // Since every person goes back and forth, we add i->j and j->i.
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < i; ++j) {
            v[i * N + j] = v[j * N + i] = 1.0 * (v[i * N + j] + v[j * N + i]);
        }
}

CantonsData readCantonsData(const char *filename) {
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
    fixMij(N, out.Mij);

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

