#include "sequential.h"
#include <cmath>
#include <sys/time.h>
#include <iostream>

double f(double x, double y) {
    return 32 * (x * (1 - x) + y * (1 - y));
}

double u(double x, double y) {
    return 16 * x * (1 - x) * y * (1 - y);
};

double getWallTime() {
    struct timeval time;
    gettimeofday(&time, NULL);
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

template <typename Func>
vector_t createVector(int n, Func f) {
    const double h = 1.0 / (n + 1);
    vector_t u(n + 2, std::vector<double>(n + 2, 0.0));
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= n; ++j)
            u[i][j] = f(i * h, j * h);
    return u;
}

double computeMedianError(vector_t v1, vector_t v2) {
    double sum = 0.0;
    for (int i=0; i<v1.size(); i++) {
        for (int j=0; j<v1.size(); j++) {
            sum += std::abs(v1[i][j] - v2[i][j]);
        }
    }
    return sum / v1.size() / v1.size();
}

double computeMaximumError(vector_t v1, vector_t v2) {
    double max = 0.0;
    for (int i=0; i<v1.size(); i++) {
        for (int j=0; j<v1.size(); j++) {
            max = std::max(max, std::abs(v1[i][j] - v2[i][j]));
        }
    }
    return max;
}

void jakobi(int nMin, int nMax) {
    std::cout << "Starte Jakobi Benchmark\n";

    for (int n = nMin; n <= nMax; n*=2) {
        auto startVector = createVector(n, [](double, double) {return 1.0;}); // Immer noch nicht zufällig
        auto analyticalResult = createVector(n, u);
        std::cout << "n=" << n << "\n";
        double h = 1.0 / (n + 1);

        double time = getWallTime();
        auto seqResult = jakobi(startVector, f, h, 0.00001);
        double seqTime = getWallTime() - time;
        double seqMedError = computeMedianError(anaResult, seqResult);
        double seqMaxError = computeMaximumError(anaResult, seqResult);
        std::cout << "Sequentiell: " << seqTime << "sek, Mittlerer Fehler: " << seqMedError << ", Maximaler Fehler: " << seqMaxError << "\n";

        //TODO parallel

        std::cout << "\n";
    }
}
