#include "sequential.h"
#include "parallel.h"
#include <cmath>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <functional>

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
    for (std::size_t i=0; i<v1.size(); i++) {
        for (std::size_t j=0; j<v1.size(); j++) {
            sum += std::abs(v1[i][j] - v2[i][j]);
        }
    }
    return sum / v1.size() / v1.size();
}

double computeMaximumError(vector_t v1, vector_t v2) {
    double max = 0.0;
    for (std::size_t i=0; i<v1.size(); i++) {
        for (std::size_t j=0; j<v1.size(); j++) {
            max = std::max(max, std::abs(v1[i][j] - v2[i][j]));
        }
    }
    return max;
}

bool vectorEquals(vector_t v1, vector_t v2) {
    for (std::size_t i=0; i<v1.size(); i++) {
        for (std::size_t j=0; j<v1.size(); j++) {
            if (v1[i][j] != v2[i][j])
                return false;
        }
    }
    return true;
}

void executeBenchmark(int n, std::function<vector_t (vector_t, double)> seqFunc,
                      std::function<vector_t (vector_t, double)> parFunc) {
    auto startVector = createVector(n, [](double, double) {return 1.0;}); // Immer noch nicht zufällig
    auto anaResult = createVector(n, u);
    double h = 1.0 / (n + 1);

    double time = getWallTime();
    vector_t seqResult = seqFunc(startVector, h);
    double seqTime = getWallTime() - time;

    double seqMedError = computeMedianError(anaResult, seqResult);
    double seqMaxError = computeMaximumError(anaResult, seqResult);
    std::cout << "Sequentiell: " << seqTime
              << "sek, Mittlerer Fehler: " << seqMedError
              << ", Maximaler Fehler: " << seqMaxError << "\n";

    time = getWallTime();
    vector_t parResult = parFunc(startVector, h);
    double parTime = getWallTime() - time;

    if (!vectorEquals(seqResult, parResult))
        std::cout << "Sequentielles und paralleles Ergebnis stimmen nicht überein!\n";
    std::cout << "Parallel: " << parTime
              << "sek, Speedup: " << seqTime / parTime << "\n\n";
}

void jakobiBenchmark(int nMin, int nMax) {
    std::cout << "Starte Jakobi Benchmark\n";

    for (int n = nMin; n <= nMax; n*=2) {
        auto seqFunc = [](vector_t v, double h){ return jakobi(v, f, h, 0.00001); };
        auto parFunc = [](vector_t v, double h){ return jakobiParallel(v, f, h, 0.00001); };
        std::cout << "n=" << n << "\n";
        executeBenchmark(n, seqFunc, parFunc);
    }
}

void gaussSeidelBenchmark(int nMin, int nMax) {
    std::cout << "Starte Gauss-Seidel Benchmark\n";

    for (int n = nMin; n <= nMax; n*=2) {
        auto seqFunc = [](vector_t v, double h){ return gaussSeidel(v, f, h, 0.00001, 1000000); };
        auto parFunc = [](vector_t v, double h){ return gaussSeidelParallel(v, f, h, 0.00001, 1000000); };
        std::cout << "n=" << n << "\n";
        executeBenchmark(n, seqFunc, parFunc);
    }
}

void mehrgitterBenchmark(int nMin, int nMax) {
    std::cout << "Starte Mehrgitter Benchmark\n";

    for (int n = nMin; n <= nMax; n*=2) {
        for (int alpha = 1; alpha <= 2; alpha++) {
            for (int z1 = 8; z1 <= 256; z1*=2) {
                for (int z2 = 8; z2 <= 256; z2*=2) {
                    auto seqFunc = [=](vector_t v, double h){
                        return mehrgitter(v, f, z1, z2, h, 4*h, alpha); };
                    auto parFunc = [=](vector_t v, double h){
                        return mehrgitterParallel(v, f, z1, z2, h, 4*h, alpha); };
                    std::cout << "n=" << n << ", alpha=" << alpha << ", z1=" << z1 << ", z2=" << z2 << "\n";
                    executeBenchmark(n, seqFunc, parFunc);
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    std::cout << std::fixed << std::setprecision(4);

    jakobiBenchmark(8, 128);
    gaussSeidelBenchmark(8, 128);
    mehrgitterBenchmark(8, 128);
}
