#include "parallel.h"
#include "sequential.h"
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>

// Grenzwerte für Benchmarkläufe
static const int    n_min                = 64;
static const int    n_max                = 1024;
static const int    alpha_min            = 1;
static const int    alpha_max            = 2;
static const int    z1_min               = 4;
static const int    z1_max               = 16;
static const int    z2_min               = 4;
static const int    z2_max               = 16;
static const int    h_max_factor_min     = 8;
static const int    h_max_factor_max     = 32;
static const double def_change_threshold = 1.0 / 10000;
static const int    def_max_iterations   = 100000;
static const int    seed                 = 0;

// Eingabefunktion
double f(double x, double y) {
    return 32 * (x * (1 - x) + y * (1 - y));
}

// Analytische Lösung
double u(double x, double y) {
    return 16 * x * (1 - x) * y * (1 - y);
};

// Erzeugt einen "Vector", wie für Jakobi, usw. benötigt.
// Der Vector wird mit der übergebenen Funktion f initialisiert.
template <typename Func>
vector_t createVector(int n, Func f) {
    const double h = 1.0 / (n + 1);
    vector_t u(n + 2, std::vector<double>(n + 2, 0.0));
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= n; ++j)
            u[i][j] = f(i * h, j * h);
    return u;
}

double computeMeanError(const vector_t &v1, const vector_t &v2) {
    double sum = 0.0;
    for (std::size_t i = 1; i < v1.size() - 1; ++i) {
        for (std::size_t j = 1; j < v1.size() - 1; ++j) {
            sum += std::abs(v1[i][j] - v2[i][j]);
        }
    }
    return sum / (v1.size() - 2) / (v1.size() - 2);
}

double computeMaximumError(const vector_t &v1, const vector_t &v2) {
    double max = 0.0;
    for (std::size_t i = 1; i < v1.size() - 1; ++i) {
        for (std::size_t j = 1; j < v1.size() - 1; ++j) {
            max = std::max(max, std::abs(v1[i][j] - v2[i][j]));
        }
    }
    return max;
}

bool operator==(const vector_t &v1, const vector_t &v2) {
    for (std::size_t i = 1; i < v1.size() - 1; ++i)
        for (std::size_t j = 1; j < v1.size() - 1; ++j)
            if (v1[i][j] != v2[i][j])
                return false;
    return true;
}

template <typename SeqFunc, typename ParFunc>
void executeBenchmark(int n, SeqFunc seqFunc, ParFunc parFunc) {
    using namespace std::chrono;
    std::default_random_engine rand(seed);
    std::uniform_real_distribution<double> dist(-10, 10); // Eingeschränkte Zufälligkeit
    auto startVector = createVector(n, [&](double, double) {return dist(rand);});
    auto anaResult = createVector(n, u);
    double h = 1.0 / (n + 1);

    auto time = high_resolution_clock::now();
    vector_t seqResult = seqFunc(startVector, h);
    auto seqTime = high_resolution_clock::now() - time;

    double seqMeanError = computeMeanError(anaResult, seqResult);
    double seqMaxError = computeMaximumError(anaResult, seqResult);
    double seqSeconds = duration_cast<duration<double>>(seqTime).count();
    std::cout << "Sequentiell: " << seqSeconds
              << " sek, Mittlerer Fehler: " << seqMeanError
              << ", Maximaler Fehler: " << seqMaxError << "\n";

    time = std::chrono::high_resolution_clock::now();
    vector_t parResult = parFunc(startVector, h);
    auto parTime = std::chrono::high_resolution_clock::now() - time;

    double parMeanError = computeMeanError(seqResult, parResult);
    double parSeconds = duration_cast<duration<double>>(parTime).count();
    std::cout << "Parallel: " << parSeconds
              << " sek, Speedup: " << seqSeconds / parSeconds
              << ", Mittlerer Fehler zu Seq: " << parMeanError <<  "\n";

    //if (seqResult != parResult)
    //    std::cout << "Debug: Die Ergebnisvektoren sind nicht äquivalent!\n";
}

void jakobiBenchmark() {
    std::cout << "Starte Jakobi Benchmark\n";

    int iter_count_seq, iter_count_par;
    auto seqFunc = [&](const vector_t &u, const double h) {
        return jakobi(u, f, h, iter_count_seq, def_change_threshold, def_max_iterations); };
    auto parFunc = [&](const vector_t &u, const double h) {
        return jakobiParallel(u, f, h, iter_count_par, def_change_threshold, def_max_iterations); };

    for (int n = n_min; n <= n_max; n*=2) {
        std::cout << "n=" << n << "\n";
        executeBenchmark(n, seqFunc, parFunc);
        std::cout << "Iterationen: Seq=" << iter_count_seq << ", Par=" << iter_count_par << "\n\n";
    }
}

void gaussSeidelBenchmark() {
    std::cout << "Starte Gauss-Seidel Benchmark\n";

    int iter_count_seq, iter_count_par;
    auto seqFunc = [&](const vector_t &u, const double h) {
        return gaussSeidel(u, f, h, iter_count_seq, def_change_threshold, def_max_iterations); };
    auto parFunc = [&](const vector_t &u, const double h) {
        return gaussSeidelParallel(u, f, h, iter_count_par, def_change_threshold, def_max_iterations); };

    for (int n = n_min; n <= n_max; n*=2) {
        std::cout << "n=" << n << "\n";
        executeBenchmark(n, seqFunc, parFunc);
        std::cout << "Iterationen: Seq=" << iter_count_seq << ", Par=" << iter_count_par << "\n\n";
    }
}

void mehrgitterBenchmark() {
    std::cout << "Starte Mehrgitter Benchmark\n";

    // Helfer: Liste schick darstellen.
    auto _s = [](const std::vector<std::pair<std::string, int>> &l) {
        std::string res = "[";
        if (l.size() > 0)
            res += l[0].first + ":" + std::to_string(l[0].second);
        for (std::size_t i = 1; i < l.size(); ++i)
            res += ", " + l[i].first + ":" + std::to_string(l[i].second);
        return res + "]";
    };

    for (int n = n_min; n <= n_max; n*=2) {
        for (int alpha = alpha_min; alpha <= alpha_max; alpha++) {
            for (int z1 = z1_min; z1 <= z1_max; z1*=2) {
                for (int z2 = z2_min; z2 <= z2_max; z2*=2) {
                    for (int h_max_factor = h_max_factor_min; h_max_factor <= h_max_factor_max; h_max_factor*=2) {
                        std::vector<std::pair<std::string, int>> iter_count_seq, iter_count_par;
                        auto seqFunc = [&](const vector_t &u, const double h) {
                            return mehrgitter(
                                u, f, z1, z2, h, h * h_max_factor, alpha,
                                iter_count_seq, def_change_threshold, def_max_iterations);
                        };
                        auto parFunc = [&](const vector_t &u, const double h) {
                            return mehrgitterParallel(
                                u, f, z1, z2, h, h * h_max_factor, alpha,
                                iter_count_par, def_change_threshold, def_max_iterations);
                        };
                        std::cout << "n=" << n << ", alpha=" << alpha << ", z1=" << z1 << ", z2=" << z2 << ", h_max_factor=" << h_max_factor << "\n";
                        executeBenchmark(n, seqFunc, parFunc);
                        std::cout << "Iterationen: " << _s(iter_count_seq) << "\n\n";
                    }
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    std::cout << std::fixed << std::setprecision(4);

    jakobiBenchmark();
    gaussSeidelBenchmark();
    mehrgitterBenchmark();
}
