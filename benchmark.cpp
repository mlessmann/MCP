#include "sequential.h"
#include "parallel.h"
#include <cmath>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <functional>
#include <stdexcept>

// Standardwerte für Abbruchkriterien
static const double def_change_threshold = 1.0 / 1000;
static const int    def_max_iterations   = 1000;

// Eingabefunktion
double f(double x, double y) {
    return 32 * (x * (1 - x) + y * (1 - y));
}

// Analytische Lösung
double u(double x, double y) {
    return 16 * x * (1 - x) * y * (1 - y);
};

double getWallTime() {
    struct timeval time;
    gettimeofday(&time, NULL);
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

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

// FK: Das Ding ist nicht das, was ich unter Median verstehe. TODO: Nochmal checken.
double computeMedianError(const vector_t &v1, const vector_t &v2) {
    double sum = 0.0;
    for (std::size_t i=0; i<v1.size(); i++) {
        for (std::size_t j=0; j<v1.size(); j++) {
            sum += std::abs(v1[i][j] - v2[i][j]);
        }
    }
    return sum / v1.size() / v1.size();
}

double computeMaximumError(const vector_t &v1, const vector_t &v2) {
    double max = 0.0;
    for (std::size_t i=0; i<v1.size(); i++) {
        for (std::size_t j=0; j<v1.size(); j++) {
            max = std::max(max, std::abs(v1[i][j] - v2[i][j]));
        }
    }
    return max;
}

bool vectorEquals(const vector_t &v1, const vector_t &v2) {
    for (std::size_t i=0; i<v1.size(); i++) {
        for (std::size_t j=0; j<v1.size(); j++) {
            if (v1[i][j] != v2[i][j])
                return false;
        }
    }
    return true;
}

template <typename SeqFunc, typename ParFunc>
void executeBenchmark(int n, SeqFunc seqFunc, ParFunc parFunc) {
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
        // TODO: Toleranz einbauen? Gauß-Seidel Parallel ist nicht exakt identisch mit GS Sequenziell.
        //throw std::logic_error("Sequentielles und paralleles Ergebnis stimmen nicht überein!\n");
        std::cout << "Sequentielles und paralleles Ergebnis stimmen nicht überein!\n";
    std::cout << "Parallel: " << parTime
              << "sek, Speedup: " << seqTime / parTime << "\n";
}

void jakobiBenchmark(int nMin, int nMax) {
    std::cout << "Starte Jakobi Benchmark\n";

    int iter_count_seq, iter_count_par;
    auto seqFunc = [&](const vector_t &u, const double h) {
        return jakobi(u, f, h, iter_count_seq, def_change_threshold, def_max_iterations); };
    auto parFunc = [&](const vector_t &u, const double h) {
        return jakobiParallel(u, f, h, iter_count_par, def_change_threshold, def_max_iterations); };

    for (int n = nMin; n <= nMax; n*=2) {
        std::cout << "n=" << n << "\n";
        executeBenchmark(n, seqFunc, parFunc);
        std::cout << "Iterationen: Seq=" << iter_count_seq << ", Par=" << iter_count_par << "\n\n";
    }
}

void gaussSeidelBenchmark(int nMin, int nMax) {
    std::cout << "Starte Gauss-Seidel Benchmark\n";

    int iter_count_seq, iter_count_par;
    auto seqFunc = [&](const vector_t &u, const double h) {
        return gaussSeidel(u, f, h, iter_count_seq, def_change_threshold, def_max_iterations); };
    auto parFunc = [&](const vector_t &u, const double h) {
        return gaussSeidelParallel(u, f, h, iter_count_par, def_change_threshold, def_max_iterations); };

    for (int n = nMin; n <= nMax; n*=2) {
        std::cout << "n=" << n << "\n";
        executeBenchmark(n, seqFunc, parFunc);
        std::cout << "Iterationen: Seq=" << iter_count_seq << ", Par=" << iter_count_par << "\n\n";
    }
}

void mehrgitterBenchmark(int nMin, int nMax) {
    std::cout << "Starte Mehrgitter Benchmark\n";

    // Helfer: Liste schick darstellen.
    auto _s = [](const std::vector<int> &l) {
        std::string res = "[";
        if (l.size() > 0)
            res += std::to_string(l[0]);
        for (std::size_t i = 1; i < l.size(); ++i)
            res += ", " + std::to_string(l[i]);
        return res + "]";
    };

    for (int n = nMin; n <= nMax; n*=2) {
        for (int alpha = 1; alpha <= 2; alpha++) {
            for (int z1 = 8; z1 <= 256; z1*=2) {
                for (int z2 = 8; z2 <= 256; z2*=2) {
                    std::vector<int> iter_count_seq, iter_count_par;
                    auto seqFunc = [&](const vector_t &u, const double h) {
                        return mehrgitter(u, f, z1, z2, h, 4*h, alpha, iter_count_seq, def_change_threshold, def_max_iterations); };
                    auto parFunc = [&](const vector_t &u, const double h) {
                        return mehrgitterParallel(u, f, z1, z2, h, 4*h, alpha, iter_count_par, def_change_threshold, def_max_iterations); };
                    std::cout << "n=" << n << ", alpha=" << alpha << ", z1=" << z1 << ", z2=" << z2 << "\n";
                    executeBenchmark(n, seqFunc, parFunc);
                    std::cout << "Iterationen: Seq=" << _s(iter_count_seq) << ", Par=" << _s(iter_count_par) << "\n\n";
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    std::cout << std::fixed << std::setprecision(4);

    jakobiBenchmark(64, 512);
    gaussSeidelBenchmark(64, 512);
    //mehrgitterBenchmark(128, 256); // Zu viele Ausgaben, daher temporär auskommentiert. (TODO)
}
