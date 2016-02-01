#pragma once

#include <cmath>
#include <string>
#include <vector>

typedef std::vector<std::vector<double>> vector_t;

// Jakobi Verfahren
template <typename Func>
vector_t jakobi(vector_t     u,                // Eingabevector, mit Rand
                Func         f,                // Eingabefunktion
                const double h,                // Feinheit des Gitters
                int          &iteration_count, // Out-Variable (profiling)
                const double change_threshold, // Abbruch, wenn Änderung kleiner Wert
                const int    max_iterations)   // Abbruch, wenn Anzahl der Iterationen erreicht
{
    iteration_count = 0;
    auto u_old = u; // Kopie
    bool running = true;
    const int size = u.size() - 1;

    while (running && iteration_count < max_iterations) {
        ++iteration_count;
        std::swap(u, u_old);
        running = false;

        for (int i = 1; i < size; ++i) {
            for (int j = 1; j < size; ++j) {
                u[i][j] = (u_old[i][j - 1] + u_old[i - 1][j]
                         + u_old[i][j + 1] + u_old[i + 1][j]
                         + h * h * f(i * h, j * h)) * 0.25;

                // Liegen noch Änderungen der Werte oberhalb des Schwellwertes?
                running = running || std::abs(u[i][j] - u_old[i][j]) > change_threshold;
            }
        }
    }

    return u;
}

// Gauß-Seidel Verfahren
template <typename Func>
vector_t gaussSeidel(vector_t     u,                // Eingabevector, mit Rand
                     Func         f,                // Eingabefunktion
                     const double h,                // Feinheit des Gitters
                     int          &iteration_count, // Out-Variable (profiling)
                     const double change_threshold, // Abbruch, wenn Änderung kleiner Wert
                     const int    max_iterations)   // Abbruch, wenn Anzahl der Iterationen erreicht
{
    iteration_count = 0;
    bool running = true;
    const int size = u.size() - 1;

    while (running && iteration_count < max_iterations) {
        ++iteration_count;
        running = false;

        for (int i = 1; i < size; ++i) {
            for (int j = 1; j < size; ++j) {
                auto u_old = u[i][j];
                u[i][j] = (u[i][j - 1] + u[i - 1][j]
                         + u[i][j + 1] + u[i + 1][j]
                         + h * h * f(i * h, j * h)) * 0.25;

                // Liegen noch Änderungen der Werte oberhalb des Schwellwertes?
                running = running || std::abs(u[i][j] - u_old) > change_threshold;
            }
        }
    }

    return u;
}

// Debug
#include <fstream>

void dump(const std::string &filename, const vector_t &v) {
    static int count = 0;
    std::ofstream fout(filename + std::to_string(count++) + ".csv");
    for (std::size_t i = 0; i < v.size(); ++i)
        for (std::size_t j = 0; j < v.size(); ++j)
            fout << i << ", " << j << ", " << v[i][j] << "\n";
}

template <typename Func>
vector_t mehrgitter(vector_t         u,  // Eingabevektor mit Rand
                    Func             f,  // Eingabefunktion
                    const int        z1, // Iterationen Phase 1
                    const int        z2, // Iterationen Phase 2
                    const double     h,  // Feinheit des Eingabegitters
                    const double     h_max, // Feinheit des gröbsten Gitters
                    const int        alpha, // Rekursionsverzweigungsbreite
                    std::vector<std::pair<std::string, int>> &iteration_count, // Out-Variable (profiling)
                    const double     change_threshold, // Abbruch, wenn Änderung kleiner Wert
                    const int        max_iterations)   // Abbruch, wenn Anzahl der Iterationen erreicht
{
    if (h >= h_max) {
        iteration_count.emplace_back("Main", 0);
        auto res = gaussSeidel(u, f, h, iteration_count.back().second,
                               change_threshold, max_iterations);
        // Debug
        dump("Mehrgitter", res);
        return res;
    }

    iteration_count.emplace_back("Down", 0);
    auto vh = gaussSeidel(u, f, h, iteration_count.back().second, change_threshold, z1);
    // Debug
    dump("Mehrgitter", vh);
    const int n_new = (u.size() - 1) / 2;
    vector_t v2h(n_new + 2, std::vector<double>(n_new + 2, 0.0));

    // Restriktion
    for (int i = 1; i <= n_new; ++i)
        for (int j = 1; j <= n_new; ++j)
            v2h[i][j] = 0.125 * (4*vh[2*i][2*j] + vh[2*i - 1][2*j] + vh[2*i + 1][2*j]
                                 + vh[2*i][2*j - 1] + vh[2+i][2*j + 1]);

    // Rekursion
    for (int i=0; i<alpha; ++i)
        if (i == 0 || 2*h < h_max) // at lowest level only recurse once
            v2h = mehrgitter(v2h, f, z1, z2, 2*h, h_max, alpha,
                             iteration_count, change_threshold, max_iterations);

    // Interpolation
    for (std::size_t i = 0; i < u.size() - 2; i+=2) {
        for (std::size_t j = 0; j < u.size() - 2; j+=2) {
            vh[i+1][j+1] = v2h[i/2+1][j/2+1];
            vh[i+1][j+2] = 0.5 * (v2h[i/2+1][j/2+1] + v2h[i/2+1][j/2+1]);
            vh[i+2][j+1] = 0.5 * (v2h[i/2+1][j/2+1] + v2h[i/2+2][j/2+1]);
            vh[i+2][j+2] = 0.25 * (v2h[i/2+1][j/2+1] + v2h[i/2+1][j/2+2] +
                                   v2h[i/2+2][j/2+1] + v2h[i/2+2][j/2+2]);
        }
    }

    iteration_count.emplace_back("Up", 0);
    vh = gaussSeidel(vh, f, h, iteration_count.back().second, change_threshold, z2);
    // Debug
    dump("Mehrgitter", vh);
    return vh;
}
