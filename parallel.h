#pragma once

#include <cmath>
#include <vector>

typedef std::vector<std::vector<double>> vector_t;

// Jakobi Verfahren
template <typename Func>
vector_t jakobiParallel(vector_t     u,                // Eingabevector, mit Rand
                        Func         f,                // Eingabefunktion
                        const double h,                // Feinheit des Gitters
                        int          &iteration_count, // Out-Variable (profiling)
                        const double change_threshold, // Abbruch, wenn Änderung kleiner Wert
                        const int    max_iterations)   // Abbruch, wenn Anzahl der Iterationen erreicht
{
    auto u_old = u; // Kopie
    bool running = true;
    iteration_count = 0;
    while (running && iteration_count++ < max_iterations) {
        std::swap(u, u_old);
        running = false;

        #pragma omp parallel for schedule(static, u.size() - 2), collapse(2)
        for (std::size_t i = 1; i < u.size() - 1; ++i) {
            for (std::size_t j = 1; j < u.size() - 1; ++j) {
                u[i][j] = (u_old[i][j - 1] + u_old[i - 1][j]
                         + u_old[i][j + 1] + u_old[i + 1][j]
                         + h * h * f(i * h, j * h)) * 0.25;

                // Liegen noch Änderungen der Werte oberhalb des Schwellwertes?
                if (std::abs(u[i][j] - u_old[i][j]) > change_threshold)
                    running = true; // Unsynchronisiert. Sollte kein Problem (bzgl. Korrektheit) sein.
            }
        }
    }
    --iteration_count; // Fix count

    return u;
}

// Gauß-Seidel Verfahren
template <typename Func>
vector_t gaussSeidelParallel(vector_t     u,                // Eingabevector, mit Rand
                             Func         f,                // Eingabefunktion
                             const double h,                // Feinheit des Gitters
                             int          &iteration_count, // Out-Variable (profiling)
                             const double change_threshold, // Abbruch, wenn Änderung kleiner Wert
                             const int    max_iterations)   // Abbruch, wenn Anzahl der Iterationen erreicht
{
    iteration_count = 0;
    // Vorlauf.
    for (std::size_t step = 0; step < u.size() - 2; ++step) {
        #pragma omp parallel for
        for (std::size_t offset = 0; offset <= step; ++offset) {
            const int i = 1 + step - offset;
            const int j = 1 + offset;

            u[i][j] = (u[i][j - 1] + u[i - 1][j]
                     + u[i][j + 1] + u[i + 1][j]
                     + h * h * f(i * h, j * h)) * 0.25;
        }
    }

    // Die erste halbe Iteration ist fast vollzogen. Mit der Diagonalen geht
    // es weiter und ab jetzt mit Anzahl der Threads = "u.size() ohne Rand".
    bool running = true;
    while (running && iteration_count++ < max_iterations) {
        running = false;

        for (std::size_t step = 0; step < u.size() - 2; ++step) {
            #pragma omp parallel for
            for (std::size_t offset = 0; offset < u.size() - 2; ++offset) {
                const int i = 1 + ((step - offset + u.size() - 2) % (u.size() - 2));
                const int j = 1 + offset;

                auto u_old = u[i][j];
                u[i][j] = (u[i][j - 1] + u[i - 1][j]
                         + u[i][j + 1] + u[i + 1][j]
                         + h * h * f(i * h, j * h)) * 0.25;

                // Liegen noch Änderungen der Werte oberhalb des Schwellwertes?
                if (std::abs(u[i][j] - u_old) > change_threshold)
                    running = true; // Unsynchronisiert. Sollte kein Problem (bzgl. Korrektheit) sein.
            }
        }
    }
    --iteration_count; // Fix count

    // Das Ergebnis ist nicht identisch mit der sequenziellen Version, da immer
    // zwei Iterationen parallel auf einem Vektor durchgeführt werden. Der
    // Vektor hat demnach nie den Zustand wie nach genau einer Iteration.
    return u;
}

template <typename Func>
vector_t mehrgitterParallel(vector_t         u,  // Eingabevektor mit Rand
                            Func             f,  // Eingabefunktion
                            const int        z1, // Iterationen Phase 1
                            const int        z2, // Iterationen Phase 2
                            const double     h,  // Feinheit des Eingabegitters
                            const double     h_max, // Feinheit des gröbsten Gitters
                            const int        alpha, // Rekursionsverzweigungsbreite
                            std::vector<int> &iteration_count, // Out-Variable (profiling)
                            const double     change_threshold, // Abbruch, wenn Änderung kleiner Wert
                            const int        max_iterations)   // Abbruch, wenn Anzahl der Iterationen erreicht
{
    if (h >= h_max) {
        iteration_count.push_back(0);
        return gaussSeidelParallel(u, f, h, iteration_count.back(),
                                   change_threshold, max_iterations);
    }

    iteration_count.push_back(0);
    auto vh = gaussSeidelParallel(u, f, h, iteration_count.back(), change_threshold, z1);
    const int n_new = (u.size() - 1) / 2;
    vector_t v2h(n_new + 2, std::vector<double>(n_new + 2, 0.0));

    // Restriktion
    for (int i = 1; i <= n_new; ++i)
        for (int j = 1; j <= n_new; ++j)
            v2h[i][j] = 0.125 * (4*vh[2*i][2*j] + vh[2*i - 1][2*j] + vh[2*i + 1][2*j]
                                 + vh[2*i][2*j - 1] + vh[2+i][2*j + 1]);

    // Rekursion
    for (int i=0; i<alpha; ++i)
        v2h = mehrgitterParallel(v2h, f, z1, z2, 2*h, h_max, alpha,
                                 iteration_count, change_threshold, max_iterations);

    // Interpolation
    for (std::size_t i = 1; i < u.size() - 1; ++i)
        for (std::size_t j = 1; j < u.size() - 1; ++j)
            vh[i][j] =  0.25 * (v2h[i/2][j/2] + v2h[i/2][j/2 + 1] +
                                v2h[i/2 + 1][j/2] + v2h[i/2 + 1][j/2 + 1]);

    iteration_count.push_back(0);
    return gaussSeidelParallel(vh, f, h, iteration_count.back(), change_threshold, z2);
}
