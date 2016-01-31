#pragma once

#include <cmath>
#include <omp.h>
#include <string>
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
    iteration_count = 0;
    auto u_old = u; // Kopie
    bool running = true;
    const int size = u.size() - 1;

    while (running && iteration_count < max_iterations) {
        ++iteration_count;
        std::swap(u, u_old);
        running = false;

        #pragma omp parallel for reduction(||:running) collapse(2)
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
vector_t gaussSeidelParallel(vector_t     u,                // Eingabevector, mit Rand
                             Func         f,                // Eingabefunktion
                             const double h,                // Feinheit des Gitters
                             int          &iteration_count, // Out-Variable (profiling)
                             const double change_threshold, // Abbruch, wenn Änderung kleiner Wert
                             const int    max_iterations)   // Abbruch, wenn Anzahl der Iterationen erreicht
{
    iteration_count = 1; // Vorlauf + Nachlauf = 1 Iteration
    const int size = u.size() - 2;

    // Vorlauf. Die erste halbe Iteration wird vorgezogen, damit anschließend
    // gleichmäßig parallel gerechnet werden kann.
    for (int step = 0; step < size; ++step) {
        #pragma omp parallel for
        for (int offset = 0; offset <= step; ++offset) {
            const int i = 1 + step - offset;
            const int j = 1 + offset;

            u[i][j] = (u[i][j - 1] + u[i - 1][j]
                     + u[i][j + 1] + u[i + 1][j]
                     + h * h * f(i * h, j * h)) * 0.25;
        }
    }

    // Ermittle, ob/wie sich der Vektor gleichmäßig aufteilen lässt.
    const int thread_count = omp_get_num_threads();
    int chunk_width = 1;
    bool simple_run = false;
    while (size / chunk_width > thread_count) {
        chunk_width *= 2;
        if (size % chunk_width != 0) {
            // Vektor lässt sich nicht (weiter) gleichmäßig unterteilen.
            // Fallback zur einfachen, parallelen Verarbeitung.
            simple_run = true;
            break;
        }
    }

    // Debug: Erzwinge simple_run. TODO: Entfernen.
    simple_run = true;

    // Wenn sich der Vektor nicht gut oder nur sehr fein aufteilen lässt,
    // verfahre nach dem einfachen, parallelen Verfahren.
    if (simple_run || chunk_width < 2) {
        bool running = true;
        while (running && iteration_count < max_iterations) {
            ++iteration_count;
            running = false;

            for (int step = 0; step < size; ++step) {
                #pragma omp parallel for reduction(||:running)
                for (int offset = 0; offset < size; ++offset) {
                    const int i = 1 + ((step - offset + size) % size);
                    const int j = 1 + offset;

                    auto u_old = u[i][j];
                    u[i][j] = (u[i][j - 1] + u[i - 1][j]
                             + u[i][j + 1] + u[i + 1][j]
                             + h * h * f(i * h, j * h)) * 0.25;

                    // Liegen noch Änderungen der Werte oberhalb des Schwellwertes?
                    running = running || std::abs(u[i][j] - u_old) > change_threshold ;
                }
            }
        }
    }
    // Wenn sich der Vektor gut aufteilen lässt, bilde in sich datenabhängige
    // Gruppen, die zu einer besseren Auslastung der Threads führen sollen.
    else {
        const int chunk_count = size / chunk_width;
        bool running = true;
        while (running && iteration_count < max_iterations) {
            ++iteration_count;
            running = false;

            for (int step = 0; step < size; step += chunk_width) {
                // Achtung: Hier kein collapse(2) einfügen!
                #pragma omp parallel for reduction(||:running)
                for (int chunk = 0; chunk < chunk_count; ++chunk) {
                    // TODO: Dreieck 1
                }

                // Achtung: Hier kein collapse(2) einfügen!
                #pragma omp parallel for reduction(||:running)
                for (int chunk = 0; chunk < chunk_count; ++chunk) {
                    // TODO: Dreieck 2
                }
            }
        }
    }

    // Nachlauf. Die letzte halbe Iteration wird noch nachgezogen, um das
    // Ergebnis identisch mit der sequenziellen Version zu halten.
    for (int step = 0; step < size - 1; ++step) {
        #pragma omp parallel for
        for (int offset = 0; offset < (size - 1) - step; ++offset) {
            const int i = size - offset;
            const int j = 2 + step + offset;

            u[i][j] = (u[i][j - 1] + u[i - 1][j]
                     + u[i][j + 1] + u[i + 1][j]
                     + h * h * f(i * h, j * h)) * 0.25;
        }
    }

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
                            std::vector<std::pair<std::string, int>> &iteration_count, // Out-Variable (profiling)
                            const double     change_threshold, // Abbruch, wenn Änderung kleiner Wert
                            const int        max_iterations)   // Abbruch, wenn Anzahl der Iterationen erreicht
{
    if (h >= h_max) {
        iteration_count.emplace_back("Main", 0);
        return gaussSeidelParallel(u, f, h, iteration_count.back().second,
                                   change_threshold, max_iterations);
    }

    iteration_count.emplace_back("Down", 0);
    auto vh = gaussSeidelParallel(u, f, h, iteration_count.back().second, change_threshold, z1);
    const int n = u.size() - 2;
    const int n_new = n / 2;
    vector_t v2h(n_new + 2, std::vector<double>(n_new + 2, 0.0));

    // Restriktion
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 1; i <= n_new; ++i)
        for (int j = 1; j <= n_new; ++j)
            v2h[i][j] = 0.125 * (4*vh[2*i][2*j] + vh[2*i - 1][2*j] + vh[2*i + 1][2*j]
                                 + vh[2*i][2*j - 1] + vh[2+i][2*j + 1]);

    // Rekursion
    for (int i=0; i<alpha; ++i)
        if (i == 0 || 2*h < h_max) // at lowest level only recurse once
            v2h = mehrgitterParallel(v2h, f, z1, z2, 2*h, h_max, alpha,
                                     iteration_count, change_threshold, max_iterations);

    // Interpolation
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 0; i < n; i+=2) {
        for (int j = 0; j < n; j+=2) {
            vh[i+1][j+1] = v2h[i/2+1][j/2+1];
            vh[i+1][j+2] = 0.5 * (v2h[i/2+1][j/2+1] + v2h[i/2+1][j/2+1]);
            vh[i+2][j+1] = 0.5 * (v2h[i/2+1][j/2+1] + v2h[i/2+2][j/2+1]);
            vh[i+2][j+2] = 0.25 * (v2h[i/2+1][j/2+1] + v2h[i/2+1][j/2+2] +
                                   v2h[i/2+2][j/2+1] + v2h[i/2+2][j/2+2]);
        }
    }

    iteration_count.emplace_back("Up", 0);
    return gaussSeidelParallel(vh, f, h, iteration_count.back().second, change_threshold, z2);
}
