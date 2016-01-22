#pragma once

#include <cmath>
#include <vector>

typedef std::vector<std::vector<double>> vector_t;

// Jakobi Verfahren
template <typename Func>
vector_t jakobi(vector_t     u,                // Eingabevector, mit Rand
                Func         f,                // Eingabefunktion
                const double h,                // Feinheit des Gitters
                const double change_threshold) // Abbruchkriterium
{
    bool running = true;
    while (running) {
        auto u_old = u; // Kopie
        running = false;

        for (std::size_t i = 1; i < u.size() - 1; ++i) {
            for (std::size_t j = 1; j < u.size() - 1; ++j) {
                u[i][j] = (u_old[i][j - 1] + u_old[i - 1][j]
                         + u_old[i][j + 1] + u_old[i + 1][j]
                         + h * h * f(i * h, j * h)) * 0.25;

                // Liegen noch Änderungen der Werte oberhalb des Schwellwertes?
                running = running || (std::abs(u[i][j] - u_old[i][j]) > change_threshold);
            }
        }
    }

    return u;
}

// Gauß-Seidel Verfahren
template <typename Func>
vector_t gauss_seidel(vector_t     u,                // Eingabevector, mit Rand
                      Func         f,                // Eingabefunktion
                      const double h,                // Feinheit des Gitters
                      const double change_threshold, // Abbruchkriterium
                      const int    max_iterations)   // Abbruchkriterium
{
    bool running = true;
    int iterations = 0;
    while (running && iterations++ < max_iterations) {
        auto u_old = u; // Kopie
        running = false;

        for (std::size_t i = 1; i < u.size() - 1; ++i) {
            for (std::size_t j = 1; j < u.size() - 1; ++j) {
                u[i][j] = (u    [i][j - 1] + u    [i - 1][j]
                         + u_old[i][j + 1] + u_old[i + 1][j]
                         + h * h * f(i * h, j * h)) * 0.25;

                // Liegen noch Änderungen der Werte oberhalb des Schwellwertes?
                running = running || (std::abs(u[i][j] - u_old[i][j]) > change_threshold);
            }
        }
    }

    return u;
}

template <typename Func>
vector_t mehrgitter(vector_t     u, // Eingabevektor mit Rand
                    Func         f, // Eingabefunktion
                    const int    z1, // Iterationen Phase 1
                    const int    z2, // Iterationen Phase 2
                    const double h, // Feinheit des Eingabegitters
                    const double h_max, // Feinheit des gröbsten Gitters
                    const int    alpha) // Rekursionsverzweigungsbreite
{
    if (h >= h_max)
    {
        return gauss_seidel(u, f, h, 0.00001, 1000000);
    }

    auto vh = gauss_seidel(u, f, h, 0.00001, z1);
    const int n_new = (u.size() - 1) / 2;
    std::vector<std::vector<double>> v2h(n_new + 2, std::vector<double>(n_new + 2, 0.0));

    // Restriktion
    for (int i = 1; i <= n_new; ++i)
        for (int j = 1; j <= n_new; ++j)
            v2h[i][j] = 0.125 * (4*vh[2*i][2*j] + vh[2*i - 1][2*j] + vh[2*i + 1][2*j]
                                 + vh[2*i][2*j - 1] + vh[2+i][2*j + 1]);

    // Rekursion
    for (int i=0; i<alpha; i++)
        v2h = mehrgitter(v2h, f, z1, z2, 2*h, h_max, alpha);

    // Interpolation
    for (std::size_t i = 1; i < u.size() - 1; ++i)
        for (std::size_t j = 1; j < u.size() - 1; ++j)
            vh[i][j] =  0.25 * (v2h[i/2][j/2] + v2h[i/2][j/2 + 1] +
                                v2h[i/2 + 1][j/2] + v2h[i/2 + 1][j/2 + 1]);

    return gauss_seidel(vh, f, h, 0.00001, z2);
}
