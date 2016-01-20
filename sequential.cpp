#include <cmath>
#include <iomanip>
#include <iostream>
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

        for (int i = 1; i < u.size() - 1; ++i) {
            for (int j = 1; j < u.size() - 1; ++j) {
                u[i][j] = (u_old[i][j - 1] + u_old[i - 1][j]
                         + u_old[i][j + 1] + u_old[i + 1][j]
                         + h * h * f(i * h, j * h)) / 4;

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

        for (int i = 1; i < u.size() - 1; ++i) {
            for (int j = 1; j < u.size() - 1; ++j) {
                u[i][j] = (u    [i][j - 1] + u    [i - 1][j]
                         + u_old[i][j + 1] + u_old[i + 1][j]
                         + h * h * f(i * h, j * h)) / 4;

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
    std::cout << "Starting with h=" << h << "\n";
    if (h >= h_max)
    {
        return gauss_seidel(u, f, h, 0.00001, 1000000);
    }

    auto vh = gauss_seidel(u, f, h, 0.00001, z1);
    std::cout << z1 << " iterations of gauss-seidel done\n";
    const int n_new = (u.size() - 1) / 2;
    std::vector<std::vector<double>> v2h(n_new + 2, std::vector<double>(n_new + 2, 0.0));

    // Restriktion
    std::cout << "n_new=" << n_new << "\n";
    for (int i = 1; i <= n_new; ++i)
        for (int j = 1; j <= n_new; ++j)
            v2h[i][j] = 0.125 * (4*vh[2*i][2*j] + vh[2*i - 1][2*j] + vh[2*i + 1][2*j]
                                 + vh[2*i][2*j - 1] + vh[2+i][2*j + 1]);
    std::cout << "Restriction done\n";

    // Rekursion
    for (int i=0; i<alpha; i++)
        v2h = mehrgitter(v2h, f, z1, z2, 2*h, h_max, alpha);
    std::cout << "Recursion done\n";

    // Interpolation
    for (int i = 1; i < u.size() - 1; ++i)
        for (int j = 1; j < u.size() - 1; ++j)
            vh[i][j] +=  0.25 * (v2h[i/2][j/2] + v2h[i/2][j/2 + 1] +
                                 v2h[i/2 + 1][j/2] + v2h[i/2 + 1][j/2 + 1]);
    std::cout << "Interpolation done\n";

    return gauss_seidel(vh, f, h, 0.00001, z2);
}

int main(int argc, char **argv) {
    // Eingabefunktion "f", Aufgabe 1.1 a)
    auto f = [](double x, double y) -> double {
        return 32 * (x * (1 - x) + y * (1 - y));
    };

    // Feinheit "h"
    const int    n = 8;
    const int    n_max = 2;
    const double h = 1.0 / n;
    const double h_max = 1.0 / n_max;

    // Iterations- und Rekursionstiefen
    const int z1 = 100;
    const int z2 = 100;
    const int alpha = 1;

    // "Zufälliger" Startvektor "u_0", Rand inklusive
    std::vector<std::vector<double>> u(n + 2, std::vector<double>(n + 2, 0.0));
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= n; ++j)
            u[i][j] = 1.0;

    // Wende Verfahren an
    std::vector<std::pair<std::string, vector_t>> results;
    const double change_threshold = 0.00001;
    results.emplace_back("Jakobi",      jakobi      (u, f, h, change_threshold));
    results.emplace_back("Gauß-Seidel", gauss_seidel(u, f, h, change_threshold, 1000000));
    results.emplace_back("Mehrgitter",  mehrgitter  (u, f, z1, z2, h, h_max, alpha));

    // Ausgabe Verfahren
    std::cout << std::fixed << std::setprecision(4);
    for (auto &res : results) {
        std::cout << "Berechnetes Ergebnis, " << res.first << ":\n";
        for (int i = 0; i < n+2; ++i) {
            for (int j = 0; j < n+2; ++j)
                std::cout << " " << res.second[i][j];
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    // Ausgabe Lösung
    std::cout << "Analytische Lösung:\n";
    auto _u = [](double x, double y) -> double {
        return 16 * x * (1 - x) * y * (1 - y);
    };
    for (int i = 0; i < n+2; ++i) {
        for (int j = 0; j < n+2; ++j) {
            std::cout << " ";
            if (i == 0 || j == 0 || i == n+1 || j == n+1)
                std::cout << 0.0;
            else
                std::cout << _u(i * h, j * h);
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    // Ausgabe Fehler
    for (auto &res : results) {
        double err = 0.0;
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j)
                err += std::abs(res.second[i][j] - _u(i * h, j * h));
        std::cout << "Summe aller Fehler, " << res.first << ": " << err << "\n";
    }

    return 0;
}
