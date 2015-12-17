#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

typedef std::vector<std::vector<double>> vector_t;

// Jakobi Verfahren
template <typename Func>
vector_t jakobi(vector_t     u,         // Eingabevector
                Func         f,         // Eingabefunktion
                const double h,         // Feinheit des Gitters
                const double tolerance) // Abbruchkriterium
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

                // Haben sich die Werte noch innerhalb der Toleranz geändert?
                running = running || (std::abs(u[i][j] - u_old[i][j]) > tolerance);
            }
        }
    }

    return u;
}

// Gauß-Seidel Verfahren
template <typename Func>
vector_t gauss_seidel(vector_t     u,         // Eingabevector
                      Func         f,         // Eingabefunktion
                      const double h,         // Feinheit des Gitters
                      const double tolerance) // Abbruchkriterium
{
    bool running = true;
    while (running) {
        auto u_old = u; // Kopie
        running = false;

        for (int i = 1; i < u.size() - 1; ++i) {
            for (int j = 1; j < u.size() - 1; ++j) {
                u[i][j] = (u    [i][j - 1] + u    [i - 1][j]
                         + u_old[i][j + 1] + u_old[i + 1][j]
                         + h * h * f(i * h, j * h)) / 4;

                // Haben sich die Werte noch innerhalb der Toleranz geändert?
                running = running || (std::abs(u[i][j] - u_old[i][j]) > tolerance);
            }
        }
    }

    return u;
}

int main(int argc, char **argv) {
    // Eingabefunktion "f", Aufgabe 1.1 a)
    auto f = [](double x, double y) -> double {
        return 32 * (x * (1 - x) + y * (1 - y));
    };

    // Feinheit "h"
    const int    n = 10;
    const double h = 1.0 / n;

    // "Zufälliger" Startvektor "u_0", Rand inklusive
    std::vector<std::vector<double>> u(n + 1, std::vector<double>(n + 1, 0.0));
    for (int i = 1; i < n; ++i)
        for (int j = 1; j < n; ++j)
            u[i][j] = 1.0;

    // Wende Verfahren an
    std::map<std::string, vector_t> results;
    const double precision = 0.00001;
    results.emplace("Jakobi",      jakobi      (u, f, h, precision));
    results.emplace("Gauß-Seidel", gauss_seidel(u, f, h, precision));

    // Ausgabe Verfahren
    std::cout << std::fixed << std::setprecision(4);
    for (auto &res : results) {
        std::cout << "Berechnetes Ergebnis, " << res.first << ":\n";
        for (int i = 0; i <= n; ++i) {
            for (int j = 0; j <= n; ++j)
                std::cout << " " << res.second[i][j];
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    // Ausgabe Lösung
    std::cout << "Gegebene Lösung:\n";
    auto _u = [](double x, double y) -> double {
        return 16 * x * (1 - x) * y * (1 - y);
    };
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            std::cout << " ";
            if (i == 0 || j == 0 || i == n || j == n)
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
        for (int i = 1; i < n; ++i)
            for (int j = 1; j < n; ++j)
                err += std::abs(res.second[i][j] - _u(i * h, j * h));
        std::cout << "Summe aller Fehler, " << res.first << ": " << err << "\n";
    }

    return 0;
}
