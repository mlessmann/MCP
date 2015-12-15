#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc, char **argv) {
    // Eingabefunktion "f"
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

    // Jakobi Verfahren
    bool keep_running = true; // TODO: Was ist ein geeignetes Abbruchkriterium?
    while (keep_running) {
        auto u_old = u; // Kopie
        keep_running = false;

        for (int i = 1; i < n; ++i) {
            for (int j = 1; j < n; ++j) {
                u[i][j] = (u_old[i][j - 1] + u_old[i - 1][j]
                         + u_old[i][j + 1] + u_old[i + 1][j]
                         + h * h * f(i * h, j * h)) / 4;

                // Haben sich die Werte noch geändert?
                keep_running = keep_running || (std::abs(u[i][j] - u_old[i][j]) > 0.00001); // Toleranz
                //keep_running = keep_running || (u[i][j] != u_old[i][j]); // Keine Toleranz!
            }
        }
    }

    // Ausgabe Jakobi
    std::cout << "Jakobi Ergebnis:\n";
    std::cout << std::fixed << std::setprecision(4);
    for (int j = 0; j <= n; ++j) {
        for (int i = 0; i <= n; ++i)
            std::cout << " " << u[i][j];
        std::cout << "\n";
    }

    // Ausgabe Lösung
    std::cout << "\nGegebene Lösung:\n";
    auto _u = [](double x, double y) -> double {
        return 16 * x * (1 - x) * y * (1 - y);
    };
    for (int j = 0; j <= n; ++j) {
        for (int i = 0; i <= n; ++i) {
            std::cout << " ";
            if (i == 0 || j == 0 || i == n || j == n)
                std::cout << 0.0;
            else
                std::cout << _u(i * h, j * h);
        }
        std::cout << "\n";
    } 

    // Ausgabe Fehler
    double err = 0.0;
    for (int i = 1; i < n; ++i)
        for (int j = 1; j < n; ++j)
            err += std::abs(u[i][j] - _u(i * h, j * h));
    std::cout << "\nSumme aller Fehler: " << err << "\n";

    return 0;
}
