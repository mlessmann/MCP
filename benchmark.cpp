#include "sequential.h"
#include <time.h>
#include <iostream>

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

vector_t createStartVector(int n) {
    vector_t u(n + 2, std::vector<double>(n + 2, 0.0));
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= n; ++j)
            u[i][j] = 1.0;
    return u;
}

void jakobi(int nMin, int nMax) {
    std::cout << "Starte Jakobi Benchmark\n";

    for (int n = nMin; n <= nMax; n*=2) {
        auto startVector = createStartVector(n);
        auto analyticalResult = analytical(startVector, u);
        std::cout << "n=" << n << "\n";

        jakobi(startVector, f, 
    }
}
