#include "CholeskyDecomposition.hpp"
#include <cmath>
#include <stdexcept>

std::vector<std::vector<double>> CholeskyDecomposition::decompose(const std::vector<std::vector<double>>& A) {
    int n = static_cast<int>(A.size());
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0.0;
            for (int k = 0; k < j; ++k) {
                sum += L[i][k] * L[j][k];
            }

            if (i == j) {
                double diag = A[i][i] - sum;
                if (diag <= 0.0) {
                    throw std::runtime_error(
                        "CholeskyDecomposition: matrix is not positive definite "
                        "(non-positive pivot encountered). Check the correlation matrix.");
                }
                L[i][j] = std::sqrt(diag);
            }
            else {
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }

    return L;
}

std::vector<double> CholeskyDecomposition::applyToVector(
    const std::vector<std::vector<double>>& L,
    const std::vector<double>& z) {

    int n = static_cast<int>(L.size());
    std::vector<double> w(n, 0.0);

    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int k = 0; k <= i; ++k) {   // L is lower triangular, so only k<=i contributes
            sum += L[i][k] * z[k];
        }
        w[i] = sum;
    }

    return w;
}