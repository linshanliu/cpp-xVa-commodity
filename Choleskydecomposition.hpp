#ifndef CHOLESKYDECOMPOSITION_HPP
#define CHOLESKYDECOMPOSITION_HPP
#include <vector>

// Computes the lower-triangular Cholesky factor L of a symmetric
// positive-definite matrix A, such that A = L * L^T.
//
// Used to correlate a vector of independent standard normal draws:
// given z ~ iid N(0,1), w = L * z produces w with covariance matrix A.
//
// This generalizes the 2x2 closed-form formula used for a single
// commodity (w1 = z1, w2 = rho*z1 + sqrt(1-rho^2)*z2) to an
// arbitrary N x N correlation/covariance matrix, needed once the
// model covers multiple commodities.
class CholeskyDecomposition {
public:
    // A must be square, symmetric, and positive definite.
    // Throws std::runtime_error if a non-positive pivot is encountered
    // (indicates A is not positive definite, e.g. an inconsistent
    // correlation matrix).
    static std::vector<std::vector<double>> decompose(const std::vector<std::vector<double>>& A);

    // Applies L to a vector z: returns w = L * z
    static std::vector<double> applyToVector(
        const std::vector<std::vector<double>>& L,
        const std::vector<double>& z);
};
#endif