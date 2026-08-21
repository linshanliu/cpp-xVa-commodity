#include "AmericanOptionPricer.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

AmericanOptionPricer::AmericanOptionPricer(const PathSimulator& simulator, double K, AmOptionType type, double r)
    : simulator_(simulator), K_(K), type_(type), r_(r) {
}

double AmericanOptionPricer::payoff(double S) const {
    if (type_ == AmOptionType::Call) {
        return std::max(S - K_, 0.0);
    }
    else {
        return std::max(K_ - S, 0.0);
    }
}

std::vector<double> AmericanOptionPricer::computeBasis(double sStd, double deltaStd) const {
    // [1, S, S^2, delta, S*delta]
    return { 1.0, sStd, sStd * sStd, deltaStd, sStd * deltaStd };
}

std::vector<double> AmericanOptionPricer::solveLinearSystem(
    std::vector<std::vector<double>> A, std::vector<double> b) const {

    int n = static_cast<int>(b.size());

    // Forward elimination with partial pivoting
    for (int col = 0; col < n; ++col) {
        // Find pivot row (largest absolute value in this column, at or below current row)
        int pivotRow = col;
        double maxAbs = std::abs(A[col][col]);
        for (int row = col + 1; row < n; ++row) {
            if (std::abs(A[row][col]) > maxAbs) {
                maxAbs = std::abs(A[row][col]);
                pivotRow = row;
            }
        }

        if (maxAbs < 1e-12) {
            // Singular / degenerate system, cannot solve reliably
            return std::vector<double>(n, 0.0);
        }

        std::swap(A[col], A[pivotRow]);
        std::swap(b[col], b[pivotRow]);

        // Eliminate this column from rows below
        for (int row = col + 1; row < n; ++row) {
            double factor = A[row][col] / A[col][col];
            for (int c = col; c < n; ++c) {
                A[row][c] -= factor * A[col][c];
            }
            b[row] -= factor * b[col];
        }
    }

    // Back substitution
    std::vector<double> x(n, 0.0);
    for (int row = n - 1; row >= 0; --row) {
        double sum = b[row];
        for (int c = row + 1; c < n; ++c) {
            sum -= A[row][c] * x[c];
        }
        x[row] = sum / A[row][row];
    }

    return x;
}

std::vector<double> AmericanOptionPricer::fitRegression(
    const std::vector<double>& sStdValues,
    const std::vector<double>& deltaStdValues,
    const std::vector<double>& yValues) const {

    int n = static_cast<int>(sStdValues.size());

    // Build normal equations: (Phi^T Phi) c = Phi^T y
    // where Phi is the n x basisDim_ design matrix
    std::vector<std::vector<double>> AtA(basisDim_, std::vector<double>(basisDim_, 0.0));
    std::vector<double> Atb(basisDim_, 0.0);

    for (int i = 0; i < n; ++i) {
        std::vector<double> phi = computeBasis(sStdValues[i], deltaStdValues[i]);
        double y = yValues[i];

        for (int r = 0; r < basisDim_; ++r) {
            Atb[r] += phi[r] * y;
            for (int c = 0; c < basisDim_; ++c) {
                AtA[r][c] += phi[r] * phi[c];
            }
        }
    }

    return solveLinearSystem(AtA, Atb);
}

double AmericanOptionPricer::evalBasis(const std::vector<double>& coeffs, double sStd, double deltaStd) const {
    std::vector<double> phi = computeBasis(sStd, deltaStd);
    double result = 0.0;
    for (int i = 0; i < basisDim_; ++i) {
        result += coeffs[i] * phi[i];
    }
    return result;
}

double AmericanOptionPricer::price(int numPaths, unsigned int seedBase, bool verbose) const {
    int numSteps = simulator_.getNumSteps();
    double dt = simulator_.getDt();

    // Step 1: simulate all paths, store full state [S, delta] at every time step
    std::vector<std::vector<double>> pathsS(numPaths);
    std::vector<std::vector<double>> pathsDelta(numPaths);

    for (int i = 0; i < numPaths; ++i) {
        RandomGenerator rng(seedBase + i);
        auto path = simulator_.simulatePath(rng);
        pathsS[i].reserve(numSteps + 1);
        pathsDelta[i].reserve(numSteps + 1);
        for (const auto& state : path) {
            pathsS[i].push_back(state[0]);
            pathsDelta[i].push_back(state[1]);
        }
    }

    // cashFlow[i] = discounted cash flow for path i, referenced to the
    // time step currently being processed in the backward induction
    std::vector<double> cashFlow(numPaths);
    for (int i = 0; i < numPaths; ++i) {
        cashFlow[i] = payoff(pathsS[i][numSteps]);
    }

    // Step 2: backward induction
    for (int j = numSteps - 1; j >= 1; --j) {

        std::vector<int> itmIndices;
        std::vector<double> sValues, deltaValues, yValues;

        for (int i = 0; i < numPaths; ++i) {
            double S = pathsS[i][j];
            double immediateExercise = payoff(S);
            if (immediateExercise > 0.0) {
                itmIndices.push_back(i);
                sValues.push_back(S);
                deltaValues.push_back(pathsDelta[i][j]);
                yValues.push_back(cashFlow[i] * std::exp(-r_ * dt));
            }
        }

        if (static_cast<int>(sValues.size()) < 2 * basisDim_) {
            // Not enough ITM paths to fit a stable regression this step
            for (int i = 0; i < numPaths; ++i) {
                cashFlow[i] *= std::exp(-r_ * dt);
            }
            continue;
        }

        // Standardize both S and delta before regression, to keep the
        // normal equations well-conditioned
        double sMean = std::accumulate(sValues.begin(), sValues.end(), 0.0) / sValues.size();
        double dMean = std::accumulate(deltaValues.begin(), deltaValues.end(), 0.0) / deltaValues.size();

        double sSqSum = 0.0, dSqSum = 0.0;
        for (double s : sValues) sSqSum += (s - sMean) * (s - sMean);
        for (double d : deltaValues) dSqSum += (d - dMean) * (d - dMean);

        double sStd = std::sqrt(sSqSum / sValues.size());
        double dStd = std::sqrt(dSqSum / deltaValues.size());
        if (sStd < 1e-8) sStd = 1.0;
        if (dStd < 1e-8) dStd = 1.0;

        std::vector<double> sStdValues, dStdValues;
        sStdValues.reserve(sValues.size());
        dStdValues.reserve(deltaValues.size());
        for (double s : sValues) sStdValues.push_back((s - sMean) / sStd);
        for (double d : deltaValues) dStdValues.push_back((d - dMean) / dStd);

        std::vector<double> coeffs = fitRegression(sStdValues, dStdValues, yValues);

        std::vector<bool> exercised(numPaths, false);
        int exerciseCount = 0;

        for (size_t k = 0; k < itmIndices.size(); ++k) {
            int i = itmIndices[k];
            double S = pathsS[i][j];
            double delta = pathsDelta[i][j];
            double immediateExercise = payoff(S);

            double sStdVal = (S - sMean) / sStd;
            double dStdVal = (delta - dMean) / dStd;
            double continuationValue = evalBasis(coeffs, sStdVal, dStdVal);

            if (immediateExercise > continuationValue) {
                cashFlow[i] = immediateExercise;
                exercised[i] = true;
                exerciseCount++;
            }
        }

        /*if (verbose && !itmIndices.empty()) {
            std::cout << "step j=" << j
                << ", ITM paths=" << itmIndices.size()
                << ", exercised=" << exerciseCount
                << " (" << (100.0 * exerciseCount / itmIndices.size()) << "%)\n";
        }*/

        for (int i = 0; i < numPaths; ++i) {
            if (!exercised[i]) {
                cashFlow[i] *= std::exp(-r_ * dt);
            }
        }
    }

    // Step 3: discount remaining cash flows from t=1 back to t=0, then average
    double sum = 0.0;
    for (int i = 0; i < numPaths; ++i) {
        sum += cashFlow[i] * std::exp(-r_ * dt);
    }

    return sum / numPaths;
}