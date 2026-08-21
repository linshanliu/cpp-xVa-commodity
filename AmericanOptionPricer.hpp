#ifndef AMERICANOPTIONPRICER_HPP
#define AMERICANOPTIONPRICER_HPP
#include "PathSimulator.hpp"
#include "RandomGenerator.hpp"
#include <vector>

enum class AmOptionType { Call, Put };

// Longstaff-Schwartz pricer for an American-style option on the
// Gibson-Schwartz spot. Does NOT implement IInstrument, because
// the price requires a global regression across all paths at each
// time step, not a pointwise (t, state) -> price evaluation.
//
// Regression basis: [1, S, S^2, delta, S*delta]
// We regress on BOTH state variables (S and delta), not just S,
// because the drift of S depends on delta ((r - delta) * S term
// in the Gibson-Schwartz SDE), so the continuation value is a
// function of both factors, not S alone.
//
// --- Debugging note ---
// An earlier version of this pricer regressed the continuation value
// on S alone (basis [1, S, S^2]), which is a natural first attempt
// but turned out to be a genuine bug rather than a benign simplification.
// Because Gibson-Schwartz is a two-factor model in which the drift of S
// depends on delta, the true continuation value is a function of both
// state variables; omitting delta from the basis produced systematically
// biased continuation-value estimates and, in turn, suboptimal early-
// exercise decisions. The bug was caught using a model-independent
// sanity check: an American option's value can never be lower than the
// corresponding European option's value, since "never exercise early"
// is always an available strategy that reproduces the European payoff
// exactly. The S-only basis violated this bound (American priced below
// European), which is a stronger signal than "the numbers look a bit
// off" -- it identifies a real defect rather than Monte Carlo noise.
// Standardizing S alone and tightening the linear solver did not fix
// the violation, confirming the issue was model misspecification
// (a missing regressor) rather than numerical conditioning. Adding
// delta and the S*delta interaction term to the basis restored the
// theoretically required American >= European inequality.
class AmericanOptionPricer {
public:
    AmericanOptionPricer(const PathSimulator& simulator, double K, AmOptionType type, double r);

    // Returns the American option price at t=0.
    // Set verbose=true to print per-step exercise diagnostics.
    double price(int numPaths, unsigned int seedBase, bool verbose = false) const;

private:
    const PathSimulator& simulator_;
    double K_;
    AmOptionType type_;
    double r_;

    static const int basisDim_ = 5;

    double payoff(double S) const;

    // Basis functions evaluated at standardized (sStd, deltaStd)
    std::vector<double> computeBasis(double sStd, double deltaStd) const;

    // Solves A x = b for a square system via Gaussian elimination
    // with partial pivoting. Returns x, or a zero vector if the
    // system is (numerically) singular.
    std::vector<double> solveLinearSystem(std::vector<std::vector<double>> A, std::vector<double> b) const;

    // Fits continuation value ~ basis(S, delta) via least squares
    // (normal equations solved with solveLinearSystem).
    std::vector<double> fitRegression(
        const std::vector<double>& sStdValues,
        const std::vector<double>& deltaStdValues,
        const std::vector<double>& yValues) const;

    double evalBasis(const std::vector<double>& coeffs, double sStd, double deltaStd) const;
};
#endif