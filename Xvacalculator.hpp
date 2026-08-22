#ifndef XVACALCULATOR_HPP
#define XVACALCULATOR_HPP
#include "CreditCurve.hpp"
#include <vector>

// Computes CVA and FVA from an Expected Exposure (EE) profile,
// as produced by ExposureEngine::computeEE.
//
// IMPORTANT: ExposureEngine's EE(t) is the expected value AT time t
// (in time-t dollars) ¡ª it is not yet discounted back to today.
// Both computeCVA and computeFVA apply the additional exp(-r*t)
// discount factor internally to bring each time step's exposure
// back to present value before integrating.
class XVACalculator {
public:
    XVACalculator(const CreditCurve& counterpartyCurve, double recoveryRate, double r);

    // CVA = (1-R) * sum_j D(t_j) * EE(t_j) * (S(t_{j-1}) - S(t_j))
    // Standard reduced-form CVA: expected loss given default, integrated
    // against the marginal default probability over each time step.
    // eeCurve[j] must correspond to time t_j = j * dt.
    double computeCVA(const std::vector<double>& eeCurve, double dt) const;

    // FVA = sum_j D(t_j) * EE(t_j) * fundingSpread * dt * S(t_j)
    // Simplified funding cost: accrues on the exposure while the
    // counterparty has survived, at a flat funding spread over the
    // risk-free rate. (A fuller treatment would also account for the
    // bank's own funding curve term structure and negative exposure /
    // funding benefit ¡ª out of scope here, noted as a simplification.)
    double computeFVA(const std::vector<double>& eeCurve, double dt, double fundingSpread) const;

private:
    const CreditCurve& counterpartyCurve_;
    double recoveryRate_;
    double r_;
};
#endif