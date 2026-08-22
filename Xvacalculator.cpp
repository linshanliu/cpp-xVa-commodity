#include "XVACalculator.hpp"
#include <cmath>

XVACalculator::XVACalculator(const CreditCurve& counterpartyCurve, double recoveryRate, double r)
    : counterpartyCurve_(counterpartyCurve), recoveryRate_(recoveryRate), r_(r) {
}

double XVACalculator::computeCVA(const std::vector<double>& eeCurve, double dt) const {
    double cva = 0.0;
    int numSteps = static_cast<int>(eeCurve.size()) - 1;

    for (int j = 1; j <= numSteps; ++j) {
        double t_j = j * dt;
        double t_prev = (j - 1) * dt;

        double survivalPrev = counterpartyCurve_.survivalProbability(t_prev);
        double survivalNow = counterpartyCurve_.survivalProbability(t_j);
        double marginalDefaultProb = survivalPrev - survivalNow;   // P(default in (t_{j-1}, t_j])

        double discountFactor = std::exp(-r_ * t_j);   // bring EE(t_j) back to present value

        cva += discountFactor * eeCurve[j] * marginalDefaultProb;
    }

    return (1.0 - recoveryRate_) * cva;
}

double XVACalculator::computeFVA(const std::vector<double>& eeCurve, double dt, double fundingSpread) const {
    double fva = 0.0;
    int numSteps = static_cast<int>(eeCurve.size()) - 1;

    for (int j = 1; j <= numSteps; ++j) {
        double t_j = j * dt;

        double discountFactor = std::exp(-r_ * t_j);
        double survivalNow = counterpartyCurve_.survivalProbability(t_j);

        fva += discountFactor * eeCurve[j] * fundingSpread * dt * survivalNow;
    }

    return fva;
}