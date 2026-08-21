#include "EuropeanCommodityOption.hpp"
#include <cmath>
#include <algorithm>

EuropeanCommodityOption::EuropeanCommodityOption(double K, double T, OptionType type,
    double kappa, double alpha, double sigmaS, double sigmaDelta, double rho, double r)
    : K_(K), T_(T), type_(type), kappa_(kappa), alpha_(alpha),
    sigmaS_(sigmaS), sigmaDelta_(sigmaDelta), rho_(rho), r_(r) {
}

double EuropeanCommodityOption::normCdf(double x) {
    // Standard normal CDF via the complementary error function
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double EuropeanCommodityOption::forwardPrice(double t, const std::vector<double>& state) const {
    double S_t = state[0];
    double delta_t = state[1];
    double tau = T_ - t;

    if (tau <= 0.0) {
        return S_t;
    }

    double expKappaTau = std::exp(-kappa_ * tau);
    double B = (1.0 - expKappaTau) / kappa_;

    double meanLogS = std::log(S_t)
        + r_ * tau
        - delta_t * B
        - alpha_ * (tau - B)
        - 0.5 * sigmaS_ * sigmaS_ * tau;  // Ito correction

    double varLogS = totalVariance(t);

    double logF = meanLogS + 0.5 * varLogS;
    return std::exp(logF);
}

double EuropeanCommodityOption::totalVariance(double t) const {
    double tau = T_ - t;
    if (tau <= 0.0) return 0.0;

    double B = (1.0 - std::exp(-kappa_ * tau)) / kappa_;

    double term1 = sigmaS_ * sigmaS_ * tau;

    double term2 = (sigmaDelta_ * sigmaDelta_ / (kappa_ * kappa_))
        * (tau - 2.0 * B + (1.0 - std::exp(-2.0 * kappa_ * tau)) / (2.0 * kappa_));

    double term3 = (2.0 * rho_ * sigmaS_ * sigmaDelta_ / kappa_) * (tau - B);

    return term1 + term2 - term3;
}

double EuropeanCommodityOption::payoff(double S_T) const {
    if (type_ == OptionType::Call) {
        return std::max(S_T - K_, 0.0);
    }
    else {
        return std::max(K_ - S_T, 0.0);
    }
}

double EuropeanCommodityOption::markToMarket(double t, const std::vector<double>& state) const {
    double tau = T_ - t;

    if (tau <= 0.0) {
        // At/after maturity: intrinsic value using current spot
        return payoff(state[0]);
    }

    double F = forwardPrice(t, state);
    double V = totalVariance(t);         // total variance of ln(S_T)
    double sqrtV = std::sqrt(V);

    double d1 = (std::log(F / K_) + 0.5 * V) / sqrtV;
    double d2 = d1 - sqrtV;

    double discountFactor = std::exp(-r_ * tau);

    if (type_ == OptionType::Call) {
        return discountFactor * (F * normCdf(d1) - K_ * normCdf(d2));
    }
    else {
        return discountFactor * (K_ * normCdf(-d2) - F * normCdf(-d1));
    }
}