#include "CommodityForward.hpp"
#include <cmath>

CommodityForward::CommodityForward(double K, double T, double kappa, double alpha,
    double sigmaS, double sigmaDelta, double rho, double r)
    : K_(K), T_(T), kappa_(kappa), alpha_(alpha),
    sigmaS_(sigmaS), sigmaDelta_(sigmaDelta), rho_(rho), r_(r), commodityIndex_(0) {
}

CommodityForward::CommodityForward(double K, double T, double kappa, double alpha,
    double sigmaS, double sigmaDelta, double rho, double r,
    int commodityIndex)
    : K_(K), T_(T), kappa_(kappa), alpha_(alpha),
    sigmaS_(sigmaS), sigmaDelta_(sigmaDelta), rho_(rho), r_(r), commodityIndex_(commodityIndex) {
}

std::vector<double> CommodityForward::extractOwnState(const std::vector<double>& fullState) const {
    // Pull out this trade's own [S, delta] pair from a (possibly larger)
    // multi-commodity state vector. For commodityIndex_=0 on a 2-element
    // state, this is just {fullState[0], fullState[1]} ¡ª identical to
    // the original single-commodity behavior.
    return { fullState[2 * commodityIndex_], fullState[2 * commodityIndex_ + 1] };
}

double CommodityForward::forwardPrice(double t, const std::vector<double>& ownState) const {
    double S_t = ownState[0];
    double delta_t = ownState[1];
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
        - 0.5 * sigmaS_ * sigmaS_ * tau;   // Ito correction

    double term1 = sigmaS_ * sigmaS_ * tau;
    double term2 = (sigmaDelta_ * sigmaDelta_ / (kappa_ * kappa_))
        * (tau - 2.0 * B + (1.0 - std::exp(-2.0 * kappa_ * tau)) / (2.0 * kappa_));
    double term3 = (2.0 * rho_ * sigmaS_ * sigmaDelta_ / kappa_) * (tau - B);
    double varLogS = term1 + term2 - term3;

    double logF = meanLogS + 0.5 * varLogS;
    return std::exp(logF);
}

double CommodityForward::markToMarket(double t, const std::vector<double>& state) const {
    std::vector<double> ownState = extractOwnState(state);

    double tau = T_ - t;
    if (tau <= 0.0) {
        return ownState[0] - K_;
    }

    double F = forwardPrice(t, ownState);
    double discountFactor = std::exp(-r_ * tau);

    return (F - K_) * discountFactor;
}