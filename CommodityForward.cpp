#include "CommodityForward.hpp"
#include <cmath>

CommodityForward::CommodityForward(double K, double T, double kappa, double alpha,
    double sigmaS, double sigmaDelta, double rho, double r)
    : K_(K), T_(T), kappa_(kappa), alpha_(alpha),
    sigmaS_(sigmaS), sigmaDelta_(sigmaDelta), rho_(rho), r_(r) {
}

double CommodityForward::forwardPrice(double t, const std::vector<double>& state) const {
    double S_t = state[0];
    double delta_t = state[1];
    double tau = T_ - t;

    if (tau <= 0.0) {
        // 已经到期,forward价格就是当前spot
        return S_t;
    }

    double expKappaTau = std::exp(-kappa_ * tau);
    double B = (1.0 - expKappaTau) / kappa_;   // 简写, 后面重复用到

    // E[ln S_T] 部分
    double meanLogS = std::log(S_t)
        + r_ * tau
        - delta_t * B
        - alpha_ * (tau - B)
        - 0.5 * sigmaS_ * sigmaS_ * tau;   // Ito correction term

    // Var[ln S_T] 部分 (凸性修正)
    double term1 = sigmaS_ * sigmaS_ * tau;

    double term2 = (sigmaDelta_ * sigmaDelta_ / (kappa_ * kappa_))
        * (tau - 2.0 * B + (1.0 - std::exp(-2.0 * kappa_ * tau)) / (2.0 * kappa_));

    double term3 = (2.0 * rho_ * sigmaS_ * sigmaDelta_ / kappa_) * (tau - B);

    double varLogS = term1 + term2 - term3;

    double logF = meanLogS + 0.5 * varLogS;

    return std::exp(logF);
}

double CommodityForward::markToMarket(double t, const std::vector<double>& state) const {
    double tau = T_ - t;
    if (tau <= 0.0) {
        // 已到期, MtM就是spot和strike的差(不需要贴现)
        return state[0] - K_;
    }

    double F = forwardPrice(t, state);
    double discountFactor = std::exp(-r_ * tau);

    return (F - K_) * discountFactor;
}