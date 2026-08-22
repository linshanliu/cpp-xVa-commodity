#ifndef COMMODITYFORWARD_HPP
#define COMMODITYFORWARD_HPP
#include "IInstrument.hpp"

class CommodityForward : public IInstrument {
public:
    // Original single-commodity constructor (commodityIndex defaults to 0,
    // so behavior on a 2-dim state vector is unchanged from before).
    CommodityForward(double K, double T, double kappa, double alpha,
        double sigmaS, double sigmaDelta, double rho, double r);

    // Multi-commodity constructor: commodityIndex selects which [S, delta]
    // pair this trade reads from a larger state vector
    // [S_0, delta_0, S_1, delta_1, ..., S_N, delta_N].
    CommodityForward(double K, double T, double kappa, double alpha,
        double sigmaS, double sigmaDelta, double rho, double r,
        int commodityIndex);

    double markToMarket(double t, const std::vector<double>& state) const override;
    double maturity() const override { return T_; }

    // forwardPrice/totalVariance operate on a single commodity's own
    // [S, delta] slice ¡ª callers of markToMarket don't need to slice
    // manually, that happens internally using commodityIndex_.
    double forwardPrice(double t, const std::vector<double>& ownState) const;

private:
    double K_;
    double T_;
    double kappa_;
    double alpha_;
    double sigmaS_;
    double sigmaDelta_;
    double rho_;
    double r_;
    int commodityIndex_;

    std::vector<double> extractOwnState(const std::vector<double>& fullState) const;
};
#endif