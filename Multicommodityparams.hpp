#ifndef MULTICOMMODITYPARAMS_HPP
#define MULTICOMMODITYPARAMS_HPP
#include <vector>

// Per-commodity parameters. Note: rho is no longer stored here since,
// with multiple commodities, correlation lives in a single full
// correlation matrix (see MultiCommodityParams below) that covers
// both intra-commodity (S vs delta) and cross-commodity correlations.
struct SingleCommodityParams {
    double S0;
    double delta0;
    double kappa;
    double alpha;
    double sigmaS;
    double sigmaDelta;
};

struct MultiCommodityParams {
    std::vector<SingleCommodityParams> commodities;

    // Full correlation matrix over all Brownian shocks, ordered as
    // [dW_S1, dW_delta1, dW_S2, dW_delta2, ..., dW_SN, dW_deltaN].
    // Size must be (2*numCommodities) x (2*numCommodities), symmetric,
    // with 1.0 on the diagonal.
    std::vector<std::vector<double>> correlationMatrix;

    double r;   // shared risk-free rate across commodities (simplification)

    int numCommodities() const { return static_cast<int>(commodities.size()); }
};
#endif