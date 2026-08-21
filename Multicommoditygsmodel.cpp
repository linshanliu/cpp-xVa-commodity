#include "MultiCommodityGSModel.hpp"
#include "CholeskyDecomposition.hpp"
#include <cmath>

MultiCommodityGSModel::MultiCommodityGSModel(const MultiCommodityParams& params)
    : params_(params), numCommodities_(params.numCommodities()) {
    // Precompute the Cholesky factor once, since the correlation matrix
    // is fixed for the lifetime of the model (avoids re-decomposing it
    // on every single step() call, which would be wasteful).
    choleskyL_ = CholeskyDecomposition::decompose(params_.correlationMatrix);
}

std::vector<double> MultiCommodityGSModel::initialState() const {
    std::vector<double> state;
    state.reserve(2 * numCommodities_);
    for (const auto& c : params_.commodities) {
        state.push_back(c.S0);
        state.push_back(c.delta0);
    }
    return state;
}

std::vector<double> MultiCommodityGSModel::step(
    const std::vector<double>& currentState, double dt, const std::vector<double>& normalDraws) const {

    // Correlate the independent draws using the precomputed Cholesky factor.
    // w is ordered the same way as the state: [w_S1, w_delta1, w_S2, w_delta2, ...]
    std::vector<double> w = CholeskyDecomposition::applyToVector(choleskyL_, normalDraws);

    std::vector<double> nextState(2 * numCommodities_);

    for (int i = 0; i < numCommodities_; ++i) {
        const SingleCommodityParams& p = params_.commodities[i];

        double S = currentState[2 * i];
        double delta = currentState[2 * i + 1];

        double wS = w[2 * i];
        double wDelta = w[2 * i + 1];

        double driftS = (params_.r - delta) * S * dt;
        double noiseS = p.sigmaS * S * std::sqrt(dt) * wS;

        double driftDelta = p.kappa * (p.alpha - delta) * dt;
        double noiseDelta = p.sigmaDelta * std::sqrt(dt) * wDelta;

        nextState[2 * i] = S + driftS + noiseS;
        nextState[2 * i + 1] = delta + driftDelta + noiseDelta;
    }

    return nextState;
}