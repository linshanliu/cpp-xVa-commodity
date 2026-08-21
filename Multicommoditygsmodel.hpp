#ifndef MULTICOMMODITYGSMODEL_HPP
#define MULTICOMMODITYGSMODEL_HPP
#include "IStochasticModel.hpp"
#include "MultiCommodityParams.hpp"
#include <vector>

// N-commodity generalization of GibsonSchwartzModel.
// State vector layout: [S_1, delta_1, S_2, delta_2, ..., S_N, delta_N]
// Each commodity follows its own Gibson-Schwartz dynamics; correlation
// across all 2N Brownian shocks (intra- and cross-commodity) is applied
// via a single Cholesky factor computed once at construction time.
class MultiCommodityGSModel : public IStochasticModel {
public:
    explicit MultiCommodityGSModel(const MultiCommodityParams& params);

    int stateDimension() const override { return 2 * numCommodities_; }

    std::vector<double> step(
        const std::vector<double>& currentState,
        double dt,
        const std::vector<double>& normalDraws
    ) const override;

    std::vector<double> initialState() const override;

private:
    MultiCommodityParams params_;
    int numCommodities_;
    std::vector<std::vector<double>> choleskyL_;   // precomputed once in the constructor
};
#endif