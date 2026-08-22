#ifndef PORTFOLIO_HPP
#define PORTFOLIO_HPP
#include "IInstrument.hpp"
#include <vector>

// A netting set: aggregates multiple trades (possibly on different
// commodities within a multi-commodity state vector) into a single
// IInstrument. markToMarket returns the SUM of all trades' MtM ¡ª
// this sum is where netting happens: a positive-MtM trade and a
// negative-MtM trade on correlated commodities partially offset
// each other, which is exactly why netted exposure is typically
// much smaller than the sum of each trade's standalone exposure.
class Portfolio : public IInstrument {
public:
    // Trades are NOT owned by the Portfolio ¡ª the caller must keep
    // the underlying IInstrument objects alive for as long as this
    // Portfolio is used (same convention as PathSimulator holding a
    // const IStochasticModel& rather than owning the model).
    void addTrade(const IInstrument& trade);

    double markToMarket(double t, const std::vector<double>& state) const override;

    // The netting set must be simulated out to the latest maturity
    // among its trades, so ExposureEngine knows how far to simulate.
    double maturity() const override;

private:
    std::vector<const IInstrument*> trades_;
};
#endif