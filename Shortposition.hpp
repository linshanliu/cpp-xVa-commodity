#ifndef SHORTPOSITION_HPP
#define SHORTPOSITION_HPP
#include "IInstrument.hpp"

// Wraps any IInstrument and negates its mark-to-market, representing
// the opposite side of the trade (e.g. having sold a forward instead
// of bought it). This lets Portfolio net together long and short
// positions on different underlyings without needing separate
// "short" versions of every instrument type.
class ShortPosition : public IInstrument {
public:
    explicit ShortPosition(const IInstrument& underlying);

    double markToMarket(double t, const std::vector<double>& state) const override;
    double maturity() const override;

private:
    const IInstrument& underlying_;
};
#endif