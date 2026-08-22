#include "ShortPosition.hpp"

ShortPosition::ShortPosition(const IInstrument& underlying)
    : underlying_(underlying) {
}

double ShortPosition::markToMarket(double t, const std::vector<double>& state) const {
    return -underlying_.markToMarket(t, state);
}

double ShortPosition::maturity() const {
    return underlying_.maturity();
}