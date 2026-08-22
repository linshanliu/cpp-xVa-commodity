#include "Portfolio.hpp"
#include <algorithm>
#include <stdexcept>

void Portfolio::addTrade(const IInstrument& trade) {
    trades_.push_back(&trade);
}

double Portfolio::markToMarket(double t, const std::vector<double>& state) const {
    double sum = 0.0;
    for (const auto* trade : trades_) {
        sum += trade->markToMarket(t, state);
    }
    return sum;
}

double Portfolio::maturity() const {
    if (trades_.empty()) {
        throw std::runtime_error("Portfolio::maturity() called on an empty portfolio.");
    }
    double maxT = trades_[0]->maturity();
    for (const auto* trade : trades_) {
        maxT = std::max(maxT, trade->maturity());
    }
    return maxT;
}