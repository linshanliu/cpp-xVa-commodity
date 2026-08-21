#ifndef EUROPEANCOMMODITYOPTION_HPP
#define EUROPEANCOMMODITYOPTION_HPP
#include "IInstrument.hpp"

enum class OptionType { Call, Put };

class EuropeanCommodityOption : public IInstrument {
public:
    EuropeanCommodityOption(double K, double T, OptionType type,
        double kappa, double alpha, double sigmaS,
        double sigmaDelta, double rho, double r);

    // Black-76 style price on top of the Gibson-Schwartz forward.
    // At t >= T, returns the intrinsic payoff (undiscounted).
    double markToMarket(double t, const std::vector<double>& state) const override;
    double maturity() const override { return T_; }

    // Forward price F(t,T) under Gibson-Schwartz (same formula as CommodityForward)
    double forwardPrice(double t, const std::vector<double>& state) const;

    // Total variance of ln(S_T) conditional on F_t, i.e. V(tau) used in Black-76
    double totalVariance(double t) const;

    double payoff(double S_T) const;

private:
    double K_, T_;
    OptionType type_;
    double kappa_, alpha_, sigmaS_, sigmaDelta_, rho_, r_;

    static double normCdf(double x);
};
#endif