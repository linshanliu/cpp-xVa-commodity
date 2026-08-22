#ifndef CREDITCURVE_HPP
#define CREDITCURVE_HPP

// Simplified counterparty credit curve using a CONSTANT hazard rate.
//
// A real desk would bootstrap a term structure of hazard rates from
// quoted CDS spreads across maturities (as noted as an explicit
// simplification in this project's scope doc: no live market data
// pipeline). Here, survival probability follows the standard
// reduced-form credit model:
//   S(t) = P(no default by time t) = exp(-hazardRate * t)
// A constant hazard rate corresponds to a flat CDS curve.
class CreditCurve {
public:
    explicit CreditCurve(double hazardRate);

    double survivalProbability(double t) const;
    double defaultProbability(double t) const;   // 1 - survivalProbability(t)

private:
    double hazardRate_;
};
#endif