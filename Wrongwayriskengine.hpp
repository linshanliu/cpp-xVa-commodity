#ifndef WRONGWAYRISKENGINE_HPP
#define WRONGWAYRISKENGINE_HPP
#include "PathSimulator.hpp"
#include "IInstrument.hpp"
#include <vector>

// Computes CVA via direct Monte Carlo, with a PATH-DEPENDENT hazard
// rate that links the counterparty's default intensity to the
// simulated commodity price on that same path:
//
//   lambda_i(t) = lambda0 * exp( -gamma * (S_i(t) - S0) / S0 )
//
// gamma > 0 models wrong-way risk: when the commodity price falls
// below its initial level (which, for a long forward, is exactly
// when exposure is elevated), the counterparty's hazard rate rises,
// so default becomes more likely precisely when exposure is high.
// Setting gamma = 0 recovers a flat hazard rate (no WWR), which lets
// this same engine be used as a Monte Carlo cross-check against
// XVACalculator::computeCVA (which assumes independence between
// exposure and credit).
//
// This engine cannot reuse the "compute EE curve, then apply a fixed
// CreditCurve" two-step approach, because WWR is precisely the
// statement that exposure and default probability are NOT
// independent across paths ¡ª the correlation has to be preserved by
// keeping exposure and hazard rate on the SAME simulated path.
class WrongWayRiskEngine {
public:
    WrongWayRiskEngine(const PathSimulator& simulator, const IInstrument& instrument,
        double lambda0, double gamma, double recoveryRate, double r);

    // commodityIndex selects which commodity's price drives the hazard
    // rate (index 0 for a single-commodity model).
    double computeCVA(int numPaths, unsigned int seedBase, int commodityIndex = 0) const;

private:
    const PathSimulator& simulator_;
    const IInstrument& instrument_;
    double lambda0_;
    double gamma_;
    double recoveryRate_;
    double r_;
};
#endif