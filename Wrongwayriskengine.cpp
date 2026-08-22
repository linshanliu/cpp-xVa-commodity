#include "WrongWayRiskEngine.hpp"
#include "RandomGenerator.hpp"
#include <cmath>
#include <algorithm>

WrongWayRiskEngine::WrongWayRiskEngine(const PathSimulator& simulator, const IInstrument& instrument,
    double lambda0, double gamma, double recoveryRate, double r)
    : simulator_(simulator), instrument_(instrument),
    lambda0_(lambda0), gamma_(gamma), recoveryRate_(recoveryRate), r_(r) {
}

double WrongWayRiskEngine::computeCVA(int numPaths, unsigned int seedBase, int commodityIndex) const {
    int numSteps = simulator_.getNumSteps();
    double dt = simulator_.getDt();

    double cvaSumAcrossPaths = 0.0;

    for (int i = 0; i < numPaths; ++i) {
        RandomGenerator rng(seedBase + i);
        auto path = simulator_.simulatePath(rng);

        double S0 = path[0][2 * commodityIndex];   // reference level for this commodity

        double survival = 1.0;   // survival probability along THIS path, starts at 1 at t=0
        double cvaThisPath = 0.0;

        for (int j = 1; j <= numSteps; ++j) {
            double t_j = j * dt;
            double S_tj = path[j][2 * commodityIndex];

            // Path-dependent hazard rate: rises as S falls below S0 (for gamma > 0)
            double lambda_tj = lambda0_ * std::exp(-gamma_ * (S_tj - S0) / S0);

            double survivalPrev = survival;
            survival = survival * std::exp(-lambda_tj * dt);   // accumulate hazard along this path
            double marginalDefaultProb = survivalPrev - survival;

            double mtm = instrument_.markToMarket(t_j, path[j]);
            double exposure = std::max(mtm, 0.0);
            double discountFactor = std::exp(-r_ * t_j);

            cvaThisPath += discountFactor * exposure * marginalDefaultProb;
        }

        cvaSumAcrossPaths += cvaThisPath;
    }

    double avgCva = cvaSumAcrossPaths / numPaths;
    return (1.0 - recoveryRate_) * avgCva;
}