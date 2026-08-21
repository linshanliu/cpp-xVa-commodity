#include "ExposureEngine.hpp"
#include "RandomGenerator.hpp"
#include <algorithm>

ExposureEngine::ExposureEngine(const PathSimulator& simulator, const IInstrument& instrument, unsigned int seedBase)
    : simulator_(simulator), instrument_(instrument), seedBase_(seedBase) {
}

std::vector<double> ExposureEngine::computeEE(int numPaths) const {
    int numSteps = simulator_.getNumSteps();
    double dt = simulator_.getDt();

    // Accumulator for sum of exposures at each time step, across all paths
    std::vector<double> sumExposure(numSteps + 1, 0.0);

    for (int i = 0; i < numPaths; ++i) {
        RandomGenerator rng(seedBase_ + i);
        auto path = simulator_.simulatePath(rng);

        for (int j = 0; j <= numSteps; ++j) {
            double t = j * dt;
            double mtm = instrument_.markToMarket(t, path[j]);
            double exposure = std::max(mtm, 0.0);
            sumExposure[j] += exposure;
        }
    }

    // Average across paths to get EE(t)
    std::vector<double> ee(numSteps + 1);
    for (int j = 0; j <= numSteps; ++j) {
        ee[j] = sumExposure[j] / numPaths;
    }

    return ee;
}