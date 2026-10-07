#include "ExposureEngine.hpp"
#include "RandomGenerator.hpp"
#include <algorithm>
#include <thread>

ExposureEngine::ExposureEngine(const PathSimulator& simulator, const IInstrument& instrument, unsigned int seedBase)
    : simulator_(simulator), instrument_(instrument), seedBase_(seedBase) {
}

void ExposureEngine::accumulateExposure(const std::vector<std::vector<double>>& path,
    double dt, std::vector<double>& sumExposure) const {
    int numSteps = static_cast<int>(path.size()) - 1;
    for (int j = 0; j <= numSteps; ++j) {
        double t = j * dt;
        double mtm = instrument_.markToMarket(t, path[j]);
        double exposure = std::max(mtm, 0.0);
        sumExposure[j] += exposure;
    }
}

std::vector<double> ExposureEngine::computeEE(int numPaths) const {
    int numSteps = simulator_.getNumSteps();
    double dt = simulator_.getDt();

    std::vector<double> sumExposure(numSteps + 1, 0.0);

    for (int i = 0; i < numPaths; ++i) {
        RandomGenerator rng(seedBase_ + i);
        auto path = simulator_.simulatePath(rng);
        accumulateExposure(path, dt, sumExposure);
    }

    std::vector<double> ee(numSteps + 1);
    for (int j = 0; j <= numSteps; ++j) {
        ee[j] = sumExposure[j] / numPaths;
    }
    return ee;
}

std::vector<double> ExposureEngine::computeEEAntithetic(int numPaths) const {
    int numSteps = simulator_.getNumSteps();
    double dt = simulator_.getDt();
    int stateDim = simulator_.getStateDimension();

    int numPairs = numPaths / 2;   // numPaths should be even; integer division drops any odd remainder

    std::vector<double> sumExposure(numSteps + 1, 0.0);
    std::vector<double> initState;   // will be set from the model via the first path's t=0 state

    for (int i = 0; i < numPairs; ++i) {
        RandomGenerator rng(seedBase_ + i);

        // Pre-generate the draws for this pair, so the SAME numbers
        // (negated) can be reused for the antithetic mirror path.
        std::vector<std::vector<double>> draws(numSteps);
        for (int j = 0; j < numSteps; ++j) {
            draws[j] = rng.generateNormals(stateDim);
        }

        // Base path uses draws as-is
        std::vector<double> startState = simulator_.getInitialState();
        auto pathA = simulator_.simulatePathFromDraws(startState, draws);

        // Antithetic path uses negated draws
        std::vector<std::vector<double>> negatedDraws(numSteps);
        for (int j = 0; j < numSteps; ++j) {
            negatedDraws[j].reserve(stateDim);
            for (int d = 0; d < stateDim; ++d) {
                negatedDraws[j].push_back(-draws[j][d]);
            }
        }
        auto pathB = simulator_.simulatePathFromDraws(startState, negatedDraws);

        accumulateExposure(pathA, dt, sumExposure);
        accumulateExposure(pathB, dt, sumExposure);
    }

    int totalPathsUsed = numPairs * 2;
    std::vector<double> ee(numSteps + 1);
    for (int j = 0; j <= numSteps; ++j) {
        ee[j] = sumExposure[j] / totalPathsUsed;
    }
    return ee;
}

std::vector<double> ExposureEngine::computeEEParallel(int numPaths, int numThreads) const {
    int numSteps = simulator_.getNumSteps();
    double dt = simulator_.getDt();

    // Each thread accumulates into its own private sum vector ¡ª no
    // shared mutable state, so no locking is needed. This is only
    // safe because Model::step() and Instrument::markToMarket() are
    // const/stateless by design (a deliberate architectural choice
    // made back in Phase 1), and each thread owns its own
    // RandomGenerator instance.
    std::vector<std::vector<double>> partialSums(numThreads, std::vector<double>(numSteps + 1, 0.0));

    auto worker = [&](int threadId, int startPath, int endPath) {
        for (int i = startPath; i < endPath; ++i) {
            RandomGenerator rng(seedBase_ + i);
            auto path = simulator_.simulatePath(rng);
            accumulateExposure(path, dt, partialSums[threadId]);
        }
        };

    std::vector<std::thread> threads;
    int pathsPerThread = numPaths / numThreads;
    int remainder = numPaths % numThreads;
    int start = 0;

    for (int t = 0; t < numThreads; ++t) {
        int count = pathsPerThread + (t < remainder ? 1 : 0);
        threads.emplace_back(worker, t, start, start + count);
        start += count;
    }
    for (auto& th : threads) {
        th.join();
    }

    std::vector<double> ee(numSteps + 1, 0.0);
    for (int j = 0; j <= numSteps; ++j) {
        double sum = 0.0;
        for (int t = 0; t < numThreads; ++t) {
            sum += partialSums[t][j];
        }
        ee[j] = sum / numPaths;
    }
    return ee;
}