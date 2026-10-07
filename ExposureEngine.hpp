#ifndef EXPOSUREENGINE_HPP
#define EXPOSUREENGINE_HPP
#include "PathSimulator.hpp"
#include "IInstrument.hpp"
#include <vector>

class ExposureEngine {
public:
    ExposureEngine(const PathSimulator& simulator, const IInstrument& instrument, unsigned int seedBase);

    // Returns EE(t) curve: EE[j] is the expected exposure at time step j (t = j * dt)
    std::vector<double> computeEE(int numPaths) const;

    // Antithetic variates: for each pair of paths, one uses draws z,
    // the other uses -z (the SAME underlying random numbers, mirrored).
    // Since the two paths in a pair are negatively correlated, averaging
    // them reduces the variance of the EE estimate compared to using
    // the same total number of independent paths.
    // numPaths should be even (numPaths/2 antithetic pairs are simulated).
    std::vector<double> computeEEAntithetic(int numPaths) const;

    // Multithreaded version of computeEE: splits numPaths across
    // numThreads worker threads. Safe because Model::step() and
    // IInstrument::markToMarket() are const/stateless ¡ª each thread
    // uses its own RandomGenerator instance, so there is no shared
    // mutable state to synchronize.
    std::vector<double> computeEEParallel(int numPaths, int numThreads) const;

private:
    const PathSimulator& simulator_;
    const IInstrument& instrument_;
    unsigned int seedBase_;

    // Shared helper: accumulate exposure for one already-simulated path
    // into a running sum vector.
    void accumulateExposure(const std::vector<std::vector<double>>& path,
        double dt, std::vector<double>& sumExposure) const;
};
#endif