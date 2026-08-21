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

private:
    const PathSimulator& simulator_;
    const IInstrument& instrument_;
    unsigned int seedBase_;
};
#endif