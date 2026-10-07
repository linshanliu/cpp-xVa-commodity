#ifndef PATHSIMULATOR_HPP
#define PATHSIMULATOR_HPP
#include "IStochasticModel.hpp"
#include "RandomGenerator.hpp"
#include <vector>

class PathSimulator {
public:
    PathSimulator(const IStochasticModel& model, int numSteps, double dt);

    std::vector<std::vector<double>> simulatePath(RandomGenerator& rng) const;
    std::vector<std::vector<double>> simulatePath(RandomGenerator& rng, const std::vector<double>& startState) const;

    // Runs the path deterministically from a pre-generated sequence of
    // normal draws, instead of pulling them from a RandomGenerator.
    // Used for antithetic variates: the SAME underlying draws (negated)
    // are reused to build a mirror path, which requires the ability to
    // "replay" a specific set of draws rather than generating fresh ones.
    // draws.size() must equal numSteps_, each draws[k] must have
    // stateDimension() entries.
    std::vector<std::vector<double>> simulatePathFromDraws(
        const std::vector<double>& startState,
        const std::vector<std::vector<double>>& draws) const;

    double getDt() const { return dt_; }
    int getNumSteps() const { return numSteps_; }
    int getStateDimension() const { return model_.stateDimension(); }
    std::vector<double> getInitialState() const { return model_.initialState(); }

private:
    const IStochasticModel& model_;
    int numSteps_;
    double dt_;
};
#endif