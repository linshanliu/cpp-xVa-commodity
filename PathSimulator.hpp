#ifndef PATHSIMULATOR_HPP
#define PATHSIMULATOR_HPP

#include "IStochasticModel.hpp"
#include "RandomGenerator.hpp"
#include <vector>

class PathSimulator {
private:
	const IStochasticModel& model_;
	int numSteps_;
	double dt_;
public:
	PathSimulator(const IStochasticModel& model, int numSteps, double dt);

	std::vector<std::vector<double>> simulatePath(RandomGenerator& rng, const std::vector<double>& startState) const;

	std::vector<std::vector<double>> simulatePath(RandomGenerator& rng) const;
};


#endif // !PATHSIMULATOR_HPP