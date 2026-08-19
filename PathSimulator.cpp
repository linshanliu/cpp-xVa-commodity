#include "PathSimulator.hpp"


PathSimulator::PathSimulator(const IStochasticModel& model, int numSteps, double dt)
	: model_(model), numSteps_(numSteps), dt_(dt) {
}


std::vector<std::vector<double>> PathSimulator::simulatePath(RandomGenerator& rng, const std::vector<double>& startState) const {
	int stateDim = model_.stateDimension();
	std::vector<double> currentState = startState;
	std::vector<std::vector<double>> path;
	path.reserve(numSteps_ + 1);
	path.push_back(currentState);
	for (int i = 0; i < numSteps_; i++) {
		std::vector<double> normals = rng.generateNormals(stateDim);
		std::vector<double> nextState = model_.step(currentState, dt_, normals);
		path.push_back(nextState);
		currentState = nextState;
	}
	return path;
}

std::vector<std::vector<double>> PathSimulator::simulatePath(RandomGenerator& rng) const {
	return simulatePath(rng, model_.initialState());
}