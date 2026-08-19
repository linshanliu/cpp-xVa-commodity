#ifndef ISTOCHASTICMODEL_HPP
#define ISTOCHASTICMODEL_HPP

#include <vector>

class IStochasticModel {
public:
	virtual ~IStochasticModel() = default;


	// return the number of state dimension
	virtual int stateDimension() const = 0;

	// given the current state, step size and the normal realization
	// return the new state
	virtual std::vector<double> step(
		const std::vector<double>& currentState,
		double dt,
		const std::vector<double>& nomalDraws
	)const = 0;


	// return the initial state
	virtual std::vector<double> initialState() const = 0;
};





#endif // !ISTOCHASTICMODEL_HPP
