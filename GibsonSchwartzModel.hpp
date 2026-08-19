#ifndef GIBSONSCHWARTZMODEL_HPP
#define GIBSONSCHWARTZMODEL_HPP

#include "IStochasticModel.hpp"
#include "ModelParameters.hpp"

#include <vector>


class GibsonSchwartzModel : public IStochasticModel {
private:
	GibsonSchwartzParams params_;

public:
	explicit GibsonSchwartzModel(const GibsonSchwartzParams& params);
	

	int stateDimension() const override { return 2; }


	std::vector<double> step(
		const std::vector<double>& currentState,
		double dt,
		const std::vector<double>& normalDraws
	)const override;


	std::vector<double> initialState() const override;

	// convert normal distribution sample to corresponding Cholesky decomposition
	std::pair<double, double> correlateNormals(double z1, double z2) const;
};



#endif // !GIBSONSCHWARTZMODEL_HPP
