#include "GibsonSchwartzModel.hpp"
#include <cmath>


GibsonSchwartzModel:: GibsonSchwartzModel(const GibsonSchwartzParams& params) : params_(params){}



std::vector<double> GibsonSchwartzModel::step(
	const std::vector<double>& currentState, double dt, const std::vector<double>& normalDraws) const {

	auto [w1, w2] = correlateNormals(normalDraws[0], normalDraws[1]);

	double driftS = (params_.r - currentState[1]) * currentState[0] * dt;
	double noiseS = params_.sigmaS * currentState[0] * sqrt(dt) * w1;

	double driftDelta = params_.kappa * (params_.alpha - currentState[1]) * dt;
	double noiseDelta = params_.sigmaDelta * sqrt(dt) * w2;

	std::vector<double> res{ currentState[0] + driftS + noiseS, currentState[1] + driftDelta + noiseDelta };

	return res;
}


std::vector<double> GibsonSchwartzModel::initialState() const {
	return { params_.S0, params_.delta0 };
}

// convert normal distribution sample to corresponding Cholesky decomposition
std::pair<double, double> GibsonSchwartzModel::correlateNormals(double z1, double z2) const {
	double w1 = z1;
	double w2 = params_.rho* z1 + sqrt(1 - params_.rho * params_.rho) * z2;
	return { w1,w2 };
}
