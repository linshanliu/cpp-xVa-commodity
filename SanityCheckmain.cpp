#include "GibsonSchwartzModel.hpp"
#include "PathSimulator.hpp"
#include "RandomGenerator.hpp"
#include "CommodityForward.hpp"
#include <iostream>
#include <numeric>
#include <cmath>

int main() {
    // Define the parameters for the Gibson-Schwartz model
    // S0=100, delta0=0.05, kappa=1.2, alpha=0.08, sigmaS=0.3, sigmaDelta=0.1, rho=-0.5, r=0.03
    GibsonSchwartzParams params{ 100.0, 0.05, 1.2, 0.08, 0.3, 0.1, -0.5, 0.03 };

    // Create the Gibson-Schwartz model
    GibsonSchwartzModel model(params);

    double T = 10.0; // Total time horizon

    int numSteps = 500;
    double dt = 10.0 / numSteps;

    // Create the path simulator
    PathSimulator simulator(model, numSteps, dt);

    int numPaths = 1000;

    // Create a vector to store the terminal deltas
    std::vector<double> terminalDeltas;
    terminalDeltas.reserve(numPaths);

    // ============================================================
    // Sanity Check 1: Mean-reversion test
    // Simulate many paths and verify that E[delta_T] converges to
    // alpha (the long-run mean of the convenience yield) as t grows.
    // ============================================================

    // Simulate multiple paths and collect the terminal deltas
    for (int i = 0; i < numPaths; ++i) {
        RandomGenerator rng(i);
        auto path = simulator.simulatePath(rng);

        // Store the terminal delta (index [1]) from the last step (.back()) of the path
        terminalDeltas.push_back(path.back()[1]);
    }

    // accumulate(begin, end, initial value) -> sums all elements starting from initial value
    double sampleMean = std::accumulate(terminalDeltas.begin(), terminalDeltas.end(), 0.0) / numPaths;
    std::cout << "=== Sanity Check 1: Mean Reversion ===\n";
    std::cout << "Sample mean of delta_T: " << sampleMean << "\n";
    std::cout << "Target (alpha):         " << params.alpha << "\n\n";

    // ============================================================
    // Sanity Check 2: Moment matching test
    // Compare the sample variance of delta_T against the analytical
    // variance of an OU process at time T.
    // Var[delta_T] = (sigmaDelta^2 / 2*kappa) * (1 - exp(-2*kappa*T))
    // ============================================================

    double sampleVariance = 0.0;
    for (const auto& delta : terminalDeltas) {
        sampleVariance += (delta - sampleMean) * (delta - sampleMean);
    }
    sampleVariance = sampleVariance / (numPaths - 1); // unbiased sample variance (Bessel's correction)

    double analyticalVar = (params.sigmaDelta * params.sigmaDelta) / (2.0 * params.kappa)
        * (1.0 - std::exp(-2.0 * params.kappa * T));

    std::cout << "=== Sanity Check 2: Variance Matching ===\n";
    std::cout << "Sample variance of delta_T:     " << sampleVariance << "\n";
    std::cout << "Analytical variance of delta_T: " << analyticalVar << "\n\n";

    // ============================================================
    // Sanity Check 3: Correlation test
    // Verify that correlateNormals() produces two correlated normal
    // draws (w1, w2) whose sample correlation matches rho.
    // ============================================================

    RandomGenerator rngTest(999);
    int numSamples = 100000; // correlation needs a large sample size to be stable
    std::vector<double> w1samples, w2samples;
    w1samples.reserve(numSamples);
    w2samples.reserve(numSamples);

    for (int i = 0; i < numSamples; ++i) {
        auto z = rngTest.generateNormals(2);
        auto [w1, w2] = model.correlateNormals(z[0], z[1]);
        w1samples.push_back(w1);
        w2samples.push_back(w2);
    }

    // Sample correlation = Cov(w1,w2) / (std(w1) * std(w2))
    // Since w1, w2 are standard normal (unit variance), covariance IS the correlation here
    double covariance = 0.0;
    double mean_w1 = std::accumulate(w1samples.begin(), w1samples.end(), 0.0) / numSamples;
    double mean_w2 = std::accumulate(w2samples.begin(), w2samples.end(), 0.0) / numSamples;
    for (int i = 0; i < numSamples; ++i) {
        covariance += (w1samples[i] - mean_w1) * (w2samples[i] - mean_w2);
    }
    covariance /= numSamples;

    std::cout << "=== Sanity Check 3: Correlation ===\n";
    std::cout << "Sample covariance of w1 and w2: " << covariance << "\n";
    std::cout << "Target (rho):                   " << params.rho << "\n\n";

    // ============================================================
    // Sanity Check 4: Forward price validation
    // Compare the analytical Gibson-Schwartz forward price F(t,T)
    // against a nested Monte Carlo estimate: simulate paths starting
    // from (t0, stateAtT0) and take the sample mean of S_T.
    // ============================================================

    double t0 = 2.0;                                  // valuation time
    std::vector<double> stateAtT0 = { 110.0, 0.06 };   // assumed state [S_t, delta_t] at t0
    double T_fwd = 5.0;                                // forward maturity

    CommodityForward fwd(100.0, T_fwd, params.kappa, params.alpha,
        params.sigmaS, params.sigmaDelta, params.rho, params.r);

    double analyticalF = fwd.forwardPrice(t0, stateAtT0);
    std::cout << "=== Sanity Check 4: Forward Price Validation ===\n";
    std::cout << "Analytical forward price at t=" << t0 << ": " << analyticalF << "\n";

    // Monte Carlo cross-check: start from (t0, stateAtT0), simulate to T_fwd,
    // and take the sample mean of the terminal spot S_T
    int numStepsShort = 150;
    double dtShort = (T_fwd - t0) / numStepsShort;
    PathSimulator shortSimulator(model, numStepsShort, dtShort);

    int numPathsShort = 5000; // price has fatter tails than delta, so use more paths
    std::vector<double> terminalSpots;
    terminalSpots.reserve(numPathsShort);

    for (int i = 0; i < numPathsShort; ++i) {
        RandomGenerator rngShort(i + 10000); // offset seed to avoid reusing earlier draws
        auto path = shortSimulator.simulatePath(rngShort, stateAtT0);
        terminalSpots.push_back(path.back()[0]); // index [0] = S
    }

    double mcMean = std::accumulate(terminalSpots.begin(), terminalSpots.end(), 0.0) / numPathsShort;
    std::cout << "Monte Carlo mean of S_T:          " << mcMean << "\n";

    double relError = std::abs(mcMean - analyticalF) / analyticalF * 100.0;
    std::cout << "Relative error: " << relError << "%\n";

    return 0;
}