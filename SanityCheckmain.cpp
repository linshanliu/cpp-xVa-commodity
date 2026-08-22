#include "GibsonSchwartzModel.hpp"
#include "PathSimulator.hpp"
#include "RandomGenerator.hpp"
#include "CommodityForward.hpp"
#include "ExposureEngine.hpp"
#include "EuropeanCommodityOption.hpp"
#include "AmericanOptionPricer.hpp"
#include "MultiCommodityGSModel.hpp"
#include "Choleskydecomposition.hpp"
#include "Portfolio.hpp"
#include "ShortPosition.hpp"
#include "CreditCurve.hpp"
#include "Xvacalculator.hpp"
#include "WrongWayRiskEngine.hpp"

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


    // ============================================================
    // Phase 1 capstone: Expected Exposure (EE) profile for a single forward
    // ============================================================

    double forwardT = 5.0;   // maturity
    double strikeK = 100.0;  // agreed delivery price

    CommodityForward mainForward(strikeK, forwardT, params.kappa, params.alpha,
        params.sigmaS, params.sigmaDelta, params.rho, params.r);

    int eeSteps = 250;
    double eeDt = forwardT / eeSteps;
    PathSimulator eeSimulator(model, eeSteps, eeDt);

    ExposureEngine engine(eeSimulator, mainForward, 42);

    int eePaths = 5000;
    std::vector<double> eeCurve = engine.computeEE(eePaths);

    std::cout << "\n=== Expected Exposure Profile ===\n";
    // Print every 25th step to keep the output readable
    for (int j = 0; j <= eeSteps; j += 25) {
        double t = j * eeDt;
        std::cout << "t = " << t << "   EE(t) = " << eeCurve[j] << "\n";
    }

    // ============================================================
// Phase 2, Step 1: European option validation (Black-76 on GS forward)
// Cross-check against nested Monte Carlo: simulate to maturity,
// average the discounted payoff.
// ============================================================

    double optK = 105.0;
    double optT = 5.0;
    EuropeanCommodityOption callOpt(optK, optT, OptionType::Call,
        params.kappa, params.alpha, params.sigmaS, params.sigmaDelta, params.rho, params.r);

    double t0_opt = 2.0;
    std::vector<double> stateAtT0_opt = { 110.0, 0.06 };

    double analyticalOptPrice = callOpt.markToMarket(t0_opt, stateAtT0_opt);
    std::cout << "\n=== European Call Option Validation ===\n";
    std::cout << "Analytical (Black-76 on GS) price at t=2: " << analyticalOptPrice << "\n";

    // Nested MC: simulate from t0_opt to optT, average discounted payoff
    int numStepsOpt = 150;
    double dtOpt = (optT - t0_opt) / numStepsOpt;
    PathSimulator optSimulator(model, numStepsOpt, dtOpt);

    int numPathsOpt = 20000;  // option payoffs are more nonlinear, need more paths
    double sumPayoff = 0.0;
    for (int i = 0; i < numPathsOpt; ++i) {
        RandomGenerator rngOpt(i + 50000);
        auto path = optSimulator.simulatePath(rngOpt, stateAtT0_opt);
        double S_T = path.back()[0];
        sumPayoff += callOpt.payoff(S_T);
    }
    double mcOptPrice = std::exp(-params.r * (optT - t0_opt)) * (sumPayoff / numPathsOpt);

    std::cout << "Monte Carlo price:                        " << mcOptPrice << "\n";
    double optRelError = std::abs(mcOptPrice - analyticalOptPrice) / analyticalOptPrice * 100.0;
    std::cout << "Relative error: " << optRelError << "%\n";


    // ============================================================
// European option EE profile ¡ª uses the analytical Black-76-on-GS
// formula directly inside ExposureEngine (O(1) per time step,
// no nested Monte Carlo needed at runtime)
// ============================================================

    PathSimulator eeSimulatorOpt(model, eeSteps, forwardT / eeSteps);  // reuse optT = forwardT = 5.0 grid
    ExposureEngine engineOpt(eeSimulatorOpt, callOpt, 77);

    std::vector<double> eeCurveOpt = engineOpt.computeEE(eePaths);

    std::cout << "\n=== Expected Exposure Profile (European Call) ===\n";
    for (int j = 0; j <= eeSteps; j += 25) {
        double t = j * (forwardT / eeSteps);
        std::cout << "t = " << t << "   EE(t) = " << eeCurveOpt[j] << "\n";
    }



    // ============================================================
// Phase 2, Step 2: American option via Longstaff-Schwartz
// Sanity check: American price should be >= European price
// (early exercise is an optional right, never a disadvantage)
// ============================================================


    //int lsmSteps = 150;   // fewer steps than the EE grid ¡ª LSM cost grows with steps
    //double lsmDt = optT / lsmSteps;
    //PathSimulator lsmSimulator(model, lsmSteps, lsmDt);

    //AmericanOptionPricer amCall(lsmSimulator, optK, AmOptionType::Call, params.r);
    //double amPrice = amCall.price(20000, 99999, true);

    //std::cout << "\n=== American Option (Longstaff-Schwartz) ===\n";
    //std::cout << "American call price: " << amPrice << "\n";
    //std::cout << "European call price: " << analyticalOptPrice << "\n";
    //std::cout << (amPrice >= analyticalOptPrice ? "PASS: American >= European\n" : "FAIL: American < European (bug!)\n");




    MultiCommodityParams mcParams;
    mcParams.r = 0.03;
    mcParams.commodities = {
        {100.0, 0.05, 1.2, 0.08, 0.3, 0.1},
        {60.0,  0.04, 1.0, 0.06, 0.35, 0.12}
    };
    mcParams.correlationMatrix = {
        {1.0,  -0.5,  0.6,  -0.2},
        {-0.5,  1.0, -0.2,   0.3},
        {0.6,  -0.2,  1.0,  -0.4},
        {-0.2,  0.3, -0.4,   1.0}
    };

    MultiCommodityGSModel mcModel(mcParams);

    std::cout << "\n=== Multi-Commodity Correlation Validation ===\n";

    // Directly test the correlated draws (not full paths), same idea as Sanity Check 3
    auto L = CholeskyDecomposition::decompose(mcParams.correlationMatrix);

    int numSamplesMC = 100000;
    std::vector<std::vector<double>> wSamples(4);  // one vector per dimension

    RandomGenerator rngMC(2024);
    for (int i = 0; i < numSamplesMC; ++i) {
        auto z = rngMC.generateNormals(4);
        auto w = CholeskyDecomposition::applyToVector(L, z);
        for (int d = 0; d < 4; ++d) {
            wSamples[d].push_back(w[d]);
        }
    }

    // Compute sample correlation between dimension 0 (S1) and dimension 2 (S2), should be ~0.6
    double mean0 = std::accumulate(wSamples[0].begin(), wSamples[0].end(), 0.0) / numSamplesMC;
    double mean2 = std::accumulate(wSamples[2].begin(), wSamples[2].end(), 0.0) / numSamplesMC;
    double cov02 = 0.0;
    for (int i = 0; i < numSamplesMC; ++i) {
        cov02 += (wSamples[0][i] - mean0) * (wSamples[2][i] - mean2);
    }
    cov02 /= numSamplesMC;

    std::cout << "Sample correlation (S1, S2): " << cov02 << "  (target: 0.6)\n";
    // Ë³ÊÖÑéÖ¤ delta1 vs delta2 (index 1, 3), target: 0.3
    double mean1 = std::accumulate(wSamples[1].begin(), wSamples[1].end(), 0.0) / numSamplesMC;
    double mean3 = std::accumulate(wSamples[3].begin(), wSamples[3].end(), 0.0) / numSamplesMC;
    double cov13 = 0.0;
    for (int i = 0; i < numSamplesMC; ++i) {
        cov13 += (wSamples[1][i] - mean1) * (wSamples[3][i] - mean3);
    }
    cov13 /= numSamplesMC;
    std::cout << "Sample correlation (delta1, delta2): " << cov13 << "  (target: 0.3)\n";





    // ============================================================
// Phase 3: Netting validation (v2 ¡ª trade2 is now a SHORT position,
// so it offsets trade1 when the two positively-correlated
// commodities move together)
// ============================================================

    CommodityForward trade1(100.0, 5.0, mcParams.commodities[0].kappa, mcParams.commodities[0].alpha,
        mcParams.commodities[0].sigmaS, mcParams.commodities[0].sigmaDelta, -0.5, mcParams.r, 0);

    CommodityForward trade2Underlying(60.0, 5.0, mcParams.commodities[1].kappa, mcParams.commodities[1].alpha,
        mcParams.commodities[1].sigmaS, mcParams.commodities[1].sigmaDelta, 0.3, mcParams.r, 1);
    ShortPosition trade2(trade2Underlying);   // now a short forward on commodity 1

    Portfolio nettedPortfolio;
    nettedPortfolio.addTrade(trade1);
    nettedPortfolio.addTrade(trade2);

    int mcEeSteps = 250;
    double mcEeDt = 5.0 / mcEeSteps;
    PathSimulator mcEeSimulator(mcModel, mcEeSteps, mcEeDt);

    ExposureEngine engineTrade1(mcEeSimulator, trade1, 111);
    ExposureEngine engineTrade2(mcEeSimulator, trade2, 222);
    ExposureEngine engineNetted(mcEeSimulator, nettedPortfolio, 333);

    int mcEePaths = 3000;
    auto eeTrade1 = engineTrade1.computeEE(mcEePaths);
    auto eeTrade2 = engineTrade2.computeEE(mcEePaths);
    auto eeNetted = engineNetted.computeEE(mcEePaths);

    std::cout << "\n=== Netting Effect Validation (long/short offsetting) ===\n";
    for (int j = 0; j <= mcEeSteps; j += 50) {
        double t = j * mcEeDt;
        double standaloneSum = eeTrade1[j] + eeTrade2[j];
        std::cout << "t=" << t
            << "  EE(trade1)=" << eeTrade1[j]
            << "  EE(trade2 short)=" << eeTrade2[j]
            << "  Sum=" << standaloneSum
            << "  EE(netted)=" << eeNetted[j] << "\n";
    }

    // ============================================================
// Phase 4: CVA / FVA on the single-commodity forward EE profile
// (reusing eeCurve computed earlier in Phase 1)
// ============================================================

    CreditCurve counterpartyCurve(0.02);   // flat hazard rate, ~2% annualized (illustrative)
    double recoveryRate = 0.4;              // standard market assumption for senior unsecured
    XVACalculator xva(counterpartyCurve, recoveryRate, params.r);

    double cva = xva.computeCVA(eeCurve, forwardT / eeSteps);
    double fundingSpread = 0.015;           // 150bps illustrative funding spread
    double fva = xva.computeFVA(eeCurve, forwardT / eeSteps, fundingSpread);

    std::cout << "\n=== XVA (single forward) ===\n";
    std::cout << "CVA: " << cva << "\n";
    std::cout << "FVA: " << fva << "\n";

    // Sensitivity check: doubling hazard rate should roughly double CVA (for small hazard rates)
    CreditCurve riskierCounterparty(0.04);
    XVACalculator xvaRiskier(riskierCounterparty, recoveryRate, params.r);
    double cvaRiskier = xvaRiskier.computeCVA(eeCurve, forwardT / eeSteps);
    std::cout << "CVA with 2x hazard rate: " << cvaRiskier << " (should be roughly ~2x original CVA)\n";


    // ============================================================
// Phase 4: Wrong-Way Risk validation
// Three-way comparison:
//   1) Analytical CVA (from XVACalculator, assumes independence)
//   2) Monte Carlo CVA with gamma=0 (should match #1, cross-check)
//   3) Monte Carlo CVA with gamma>0 (WWR: should be HIGHER than #1/#2)
// ============================================================

    double lambda0 = 0.02;
    double gammaWWR = -3.0;   // negative sign because mainForward is a LONG position
    // (exposure rises with S, so hazard rate must also rise with S for true WWR)


    WrongWayRiskEngine engineNoWWR(eeSimulator, mainForward, lambda0, 0.0, recoveryRate, params.r);
    WrongWayRiskEngine engineWithWWR(eeSimulator, mainForward, lambda0, gammaWWR, recoveryRate, params.r);


    int wwrPaths = 20000;
    double cvaMcNoWWR = engineNoWWR.computeCVA(wwrPaths, 5000, 0);
    double cvaMcWithWWR = engineWithWWR.computeCVA(wwrPaths, 5000, 0);

    std::cout << "\n=== Wrong-Way Risk Validation ===\n";
    std::cout << "Analytical CVA (independence assumption): " << cva << "\n";
    std::cout << "Monte Carlo CVA (gamma=0, cross-check):    " << cvaMcNoWWR << "\n";
    std::cout << "Monte Carlo CVA (gamma=3, with WWR):       " << cvaMcWithWWR << "\n";
    std::cout << "WWR uplift: " << ((cvaMcWithWWR / cvaMcNoWWR - 1.0) * 100.0) << "%\n";


    std::cout << "\n=== WWR Sensitivity to gamma ===\n";
    for (double g : {0.0, -1.0, -2.0, -3.0, -5.0}) {
        WrongWayRiskEngine engineG(eeSimulator, mainForward, lambda0, g, recoveryRate, params.r);
        double cvaG = engineG.computeCVA(wwrPaths, 5000, 0);
        std::cout << "gamma=" << g << "  CVA=" << cvaG << "\n";
    }
    return 0;
}