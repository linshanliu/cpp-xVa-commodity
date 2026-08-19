#ifndef GIBSONSCHWARTZPARAMS_HPP
#define GIBSONSCHWARTZPARAMS_HPP




struct GibsonSchwartzParams {
    double S0;          // underlying 
    double delta0;      // convenience yield at 0
    double kappa;       // mean reverting speed
    double alpha;       // long term average convenience yield
    double sigmaS;      // underlying vol
    double sigmaDelta;  // convenience vol
    double rho;         // covariance for Brownian motion
    double r;           // risk free rates
};


#endif // !GIBSONSCHWARTZPARAMS_HPP