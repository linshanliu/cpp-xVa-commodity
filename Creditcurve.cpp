#include "CreditCurve.hpp"
#include <cmath>

CreditCurve::CreditCurve(double hazardRate) : hazardRate_(hazardRate) {}

double CreditCurve::survivalProbability(double t) const {
    return std::exp(-hazardRate_ * t);
}

double CreditCurve::defaultProbability(double t) const {
    return 1.0 - survivalProbability(t);
}