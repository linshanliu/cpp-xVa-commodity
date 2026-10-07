#ifndef IINSTRUMENT_HPP
#define IINSTRUMENT_HPP
#include <vector>

class IInstrument {
public:
    virtual ~IInstrument() = default;

	// given the current time t and the current state of the model, return the mark-to-market value of this instrument
    virtual double markToMarket(double t, const std::vector<double>& state) const = 0;


	// getter for the maturity of the instrument, ExposureEngine needs to know when the instrument has matured
    // so it can stop calculating exposure after that date
    virtual double maturity() const = 0;
};
#endif