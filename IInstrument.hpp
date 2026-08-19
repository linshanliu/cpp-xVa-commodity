#ifndef IINSTRUMENT_HPP
#define IINSTRUMENT_HPP
#include <vector>

class IInstrument {
public:
    virtual ~IInstrument() = default;

    // 给定当前时刻t、当前市场状态[S_t, delta_t], 计算这笔交易此刻的MtM
    virtual double markToMarket(double t, const std::vector<double>& state) const = 0;

    // 这笔交易的到期时间,ExposureEngine要知道过了到期日就不用再算了
    virtual double maturity() const = 0;
};
#endif