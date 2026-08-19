#ifndef COMMODITYFORWARD_HPP
#define COMMODITYFORWARD_HPP
#include "IInstrument.hpp"

class CommodityForward : public IInstrument {
public:
    // K: 约定交割价, T: 到期时间
    // kappa, alpha, sigmaS, sigmaDelta, rho, r: 定价需要的模型参数
    CommodityForward(double K, double T, double kappa, double alpha,
        double sigmaS, double sigmaDelta, double rho, double r);

    double markToMarket(double t, const std::vector<double>& state) const override;
    double maturity() const override { return T_; }

    // 暴露出来方便单独测试/对照MonteCarlo
    double forwardPrice(double t, const std::vector<double>& state) const;

private:
    double K_;
    double T_;
    double kappa_;
    double alpha_;
    double sigmaS_;
    double sigmaDelta_;
    double rho_;
    double r_;
};
#endif