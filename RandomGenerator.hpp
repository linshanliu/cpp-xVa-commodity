#ifndef RANDOMGENERATOR_HPP
#define RANDOMGENERATOR_HPP
#include <vector>
#include <random>

class RandomGenerator {
public:
    explicit RandomGenerator(unsigned int seed);

    std::vector<double> generateNormals(int n);

private:
    std::mt19937 engine_;
    std::normal_distribution<double> dist_;
};
#endif