#include "RandomGenerator.hpp"

RandomGenerator::RandomGenerator(unsigned int seed) : engine_(seed), dist_(0.0, 1.0) {}

std::vector<double> RandomGenerator::generateNormals(int n) {
    std::vector<double> result;
    result.reserve(n);
    for (int i = 0; i < n; ++i) {
        result.push_back(dist_(engine_));
    }
    return result;
}
