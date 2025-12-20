#include <random>

#ifndef NOUTURN_HPP_UNIFORMGENERATOR_H
#define NOUTURN_HPP_UNIFORMGENERATOR_H

// Uniform random number generator used for NUTS tree construction.
// This is a wrapper around std::mt19937 and always produces
// genuine uniform random variable in (0, 1).

class UniformGenerator {
public:
  explicit UniformGenerator(int seed)
    : generator(seed),
      distribution(0.0, 1.0) {}
  
  inline double getUniform() {
    return distribution(generator);
  }
  
private:
  std::mt19937 generator;
  std::uniform_real_distribution<double> distribution;
};

#endif // NOUTURN_HPP_UNIFORMGENERATOR_H
