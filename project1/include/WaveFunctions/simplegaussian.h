#pragma once

#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(System* system, double alpha);
    double evaluate(std::vector<Particle*> particles);
    double computeDoubleDerivative(std::vector<Particle*> particles);
private:
    double calculate_r_squared(std::vector<Particle*> particles);
};
