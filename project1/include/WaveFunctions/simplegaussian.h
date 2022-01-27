#pragma once

#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(System* system, double alpha);
    double evaluate(std::vector<Particle*> particles);
    double computeDoubleDerivative(std::vector<Particle*> particles);
};
