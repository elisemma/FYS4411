#pragma once

#include "wavefunction.h"

class Interactive : public WaveFunction {
public:
    Interactive(System* system, double alpha_, double gamma_squared_, double a_);
    double evaluate(std::vector<Particle*> particles);
    double computeDoubleDerivativeAnalytical(std::vector<Particle*> particles);
    double computeDoubleDerivativeNumerical(std::vector<Particle*> particles);
    std::vector<double> computeQuantumForceAnalytical(Particle* particles);
    double computeDerivative(std::vector<class Particle*> particles);
private:
    double a;
    double gamma_squared;
};
