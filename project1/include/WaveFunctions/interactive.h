#pragma once

#include "wavefunction.h"

class Interactive : public WaveFunction {
public:
    Interactive(System* system, double alpha_, double beta_, double a_);
    double evaluate(std::vector<Particle*> particles);
    double computeDoubleDerivativeAnalytical(std::vector<Particle*> particles);
    // double computeDoubleDerivativeNumerical(std::vector<Particle*> particles);
    std::vector<double> computeQuantumForceAnalytical(Particle* particles);
    double computeDerivative(std::vector<class Particle*> particles);
private:
    double getDistance(std::vector<class Particle*> particles, int i, int j);
    double a;
    double beta;
};
