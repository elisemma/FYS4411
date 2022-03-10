#pragma once

#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(System* system, double alpha);
    double evaluate(std::vector<Particle*> particles);
    // double computeDoubleDerivative(std::vector<Particle*> particles, bool analytical);
    double computeDoubleDerivativeAnalytical(std::vector<Particle*> particles);
    // double computeDoubleDerivativeNumerical(std::vector<Particle*> particles);
    std::vector<double> computeQuantumForceAnalytical(Particle* particle);
    double computeDerivative(std::vector<class Particle*> particles);
// private:
//     void move_particles(std::vector<Particle*> particles, double step_length, int dim);
};
