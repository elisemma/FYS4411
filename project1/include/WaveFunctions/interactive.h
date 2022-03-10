#pragma once

#include "wavefunction.h"

class Interactive : public WaveFunction {
public:
    Interactive(System* system, double alpha_, double beta_, double a_);
    double evaluate(std::vector<Particle*> particles);
    double computeDoubleDerivativeAnalytical(std::vector<Particle*> particles);
    std::vector<double> computeQuantumForceAnalytical(Particle* particles);
    double computeDerivative(std::vector<Particle*> particles);
private:
    double r_beta(std::vector<Particle*> particles, int k);
    double r_beta_squared(std::vector<Particle*> particles, int k);
    double getDistance(std::vector<Particle*> particles, int i, int j);
    double computePhiLaplacian();
    double computeDistance(int k, int j, std::vector<Particle*> particles);
    double computeU(double r_kj);
    double computeUPrime(double r_kj);
    double computeUDoublePrime(double r_kj);
    double computeLaplacianK(std::vector<Particle*> particles, std::size_t k);
    double a;
    double beta;
};
