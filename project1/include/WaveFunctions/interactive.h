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
    double computePhiLaplacian();
    double computeDistance(int k, int j);
    double computeUPrime(int k, int i);
    double computeNaplaK(std::vector<Particle*> particles);
    double computeLaplacianK(std::vector<Particle*> particles, int k);
    double computeUDoublePrime(int k, int i);
    // double computePhiNabla();
    double a;
    double beta;
};
