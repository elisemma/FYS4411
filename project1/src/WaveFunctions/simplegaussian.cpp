#include "WaveFunctions/simplegaussian.h"
#include <cmath>
#include <cassert>
#include "WaveFunctions/wavefunction.h"
#include "project1/system.h"
#include "project1/particle.h"

SimpleGaussian::SimpleGaussian(System* system, double alpha_) : WaveFunction(system) {
    assert(alpha_ >= 0);
    this->alpha = alpha_;
}

// TODO: Should this value be saved/cached?
double SimpleGaussian::calculate_r_squared(std::vector<Particle*> particles) {
    double r_squared = 0;
    for (Particle* particle : particles) {
        for (double pos_i : particle->getPosition()) {
            r_squared += pow(pos_i, 2);
        }
    }
    return r_squared;
}

double SimpleGaussian::evaluate(std::vector<Particle*> particles) {
    double r_squared = calculate_r_squared(particles);
    return exp(-alpha * r_squared);
}


double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */
    double r_squared = calculate_r_squared(particles);
    return 0; // the analytical expression goes here
}
