// TODO: Remove this iostream
#include <iostream>
#include "WaveFunctions/simplegaussian.h"
#include <cmath>
#include <cassert>
#include "WaveFunctions/wavefunction.h"
#include "project1/system.h"
#include "project1/particle.h"

using namespace std;

SimpleGaussian::SimpleGaussian(System* system, double alpha_) : WaveFunction(system) {
    assert(alpha_ >= 0);
    this->alpha = alpha_;
}

double SimpleGaussian::evaluate(vector<Particle*> particles) {
    double r_squared = m_system->calculate_r_squared(particles);
    return exp(-alpha * r_squared);
}

// double SimpleGaussian::computeDoubleDerivative(vector<Particle*> particles, bool analytical) {
//     if (analytical) {
//       return computeDoubleDerivativeAnalytical(particles);
//     } else {
//       return computeDoubleDerivativeNumerical(particles);
//     }
// }

double SimpleGaussian::computeDoubleDerivativeAnalytical(vector<Particle*> particles) {

    int numberOfDimensions = m_system->getNumberOfDimensions();
    int numberOfParticles = m_system->getNumberOfParticles();
    double r_squared = m_system->calculate_r_squared(particles);

    return -2*numberOfParticles*numberOfDimensions*alpha + 4*pow(alpha,2)*r_squared;
}

void move_particles(vector<Particle*> particles, double step_length, int dim) {
    for (auto particle : particles) {
        particle->adjustPosition(step_length, dim);
    }
}

double SimpleGaussian::computeDoubleDerivativeNumerical(vector<Particle*> particles) {
  //TODO: Ta med step length i rapport, og numerisk feil (mindre og tregere enn analytisk)
  //TODO: Test for ulike step length og undersøk CPU tid
    double double_derivative = 0, step_length = 1e-6;

    for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
        move_particles(particles, -step_length, dim);
        double_derivative += m_system->getWaveFunction()->evaluate(particles);
        move_particles(particles, 2*step_length, dim);
        double_derivative += m_system->getWaveFunction()->evaluate(particles);
        move_particles(particles, -step_length, dim);
    }
    double_derivative -= 2*m_system->getNumberOfDimensions()*m_system->getWaveFunction()->evaluate(particles);

    double evaluated_wave = m_system->getWaveFunction()->evaluate(particles);
    double_derivative /= (step_length * step_length * evaluated_wave);

    return double_derivative;
}

vector<double> SimpleGaussian::computeQuantumForceAnalytical(Particle* particle) {
    vector<double> force;
    for (auto dimension_pos : particle->getPosition()) {
        force.push_back(-4*alpha*dimension_pos);
    }
    return force;
}


double SimpleGaussian::computeDerivative(vector<Particle*> particles) {
    // TODO: Should we divide by the wave function giving -alpha*r_squared?
    // And do we want the -2? (mostly because I have seen them a few times before, and not particularly because I think they would fit)

    double r_squared = m_system->calculate_r_squared(particles);
    return -r_squared;
    // return -2*alpha * r_squared;
    // double minus_alpha_r_squared = - alpha * m_system->calculate_r_squared(particles);
    //
    // return minus_alpha_r_squared * exp(minus_alpha_r_squared);
}
