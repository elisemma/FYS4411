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
    // double r_squared = m_system->calculate_r_squared(particles);
    double r_squared = m_system->getRSquared();
    // assert(abs(m_system->getRSquared() - m_system->calculate_r_squared(particles)) < 1e-13);
    // if ((m_system->getRSquared() != m_system->calculate_r_squared(particles))) {
    //     exit(1);
    // }
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
    // double r_squared = m_system->calculate_r_squared(particles);
    double r_squared = m_system->getRSquared();

    return -2*numberOfParticles*numberOfDimensions*alpha + 4*pow(alpha,2)*r_squared;
}

vector<double> SimpleGaussian::computeQuantumForceAnalytical(Particle* particle) {
    vector<double> force;
    for (auto dimension_pos : particle->getPosition()) {
        force.push_back(-4*alpha*dimension_pos);
    }
    return force;
}


double SimpleGaussian::computeDerivative(vector<Particle*> particles) {
    // And do we want the -2? (mostly because I have seen them a few times before, and not particularly because I think they would fit)

    // double r_squared = m_system->calculate_r_squared(particles);
    double r_squared = m_system->getRSquared();
    return -r_squared;
    // return -2*alpha * r_squared;
    // double minus_alpha_r_squared = - alpha * m_system->calculate_r_squared(particles);
    //
    // return minus_alpha_r_squared * exp(minus_alpha_r_squared);
}
