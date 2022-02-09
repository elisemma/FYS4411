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

double SimpleGaussian::computeDoubleDerivative(vector<Particle*> particles, bool analytical) {
    if (analytical) {
      return computeDoubleDerivativeAnalytical(particles);
    } else {
      return computeDoubleDerivativeNumerical(particles);
    }
}

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
    double double_derivative = 0, step_length = m_system->getStepLength();

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
