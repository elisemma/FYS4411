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
    double r_squared = m_system->calculate_r_squared(particles);
    // cout << "r2: " << r_squared << endl;
    return 2*alpha*exp(-alpha*r_squared)*(2*alpha*r_squared - 1);
}

// double SimpleGaussian::computeDoubleDerivativeNumerical(vector<class Particle*> particles){
//   return 0;
// }

void move_particles(vector<Particle*> particles, double step_length, int dim) {
    for (auto particle : particles) {
        particle->adjustPosition(step_length, dim);
    }
}

double SimpleGaussian::computeDoubleDerivativeNumerical(vector<Particle*> particles) {
    double double_derivative = 0, step_length = m_system->getStepLength();

    for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
        // TODO: This is the same each time, so maybe move outside of the loop if working?
        double_derivative -= 2*m_system->getWaveFunction()->evaluate(particles);

        move_particles(particles, -step_length, dim);

        double_derivative += m_system->getWaveFunction()->evaluate(particles);

        move_particles(particles, 2*step_length, dim);

        double_derivative += m_system->getWaveFunction()->evaluate(particles);

        move_particles(particles, -step_length, dim);
    }

    double evaluated_wave = m_system->getWaveFunction()->evaluate(particles);
    double_derivative /= (step_length * step_length * evaluated_wave);

    return double_derivative;
}
