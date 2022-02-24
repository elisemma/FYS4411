#include "WaveFunctions/wavefunction.h"
#include "project1/system.h"

using namespace std;

WaveFunction::WaveFunction(System* system) {
    m_system = system;
}

void WaveFunction::move_particles(std::vector<Particle*> particles, double step_length, int dim) {
    for (auto particle : particles) {
        particle->adjustPosition(step_length, dim);
    }
}



double WaveFunction::computeDoubleDerivativeNumerical(vector<Particle*> particles) {
    //TODO: Set the step length in a better and more general way
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

