#include "WaveFunctions/interactive.h"
#include "WaveFunctions/wavefunction.h"
#include "project1/system.h"
#include "project1/particle.h"
#include <cmath>
#include <cassert>

using namespace std;

Interactive::Interactive(System* system, double alpha_, double beta_, double a_) : WaveFunction(system) {
    assert(alpha_ >= 0);
    assert(beta_ >= 0);

    this->alpha = alpha_;
    this->beta = beta_;
    this->a = a_;
}

double Interactive::evaluate(vector<Particle*> particles) {
    double r_squared = 0;

    for (auto particle : particles) {
        vector<double> particle_pos = particle->getPosition();
        for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
            if (dim == 2) {
                r_squared += particle_pos[2] * particle_pos[2] * beta;
            } else {
                r_squared += particle_pos[dim] * particle_pos[dim];
            }
        }
    }

    double g_product = exp(-alpha * r_squared);

    int number_of_particles = particles.size();
    double f_product = 1;
    for (int j = 0; j < number_of_particles - 1; j++) {
        for (int k = j + 1; k < number_of_particles; k++) {
            double distance = getDistance(particles, j, k);

            // TODO: Is this one named Mr. Jastrow?
            if (distance > a) {
                f_product *= 1 - (a / distance);
            } else {
                f_product *= 0;
                break;
            }
        }
    }

    return g_product * f_product;
}

// TODO: Calculate all the distances at once for improved speed
double Interactive::getDistance(std::vector<class Particle*> particles, int i, int j) {
    vector<double> pos_i = particles[i]->getPosition();
    vector<double> pos_j = particles[j]->getPosition();

    double distance = 0;
    for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
        distance += abs(pos_i[dim] - pos_j[dim]);
    }

    return distance;
}

// TODO: THESE ARE LIKE BIG TIME WRONG, PLEASE DONT EVER FOREVER USE THEM!!
double Interactive::computeDoubleDerivativeAnalytical(vector<Particle*> particles) {
    int numberOfDimensions = m_system->getNumberOfDimensions();
    int numberOfParticles = m_system->getNumberOfParticles();
    double r_squared = m_system->calculate_r_squared(particles);
    return -2*numberOfParticles*numberOfDimensions*alpha + 4*pow(alpha,2)*r_squared;
}

vector<double> Interactive::computeQuantumForceAnalytical(Particle* particle) {
    vector<double> force;
    for (auto dimension_pos : particle->getPosition()) {
        force.push_back(-4*alpha*dimension_pos);
    }
    return force;
}


double Interactive::computeDerivative(vector<Particle*> particles) {
    // TODO: Should we divide by the wave function giving -alpha*r_squared?
    // And do we want the -2? (mostly because I have seen them a few times before, and not particularly because I think they would fit)

    double r_squared = m_system->calculate_r_squared(particles);
    return -r_squared;
    // return -2*alpha * r_squared;
    // double minus_alpha_r_squared = - alpha * m_system->calculate_r_squared(particles);
    //
    // return minus_alpha_r_squared * exp(minus_alpha_r_squared);
}
