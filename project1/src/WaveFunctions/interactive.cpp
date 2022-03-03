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
    double f_product = 1; // f_product seems to be the jastrow factor, nice to meet you!
    for (int j = 0; j < number_of_particles - 1; j++) {
        for (int k = j + 1; k < number_of_particles; k++) {
            double distance = m_system->getDistance(j, k);

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

// TODO: THESE ARE LIKE BIG TIME WRONG, PLEASE DONT EVER FOREVER USE THEM!!
double Interactive::computeDoubleDerivativeAnalytical(vector<Particle*> particles) {
    double laplacian = 0;

    for (int k = 0; k < m_system->getNumberOfParticles(); k++) {
        laplacian += computeLaplacianK(particles, k);
    }

    return laplacian;
}

double Interactive::computeLaplacianK(vector<Particle*> particles, int k) {
    // Term 1
    double term1 = 0;
    double r_squared_with_beta = 0;
    for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
        r_squared_with_beta += particles[k]->getPosition()[dim] * particles[k]->getPosition()[dim];
        if (dim == 2) {
            r_squared_with_beta *= beta * beta;
        }
    }
    term1 = 2 * alpha * (2 * alpha * r_squared_with_beta - (2 + beta));
    
    // Term 2
    double term2 = 0;
    for (int i = 0; i < m_system->getNumberOfParticles(); i++) {
        if (i == k) continue;

        double r_ki = computeDistance(k, i);

        for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
            double nabla_phi = - 4 * alpha * particles[k]->getPosition()[dim];
            if (dim == 2) nabla_phi *= beta;

            double r_k_r_i = particles[k]->getPosition()[dim] - particles[i]->getPosition()[dim];

            term2 += - 2 * nabla_phi * (r_k_r_i / r_ki)  * computeUPrime(i, k);
        }
    }

    // Term 3
    double term3 = 0;

    for (int i = 0; i < m_system->getNumberOfParticles(); i++) {
        if (i == k) continue;

        for (int j = 0; j < m_system->getNumberOfParticles(); j++) {
            if (j == k) continue;
            
            double r_ki = computeDistance(k, j);
            double r_kj = computeDistance(k, i);

            for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
                double r_k_r_i = particles[k]->getPosition()[dim] - particles[i]->getPosition()[dim];
                double r_k_r_j = particles[k]->getPosition()[dim] - particles[j]->getPosition()[dim];

                term3 += (r_k_r_i * r_k_r_j) / (r_ki * r_kj) * computeUPrime(k, i) * computeUPrime(k, j);
            }
        }
    }

    // Term 4
    double term4 = 0;
    for (int i = 0; i < m_system->getNumberOfParticles(); i++) {
        if (i == k) continue;

        double r_ki = computeDistance(k, i);

        // TODO: i term 4 flytter vi på leddene
        term4 += 2 / r_ki * computeUPrime(k, i) + computeUDoublePrime(k, i);
    }

    return term1 + term2 + term3 + term4;
}

double Interactive::computeUPrime(int k, int i) {
    double r_ki = computeDistance(k, i);
    return a / (r_ki * (r_ki - a));
}

double Interactive::computeUDoublePrime(int i, int k) {
    double r_ki = computeDistance(k, i);
    return a * (a - 2*r_ki) / ((r_ki * r_ki) * (r_ki-a) * (r_ki-a));
}

double Interactive::computeDistance(int k, int j) {
    double distance = 0;
    for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
        double single_distance = m_system->getParticles()[k]->getPosition()[dim] - m_system->getParticles()[j]->getPosition()[dim];
        distance += single_distance * single_distance;
    }
    return sqrt(distance);
}

vector<double> Interactive::computeQuantumForceAnalytical(Particle* particle) {
    vector<double> force;
    for (auto dimension_pos : particle->getPosition()) {
        force.push_back(-4*alpha*dimension_pos);
    }
    return force;
}


double Interactive::computeDerivative(vector<Particle*> particles) {
    // And do we want the -2? (mostly because I have seen them a few times before, and not particularly because I think they would fit)

    double r_squared = m_system->calculate_r_squared(particles);
    return -r_squared;
    // return -2*alpha * r_squared;
    // double minus_alpha_r_squared = - alpha * m_system->calculate_r_squared(particles);
    //
    // return minus_alpha_r_squared * exp(minus_alpha_r_squared);
}
