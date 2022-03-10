#include <iostream>
#include "WaveFunctions/interactive.h"
#include "WaveFunctions/wavefunction.h"
#include "project1/system.h"
#include "project1/particle.h"
#include <cmath>
#include <cassert>

using namespace std;

#define px particle->getPosition()[0]
#define py particle->getPosition()[1]
#define pz particle->getPosition()[2]
#define rx(i) particles[i]->getPosition()[0]
#define ry(i) particles[i]->getPosition()[1]
#define rz(i) particles[i]->getPosition()[2]

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
            // TODO: Make sure this should not include beta
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

double Interactive::r_beta(vector<Particle*> particles, int k) {
    return rx(k) + ry(k) + beta*rz(k);
}

double Interactive::r_beta_squared(vector<Particle*> particles, int k) {
    return rx(k)*rx(k) + ry(k)*ry(k) + beta*rz(k)*rz(k);
}

// TODO: THESE ARE LIKE BIG TIME WRONG, PLEASE DONT EVER FOREVER USE THEM!!
double Interactive::computeDoubleDerivativeAnalytical(vector<Particle*> particles) {
    // TODO: Remove this for the final version
    // int N = m_system->getNumberOfParticles();
    //
    // double laplacian = 0, term_2_sum[3] = {0, 0, 0}, term3_sum = 0, term_4_sum = 0;
    //
    // double *r_kj = (double*) malloc(N * N * sizeof(double));
    // double *r_k_r_j = (double*) malloc(3* N * N * sizeof(double));
    // double *u_prime_kj = (double*) malloc(N * N * sizeof(double));
    // double *term_2_sum_terms = (double*) malloc(3 * N * N * sizeof(double));
    // double *term_4_terms = (double*) malloc(N * N * sizeof(double));
    //
    // for (int k = 0; k < N; k++) {
    //     // TODO: Make sure the term for k=N (N-1) is also included
    //     for (int j = k+1; j < N; j++) {
    //         r_kj[k*N + j] = computeDistance(k, j, particles);
    //
    //         r_k_r_j[3*(k*N + j) + 0] = particles[k]->getPosition()[0] - particles[j]->getPosition()[0];
    //         r_k_r_j[3*(k*N + j) + 1] = particles[k]->getPosition()[1] - particles[j]->getPosition()[1];
    //         r_k_r_j[3*(k*N + j) + 2] = particles[k]->getPosition()[2] - particles[j]->getPosition()[2];
    //
    //         u_prime_kj[k*N + j] = computeUPrime(r_kj[k*N + j]);
    //
    //         // Term 2
    //         // TODO: Include the term for k=N (N-1)
    //         term_2_sum_terms[3*k + 0] = r_k_r_j[3*(k*N + j) + 0] / abs(r_kj[k*N + j]) * u_prime_kj[k*N + j];
    //         term_2_sum_terms[3*k + 1] = r_k_r_j[3*(k*N + j) + 1] / abs(r_kj[k*N + j]) * u_prime_kj[k*N + j];
    //         term_2_sum_terms[3*k + 2] = r_k_r_j[3*(k*N + j) + 2] / abs(r_kj[k*N + j]) * u_prime_kj[k*N + j];
    //
    //         term_2_sum[0] += term_2_sum_terms[3*k + 0];
    //         term_2_sum[1] += term_2_sum_terms[3*k + 1];
    //         term_2_sum[2] += term_2_sum_terms[3*k + 2];
    //
    //         // Term 4
    //         term_4_terms[k*N + j] += computeUDoublePrime(r_kj[k*N + j]) + 2 / r_kj[k*N + j] * u_prime_kj[k*N + j];
    //         term_4_sum += term_4_terms[k*N + j];
    //     }
    // }
    //
    // for (int j = 0; j < N; j++) {
    //     for (int i = 0; i < N; i++) {
    //
    //     }
    // }
    //
    // for (size_t k = 0; k < N; k++) {
    //     // Term 1
    //     // TODO: Are these squared?
    //     double r_k_beta_squared = r_beta_squared(particles, k);
    //     laplacian += 2 * alpha * (2 * alpha * r_k_beta_squared - (2 + beta));
    // }

    double laplacian = 0;

    // double *r_kj = (double*) malloc(particles.size() * particles.size() * sizeof(double));
    // double *u_kj_prime = (double*) malloc(particles.size() * particles.size() * sizeof(double));

    // TODO: Can we get away with only filling the upper triangle?
    // for (size_t k = 0; k < particles.size(); k++) {
    //     r_kj[k * N+k] = 0;
    //     for (size_t j = k+1; j < particles.size(); j++) {
    //         r_kj[k*N + j] = computeDistance(k, j, particles);
    //         r_kj[j*N + k] = r_kj[k*N + j];
    //
    //
    //         u_kj_prime[k*N + j] = computeUPrime(k, j, r_kj[k*N + j]);
    //     }
    // }

    for (size_t k = 0; k < particles.size(); k++) {
        laplacian += computeLaplacianK(particles, k);
    }

    return laplacian;
}

// double Interactive::computeLaplacianK(vector<Particle*> particles, int k) {
double Interactive::computeLaplacianK(vector<Particle*> particles, size_t k) {
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
    for (size_t i = 0; i < particles.size(); i++) {
        if (i == k) continue;

        double r_ki = computeDistance(k, i, particles);

        for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
            double nabla_phi = - 2 * 2 * alpha * particles[k]->getPosition()[dim];
            if (dim == 2) nabla_phi *= beta;

            double r_k_r_i = particles[k]->getPosition()[dim] - particles[i]->getPosition()[dim];

            term2 += - 2 * nabla_phi * (r_k_r_i / r_ki)  * computeUPrime(computeDistance(k, i, particles));
        }
    }

    // Term 3
    double term3 = 0;

    for (size_t i = 0; i < particles.size(); i++) {
        if (i == k) continue;

        for (size_t j = 0; j < particles.size(); j++) {
            if (j == k) continue;
            
            double r_ki = computeDistance(k, j, particles);
            double r_kj = computeDistance(k, i, particles);

            for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
                double r_k_r_i = particles[k]->getPosition()[dim] - particles[i]->getPosition()[dim];
                double r_k_r_j = particles[k]->getPosition()[dim] - particles[j]->getPosition()[dim];

                term3 += (r_k_r_i * r_k_r_j) / (r_ki * r_kj) * computeUPrime(computeDistance(k, i, particles)) * computeUPrime(computeDistance(k, j, particles));
            }
        }
    }

    // Term 4
    double term4 = 0;
    for (size_t i = 0; i < particles.size(); i++) {
        if (i == k) continue;

        double r_ki = computeDistance(k, i, particles);

        // TODO: i term 4 flytter vi på leddene
        term4 += 2 / r_ki * computeUPrime(computeDistance(k, i, particles)) + computeUDoublePrime(computeDistance(k, i, particles));
    }

    return term1 + term2 + term3 + term4;
}

double Interactive::computeU(double r_kj) {
    // TODO: Do I need to check if it is smaller than a?
    if (r_kj <= a) {
        return 0;
    } else {
        return 1 - a/r_kj;
    }
}


double Interactive::computeUPrime(double r_kj) {
    // double r_ki = computeDistance(k, i, particles);
    return a / (r_kj * (r_kj - a));
}

double Interactive::computeUDoublePrime(double r_kj) {
    return a * (a - 2*r_kj) / ((r_kj * r_kj) * (r_kj-a) * (r_kj-a));
}

double Interactive::computeDistance(int k, int j, vector<Particle*> particles) {
    // double distance = 0;
    // for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
    //     double single_distance = particles[k]->getPosition()[dim] - particles[j]->getPosition()[dim];
    //     distance += single_distance * single_distance;
    // }
    return sqrt((rx(k)-rx(j))*(rx(k)-rx(j)) + (ry(k)-ry(j))*(ry(k)-ry(j)) + (rz(k)-rz(j))*(rz(k)-rz(j)));
}


vector<double> Interactive::computeQuantumForceAnalytical(Particle* particle) {
    vector<double> force = {
        -4 * alpha *        particle->getPosition()[0],
        -4 * alpha *        particle->getPosition()[1],
        -4 * alpha * beta * particle->getPosition()[2],
    };

    return force;
}


double Interactive::computeDerivative(vector<Particle*> particles) {
    double nabla = 0;

    double r_squared = 0;

    for (auto particle : particles) {
        vector<double> particle_pos = particle->getPosition();
        r_squared += px*px + py*py + beta*pz*pz;
    }

    nabla = - r_squared * evaluate(particles);

    return nabla;
}
