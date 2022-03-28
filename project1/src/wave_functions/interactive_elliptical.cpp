#include <iostream>
#include <cmath>
#include <cassert>
#include "project1/wave_functions/interactive_elliptical.h"
#include "project1/constants.h"

InteractiveElliptical::InteractiveElliptical(int num_particles, int dimensions, double alpha, double a): num_particles(num_particles), alpha(alpha), a(a) {
    assert(dimensions == 3 && "Elliptical wave function only works in 3 dimensions");
}

double InteractiveElliptical::evaluate(double *particles) {
    double r_squared = compute_r_squared(particles);

    double g_product = exp(-alpha * r_squared), f_product = 1;

    for (int particle_k = 0; particle_k < num_particles; particle_k++) {
        for (int particle_j = particle_k + 1; particle_j < num_particles; particle_j++) {
            double distance = compute_distance(particles, particle_k, particle_j);

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


double InteractiveElliptical::compute_laplacian(double *particles) {
    double term1 = 0, term2 = 0, term3 = 0, term4 = 0;
    double r_kx, r_ky, r_kz, r_jx, r_jy, r_jz;

    double r_kj, u_gradient_correlations[3] = {0, 0, 0};

    for (int k = 0; k < num_particles; k++) {
        r_kx = particles[3 * k + 0];
        r_ky = particles[3 * k + 1];
        r_kz = particles[3 * k + 2];

        // Compute term 1
        double r_squared_with_beta = r_kx*r_kx + r_ky*r_ky + constants::beta*constants::beta*r_kz*r_kz;
        term1 += 2 * alpha * (2 * alpha * r_squared_with_beta - (2 + constants::beta));

        u_gradient_correlations[0] = 0;
        u_gradient_correlations[1] = 0;
        u_gradient_correlations[2] = 0;
        for (int j = 0; j < num_particles; j++) {
            r_jx = particles[3 * j + 0];
            r_jy = particles[3 * j + 1];
            r_jz = particles[3 * j + 2];

            if (k == j) continue;
            r_kj = compute_distance(particles, k, j);

            double u_prime_over_r_kj = compute_u_prime(r_kj) / r_kj;

            u_gradient_correlations[0] += (r_kx - r_jx) * u_prime_over_r_kj;
            u_gradient_correlations[1] += (r_ky - r_jy) * u_prime_over_r_kj;
            u_gradient_correlations[2] += (r_kz - r_jz) * u_prime_over_r_kj;

            // Compute term 4
            term4 += 2 * u_prime_over_r_kj;
            term4 += compute_u_double_prime(r_kj);
        }

        // Compute term 2
        term2 += 2 * (-2) * alpha *                   r_kx * u_gradient_correlations[0];
        term2 += 2 * (-2) * alpha *                   r_ky * u_gradient_correlations[1];
        term2 += 2 * (-2) * alpha * constants::beta * r_kz * u_gradient_correlations[2];

        // Compute term 3
        term3 += u_gradient_correlations[0] * u_gradient_correlations[0];
        term3 += u_gradient_correlations[1] * u_gradient_correlations[1];
        term3 += u_gradient_correlations[2] * u_gradient_correlations[2];
    }

    return term1 + term2 + term3 + term4;
}

double InteractiveElliptical::compute_r_squared(double *particles) {
    double r_squared = 0;

    for (int particle = 0; particle < num_particles; particle++) {
        r_squared +=
            particles[3 * particle + 0] * particles[3 * particle + 0]
                +                   particles[3 * particle + 1] * particles[3 * particle + 1]
                + constants::beta * particles[3 * particle + 2] * particles[3 * particle + 2];
    }

    return r_squared;
}

double InteractiveElliptical::compute_distance(double *particles, int particle_k, int particle_i) {
    double distance = 
        + (particles[3 * particle_k + 0] - particles[3 * particle_i + 0])*(particles[3 * particle_k + 0] - particles[3 * particle_i + 0])
        + (particles[3 * particle_k + 1] - particles[3 * particle_i + 1])*(particles[3 * particle_k + 1] - particles[3 * particle_i + 1])
        + (particles[3 * particle_k + 2] - particles[3 * particle_i + 2])*(particles[3 * particle_k + 2] - particles[3 * particle_i + 2]);

    return sqrt(distance);
}

double InteractiveElliptical::compute_nabla_alpha(double *particles) {
    return - compute_r_squared(particles);
}

double InteractiveElliptical::compute_u_prime(double r_kj) {
    return a/(r_kj * (r_kj - a));
}

double InteractiveElliptical::compute_u_double_prime(double r_kj) {
    return a * (a - 2*r_kj) / ((r_kj * r_kj) * (r_kj-a) * (r_kj-a));
}

double InteractiveElliptical::compute_local_energy(double *particles) {
    double r_squared_gamma_squared = 0, beta_squared = constants::beta * constants::beta;

    for (int particle = 0; particle < num_particles; particle++) {
        r_squared_gamma_squared += 
                             particles[3 * particle + 0] * particles[3 * particle + 0]
            +                particles[3 * particle + 1] * particles[3 * particle + 1]
            + beta_squared * particles[3 * particle + 2] * particles[3 * particle + 2];
    }

    double potential_energy = 0.5 * constants::omega * r_squared_gamma_squared;

    double kinetic_energy = - 0.5 * constants::omega * compute_laplacian(particles);

    return potential_energy + kinetic_energy;
}

void InteractiveElliptical::compute_quantum_force(double *particles, int particle, double *force) {
    std::cout << "Trying to compute the quantum force for the elliptical case.." << std::endl;
    std::cout << "Particles " << &particles << std::endl;
    std::cout << "Particle " << particle << std::endl;
    std::cout << "Force " << &force << std::endl;
    throw std::logic_error("Not implemented");
}
