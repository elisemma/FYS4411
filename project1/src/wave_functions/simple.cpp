#include <cmath>
#include "project1/wave_functions/simple_analytical.h"
#include "project1/constants.h"

Simple::Simple(int num_particles, int dimensions, double alpha): num_particles(num_particles), dimensions(dimensions), alpha(alpha) {}

double Simple::compute_r_squared(double *particles) {
    double r_squared = 0;

    for (int i = 0; i < num_particles * dimensions; i++) {
        r_squared += particles[i] * particles[i];
    }

    return r_squared;
}

double Simple::evaluate(double *particles) {
    return exp(-alpha * compute_r_squared(particles));
}

double Simple::compute_nabla_alpha(double *particles) {
    return -compute_r_squared(particles);
}

double Simple::compute_laplacian(double *particles) {
    return -2 * num_particles * dimensions * alpha + 4 * alpha * alpha * compute_r_squared(particles);
}

double Simple::compute_local_energy(double *particles) {
    double r_squared = compute_r_squared(particles);

    double potential_energy = 0.5 * constants::omega * constants::omega * r_squared;

    double kinetic_energy = - 0.5 * compute_laplacian(particles);

    return potential_energy + kinetic_energy;
}

void Simple::compute_quantum_force(double *particles, int particle, double *force) {
    for (int dim = 0; dim < dimensions; dim++) {
        force[dim] = - 4 * alpha * particles[particle * dimensions + dim];
    }
}
