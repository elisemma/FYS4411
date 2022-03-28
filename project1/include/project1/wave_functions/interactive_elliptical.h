#pragma once

#include "project1/wave_functions/wave_function.h"

class InteractiveElliptical : public WaveFunction {
    public:
        InteractiveElliptical(int num_particles, int dimensions, double alpha, double a);
        double evaluate(double *particles);
        double compute_laplacian(double *particles);
        double compute_nabla_alpha(double *particles);
        double compute_local_energy(double *particles);
        void compute_quantum_force(double *particles, int particle, double *force);
    private:
        double compute_r_squared(double *particles);
        double particle_squared(double *particles, int particle);
        double compute_u_prime(double r_kj);
        double compute_u_double_prime(double r_kj);

        // Helper methods
        double compute_distance(double *particles, int k, int j);

        // Parameters
        int num_particles;
        double alpha;
        double a;
};
