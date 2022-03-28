#pragma once

#include "project1/wave_functions/wave_function.h"

class Simple : public WaveFunction {
    public:
        Simple(int num_particles, int dimensions, double alpha);
        double evaluate(double *particles);
        double compute_nabla_alpha(double *particles);
        double compute_local_energy(double *particles);
        void compute_quantum_force(double *particles, int particle, double *force);
        virtual double compute_laplacian(double *particles) = 0;
    protected:
        double compute_r_squared(double *particles);
        int num_particles;
        int dimensions;
        double alpha;
};
