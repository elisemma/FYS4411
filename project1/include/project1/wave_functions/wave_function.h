#pragma once

class WaveFunction {
    public:
        virtual ~WaveFunction() = default;
        virtual double evaluate(double *particles) = 0;
        virtual double compute_nabla_alpha(double *particles) = 0;
        virtual double compute_local_energy(double *particles) = 0;
        virtual void compute_quantum_force(double *particles, int particle, double *force) = 0;
};
