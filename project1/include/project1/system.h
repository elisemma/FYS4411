#pragma once

#include <vector>
#include <random>
#include <fstream>

#include "project1/random_engine.h"
#include "project1/wave_functions/wave_function.h"
#include "project1/parameters.h"

class System {
    public:
        System(Parameters parameters);
        ~System();
        void run_simulation();

        // Variables for output
        double accepted_ratio, energy_expectation, energy_variance, alpha_derivative;

        void write_system_data_to_file();
    private:
        // Functions
        bool metropolis_step();
        bool importance_metropolis_step();
        void initialize_particles(int num_particles, int dimensions);

        // Other classes
        WaveFunction* wave_function;
        RandomEngine random_engine;

        // General variables
        double *particles;
        double *one_body_densities;
        int num_particles, dimensions;
        double delta_t;
        double *energies;
        bool importance_sampling;
        std::string system_data_filename;
};
