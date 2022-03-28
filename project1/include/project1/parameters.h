#pragma once

#include <fstream>
#include "project1/system_enum.h"

void write_heading(std::ofstream &file);

struct Parameters {
    Parameters(int num_particles, int dimensions, double alpha, double delta_t, SystemEnum system_enum);
    Parameters(int num_particles, int dimensions, double alpha, double delta_t, SystemEnum system_enum, std::string energy_filename);

    int num_particles, dimensions;
    double alpha, delta_t;
    bool importance_sampling;
    SystemEnum system_enum;
    std::string energy_filename;

    // Output parameters
    double energy_expectation, energy_variance, accepted_ratio;
    int time_ms;

    void add_to_file(std::ofstream &file);
};
