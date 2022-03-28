#include "project1/parameters.h"
#include <iomanip>

void write_heading(std::ofstream &file) {
    file << std::setprecision(15);
    file <<
         "num_particles" << "\t" <<
         "dimensions" << "\t" <<
         "alpha" << "\t" <<
         "delta_t" << "\t" <<
         "system_enum" << "\t" <<
         "energy_expectation" << "\t" <<
         "energy_variance" << "\t" <<
         "accepted_ratio" << "\t" <<
         "time_ms" << "\n";
}

Parameters::Parameters(int num_particles, int dimensions, double alpha, double delta_t, SystemEnum system_enum): 
    num_particles(num_particles),
    dimensions(dimensions),
    alpha(alpha),
    delta_t(delta_t),
    system_enum(system_enum),
    energy_filename("")
{}

Parameters::Parameters(int num_particles, int dimensions, double alpha, double delta_t, SystemEnum system_enum, std::string energy_filename):
    num_particles(num_particles),
    dimensions(dimensions),
    alpha(alpha),
    delta_t(delta_t),
    system_enum(system_enum),
    energy_filename(energy_filename)
{}

void Parameters::add_to_file(std::ofstream &file) {
    file <<
         num_particles << "\t" <<
         dimensions << "\t" <<
         alpha << "\t" <<
         delta_t << "\t" <<
         system_enum_to_string(system_enum) << "\t" <<
         energy_expectation << "\t" <<
         energy_variance << "\t" <<
         accepted_ratio << "\t" <<
         time_ms << "\n";
}
