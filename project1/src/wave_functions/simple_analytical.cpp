#include <cmath>
#include "project1/wave_functions/simple_analytical.h"
#include "project1/constants.h"

SimpleAnalytical::SimpleAnalytical(int num_particles, int dimensions, double alpha): Simple(num_particles, dimensions, alpha) {}

double SimpleAnalytical::compute_laplacian(double *particles) {
    return -2 * num_particles * dimensions * alpha + 4 * alpha * alpha * compute_r_squared(particles);
}
