#include <cmath>
#include "project1/wave_functions/simple_numerical.h"
#include "project1/constants.h"

SimpleNumerical::SimpleNumerical(int num_particles, int dimensions, double alpha): Simple(num_particles, dimensions, alpha) {}

double SimpleNumerical::compute_laplacian(double *particles) {
    double double_derivative = 0;

    for (int particle = 0; particle < num_particles; particle++) {
        for (int dimension = 0; dimension < dimensions; dimension++) {
            particles[particle * dimensions + dimension] += constants::numerical_h;
            double_derivative += evaluate(particles);

            particles[particle * dimensions + dimension] -= 2 * constants::numerical_h;
            double_derivative += evaluate(particles);

            particles[particle * dimensions + dimension] += constants::numerical_h;
        }
    }
    double evaluated = evaluate(particles);
    double_derivative -= 2 * dimensions * num_particles * evaluated;
    double_derivative /= constants::numerical_h * constants::numerical_h * evaluated;

    return double_derivative;
}
