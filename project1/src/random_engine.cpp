#include "project1/random_engine.h"
#include "project1/constants.h"

RandomEngine::RandomEngine(int num_particles) : num_particles(num_particles) {
    random_engine = std::mt19937_64(constants::seed);
}

int RandomEngine::random_particle() {
    return std::uniform_int_distribution<int>(0, num_particles - 1)(random_engine);
}

double RandomEngine::next_double() {
    // Produces uniformly distributed random floating-point values in the
    // half-open interval [0, 1).

    return std::uniform_real_distribution<double>(0.0, 1.0)(random_engine);
}

double RandomEngine::next_double(double a, double b) {
    // Produces uniformly distributed random floating-point values in the
    // half-open interval [a, b).

    return std::uniform_real_distribution<double>(a, b)(random_engine);
}

double RandomEngine::next_normal(double mean, double standardDeviation) {
    // Produces normal distributed random floating-point values with mean
    // ``mean`` and standard deviation ``standardDeviation``.

    return std::normal_distribution<double>(mean, standardDeviation)(random_engine);
}
