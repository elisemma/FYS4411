#pragma once

#include <random>

class RandomEngine {

public:
    RandomEngine(int num_particles);

    int random_particle();
    double next_double();
    double next_double(double a, double b);
    double next_normal(double mean, double standardDeviation);
private:
    std::mt19937_64 random_engine;
    int num_particles;
};

