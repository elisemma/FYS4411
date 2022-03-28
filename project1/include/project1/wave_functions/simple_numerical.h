#pragma once

#include "project1/wave_functions/simple.h"

class SimpleNumerical : public Simple {
    public:
        SimpleNumerical(int num_particles, int dimensions, double alpha);
    private:
        double compute_laplacian(double *particles);
};
