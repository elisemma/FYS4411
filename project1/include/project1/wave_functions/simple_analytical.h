#pragma once

#include "project1/wave_functions/wave_function.h"
#include "project1/wave_functions/simple.h"

class SimpleAnalytical : public Simple {
    public:
        SimpleAnalytical(int num_particles, int dimensions, double alpha);
    private:
        double compute_laplacian(double *particles);
};
