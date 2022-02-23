#include "WaveFunctions/interactive.h"
#include "WaveFunctions/wavefunction.h"
#include "project1/system.h"
#include "project1/particle.h"
#include <cmath>
#include <cassert>

using namespace std;

Interactive::Interactive(System* system, double alpha_, double gamma_squared_, double a_) {
    assert(alpha_ >= 0);
    assert(gamma_squared_ >= 0);

    this->alpha = alpha_;
    this->gamma_squared = gamma_squared_;
    this->a = a_;
}

double Interactive::evaluate(vector<Particle*> particles) {
    double r_squared = 0;

    for (auto particle : particles) {
        vector<double> particle_pos = particle->getPosition();
        for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
            if (dim == 2) {
                r_squared += particle_pos[2] * particle_pos[2] * gamma_squared;
            } else {
                r_squared += particle_pos[dim] * particle_pos[dim];
            }
        }
    }
}
