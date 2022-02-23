#include "Hamiltonians/elliptical.h"
#include <cassert>
#include "project1/system.h"
#include "project1/particle.h"
#include "WaveFunctions/wavefunction.h"

using namespace std;

Elliptical::Elliptical(System* system, double omega, double gamma, double a) : Hamiltonian(system) {
    assert(omega > 0);

    m_omega = omega;
    m_gamma_squared = gamma * gamma;
    m_a = a;
}

double Elliptical::computeLocalEnergy(std::vector<Particle*> particles) {
    double r_squared = 0;

    for (auto particle : particles) {
        vector<double> particle_pos = particle->getPosition();
        for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
            if (dim == 2) {
                r_squared += particle_pos[2] * particle_pos[2] * m_gamma_squared;
            } else {
                r_squared += particle_pos[dim] * particle_pos[dim];
            }
        }
    }

    // TODO: Switch between numerical and analytical in a nicer way :}
    double nabla_squared = m_system->getWaveFunction()->computeDoubleDerivativeNumerical(particles);

    // TODO: Give is a better name
    double external_idn = 0.5 * (-nabla_squared + r_squared);

    int number_of_particles = particles.size();
    double v_int_sum = 0;
    for (int j = 0; j < number_of_particles - 1; j++) {
        for (int k = j + 1; k < number_of_particles; k++) {
            double distance = m_system->getDistance(j, k);

            // TODO: Increase this number plz:)
            if (distance <= m_a) {
                v_int_sum = 1e80;
                break;
            }
        }
    }
    
    return external_idn + v_int_sum;
}

// // TODO: Move getDistance to a nicer place, and use the same in interactive.cpp
// double Elliptical::getDistance(std::vector<class Particle*> particles, int i, int j) {
//     vector<double> pos_i = particles[i]->getPosition();
//     vector<double> pos_j = particles[j]->getPosition();
//
//     double distance = 0;
//     for (int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
//         distance += abs(pos_i[dim] - pos_j[dim]);
//     }
//
//     return distance;
// }
