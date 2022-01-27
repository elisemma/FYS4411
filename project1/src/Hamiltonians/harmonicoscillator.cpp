#include "Hamiltonians/harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "project1/system.h"
#include "project1/particle.h"
#include "WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}

// double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */
//
//     double potentialEnergy = 0;
//     double kineticEnergy   = 0;
//     return kineticEnergy + potentialEnergy;
// }



double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    double r_squared = 0, double_derivative = 0;

    // USE THE SAME AS IN simplegaussian.cpp
    for (Particle* particle : particles) {
        for (double pos_i : particle->getPosition()) {
            r_squared += pow(pos_i, 2);
        }
    }
    // TODO: m_omega eller m_omega * m_omega, we do be lovin us some pow action!<3
    double potentialEnergy = 0.5 * pow(m_omega, 2) * r_squared;


    // TODO: Check if numeric or analytical
    // TODO: the comments say something about analytical or numerical, here we only use analytical
    double_derivative = m_system->getWaveFunction()->computeDoubleDerivative(particles);

    double kineticEnergy = -0.5 * double_derivative;

    return kineticEnergy + potentialEnergy;
}
