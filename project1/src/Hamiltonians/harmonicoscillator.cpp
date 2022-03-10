#include <iostream>
#include "Hamiltonians/harmonicoscillator.h"
#include <math.h>
#include <cassert>
#include "project1/system.h"
#include "project1/particle.h"
#include "WaveFunctions/wavefunction.h"

HarmonicOscillator::HarmonicOscillator(System* system, double omega, bool use_numerical) :Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
    m_use_numerical = use_numerical;
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
    /* You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    double double_derivative; //, r_squared = m_system->calculate_r_squared(particles);

    double potentialEnergy = 0.5 * (m_omega * m_omega) * m_system->getRSquared();

    if (m_use_numerical) {
        double_derivative = m_system->getWaveFunction()->computeDoubleDerivativeNumerical(particles);
    } else {
        double_derivative = m_system->getWaveFunction()->computeDoubleDerivativeAnalytical(particles);
    }

    //double_derivative = m_system->getWaveFunction()->computeDoubleDerivative(particles, false);

    double kineticEnergy = -0.5 * double_derivative;

    // cout << "r2: " << r_squared << endl;
    // cout << "potentialEnergy: " << potentialEnergy << endl;
    // cout << "kineticEnergy: "   << kineticEnergy   << endl;

    return kineticEnergy + potentialEnergy;
}
