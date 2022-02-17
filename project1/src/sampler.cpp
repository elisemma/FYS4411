#include <iostream>
#include <cmath>
#include <vector>
#include "project1/sampler.h"
#include "project1/system.h"
#include "project1/particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

// TODO: What to do with accepted step?
void Sampler::sample(bool acceptedStep) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
    }

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */
    double localEnergy = m_system->getHamiltonian()->
                         computeLocalEnergy(m_system->getParticles());

    m_cumulativeEnergy  += localEnergy;

    double trail_derivative = m_system->getWaveFunction()->computeDerivative(m_system->getParticles());

    m_cululativeTrailDerivative += trail_derivative;
    m_cumulativeTrailEnergyDerivative += trail_derivative * localEnergy;

    m_stepNumber++;
}

double Sampler::getAlphaDerivativeChange() {
    double steps = m_system->getNumberOfMetropolisSteps();

    m_energy = m_cumulativeEnergy / steps;

    // Term 1
    double m_expectedTrailEnergyDerivative = m_cumulativeTrailEnergyDerivative / steps;

    // Term 2 (part 1)
    double m_expectedTrailDerivative = m_cululativeTrailDerivative / steps;

    // Final
    // TODO: We need eta
    return 0.01 * 2 * (m_expectedTrailEnergyDerivative - m_expectedTrailDerivative * m_energy);
}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    // int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    // std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << "Alpha: " << m_system->getWaveFunction()->getAlpha() << endl;
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << endl;
}
