#include "project1/system.h"
#include <cassert>
#include "project1/sampler.h"
#include "project1/particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <math.h>


#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/elliptical.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/interactive.h"
#include "InitialStates/randomuniform.h"

using namespace std;

System::System(double omega, double alpha, int numberOfDimensions, int numberOfParticles, double equilibration, double stepLength, bool useNumerical, int seed) {
    m_random = new Random(seed);
    setHamiltonian                     (new HarmonicOscillator(this, omega, useNumerical));
    setWaveFunction                    (new SimpleGaussian(this, alpha));
    setInitialState                    (new RandomUniform(this, numberOfDimensions, numberOfParticles));
    setEquilibrationFraction           (equilibration);
    setStepLength                      (stepLength);
}

System::System(double omega, double alpha, double beta, double gamma, double a, int numberOfDimensions, int numberOfParticles, double equilibration, double stepLength, int seed) {
    assert(beta == gamma);

    m_random = new Random(seed);

    setHamiltonian                     (new Elliptical(this, omega, beta, a));
    setWaveFunction                    (new Interactive(this, alpha, beta, a));
    // TODO: Check that the particles are not overlapping for the interactive case, using 
        // for (int j = 0; j < number_of_particles - 1; j++) {
        //     for (int k = j + 1; k < number_of_particles; k++) {
        //         double distance = m_system->getDistance(j, k);
        //
        //         if (distance <= m_a) {
        //             RAISE ERROR, "The particles are overlapping!";
        //             REMAKE THE SYSTEM IN ANOTHER CONFIGURATION, TO AVOID THE OVERLAPPING PARTICLES
        //         }
        //     }
        // }
    setInitialState                    (new RandomUniform(this, numberOfDimensions, numberOfParticles));
    setEquilibrationFraction           (equilibration);
    setStepLength                      (stepLength);
}

// TODO: Only update this value
double System::calculate_r_squared(std::vector<Particle*> particles) {
    double r_squared = 0;
    for (Particle* particle : particles) {
        for (double pos_i : particle->getPosition()) {
            r_squared += pos_i * pos_i;
        }
    }
    return r_squared;
}

bool System::metropolisStep(bool importance, double delta_t) {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

    Particle *particle = m_particles[m_random->nextInt(m_numberOfParticles - 1)];

    double movement[m_numberOfDimensions], wave_before = m_waveFunction->evaluate(m_particles);
    std::vector<double> quantum_before;

    if (importance) quantum_before = m_waveFunction->computeQuantumForceAnalytical(particle);

    for (int i = 0; i < m_numberOfDimensions; i++) {
        movement[i] = m_stepLength*(m_random->nextDouble() - 0.5);
        particle->adjustPosition(movement[i], i);
    }

    double wave_after = m_waveFunction->evaluate(m_particles);

    double ratio = wave_after*wave_after/wave_before*wave_before;

    // TODO: Check the the implementation is efficient
    if (importance) {
        std::vector<double> quantum_after = m_waveFunction->computeQuantumForceAnalytical(particle);
        // do some green shroom ratio or smtn

        //double delta_t = 0.01;
        double D = 0.5;
        // double F_after = m_waveFunction->computeQuantumForceAnalytical(particle);
        double N = getNumberOfParticles();
        double d = getNumberOfDimensions();
        double G_yx = 0, G_xy = 0;

        for (int j = 0; j < m_numberOfDimensions; j++){
            G_yx += pow(1/pow(4*M_PI*D*delta_t, d*N/2)*exp(-pow( movement[j] - D*delta_t*quantum_before[j],2)/(4*D*delta_t)),2);
            G_xy += pow(1/pow(4*M_PI*D*delta_t, d*N/2)*exp(-pow(-movement[j] - D*delta_t*quantum_after[j] ,2)/(4*D*delta_t)),2);
        }

        //double q = sqrt(G_yx/G_xy)*ratio;
        ratio *= sqrt(G_yx/G_xy);

    }

    if (m_random->nextDouble() < ratio) return true;

    for (int i = 0; i < m_numberOfDimensions; i++) particle->adjustPosition(-movement[i], i);

    return false;


}

void System::runMetropolisSteps(int numberOfMetropolisSteps, double delta_t, bool importanceSampling) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep(importanceSampling, delta_t);

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */
        m_sampler->sample(acceptedStep);
    }
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();

    m_alphaDerivativeChange = m_sampler->getAlphaDerivativeChange();
    m_energy = m_sampler->getEnergy();
    m_energyVariance = m_sampler->getEnergyVariance();
    m_acceptedMetropolisStepRatio = m_sampler->getAcceptedMetropolisStepRatio();

    // return m_sampler->getEnergy();
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}

// TODO: Calculate all the distances at once for improved speed
// double System::getDistance(std::vector<Particle*> particles, int i, int j) {
double System::getDistance(int i, int j) {
    vector<double> pos_i = m_particles[i]->getPosition();
    vector<double> pos_j = m_particles[j]->getPosition();

    double distance = 0;
    for (int dim = 0; dim < getNumberOfDimensions(); dim++) {
        distance += abs(pos_i[dim] - pos_j[dim]);
    }

    return distance;
}
