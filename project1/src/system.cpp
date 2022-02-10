#include "project1/system.h"
#include <cassert>
#include "project1/sampler.h"
#include "project1/particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <math.h>


System::System() {
    m_random = new Random();
}

System::System(int seed) {
    m_random = new Random(seed);
}

// TODO: Should this value be saved/cached?
double System::calculate_r_squared(std::vector<Particle*> particles) {
    double r_squared = 0;
    for (Particle* particle : particles) {
        for (double pos_i : particle->getPosition()) {
            r_squared += pow(pos_i, 2);
        }
    }
    return r_squared;
}

bool System::metropolisStep(bool importance) {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

    Particle *particle = m_particles[m_random->nextInt(m_numberOfParticles - 1)];

    double movement[m_numberOfDimensions], wave_before = m_waveFunction->evaluate(m_particles);
    std::vector<double> quantum_before;

    // TODO: Check that importance step is supposed to only do one particle at a time
    // TODO: We want some monte carlo action <3??
    // TODO: Check that the quantum_before is actually updated here!
    if (importance) quantum_before = m_waveFunction->computeQuantumForceAnalytical(particle); // TODO: GET THE QUANTUM FORCE!

    for (int i = 0; i < m_numberOfDimensions; i++) {
        // TODO: should this be -.5 or *2-1?
        movement[i] = m_stepLength*(m_random->nextDouble() - 0.5);
        particle->adjustPosition(movement[i], i);
    }

    double wave_after = m_waveFunction->evaluate(m_particles);

    double ratio = wave_after*wave_after/wave_before*wave_before;

    // TODO: Check the the implementation is efficient
    if (importance) {
        std::vector<double> quantum_after = m_waveFunction->computeQuantumForceAnalytical(particle);
        // TODO: get the quantum force after the move
        // do some green shroom ratio or smtn

        double delta_t = 0.1;
        double D = 0.5;
        // double F_after = m_waveFunction->computeQuantumForceAnalytical(particle);
        //TODO: N = numberOfParticles???
        double N = 1;
        //TODO: Check minus signs!!!
        double G = 1/pow(4*M_PI*D*delta_t, 3*N/2)*exp(pow(movement + D*delta_t*quantum_after,2)/(4*D*delta_t));

    }

    // should we check if ratio is more than 1?
    if (m_random->nextDouble() < ratio) return true;

    for (int i = 0; i < m_numberOfDimensions; i++) particle->adjustPosition(-movement[i], i);

    return false;


}

double System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        // TODO: Check if metropolis step is used..
        bool acceptedStep = metropolisStep(true);

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

    return m_sampler->getEnergy();
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
