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

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

    Particle *particle = m_particles[m_random->nextInt(m_numberOfParticles - 1)];

    double movement[m_numberOfDimensions], wave_before = m_waveFunction->evaluate(m_particles);

	for (int i = 0; i < m_numberOfDimensions; i++) {
        // TODO: should this be -.5 or *2-1?
        movement[i] = m_stepLength*(m_random->nextDouble() - 0.5);
        particle->adjustPosition(movement[i], i);
	}

	double wave_after = m_waveFunction->evaluate(m_particles);
    // TODO: Do we need them pows?
	double ratio = pow(wave_after, 2)/pow(wave_before, 2);

    // should we check if ratio is more than 1?
    if (m_random->nextDouble() < ratio) return true;

    for (int i = 0; i < m_numberOfDimensions; i++) particle->adjustPosition(-movement[i], i);

    return false;
    // FOKKER_PLANCK:
    // Vi tror at y er den nye og x er den naavaerende posisjonen
    //for (int i = 0; i < m_numberOfDimensions; i++) {
          // TODO: should this be -.5 or *2-1?
          //movement[i] = m_stepLength*(m_random->nextDouble() - 0.5);
          //particle->adjustPosition(movement[i], i);
          //x =
          //y =
  	}
    //double F_x = -4*alpha*r
    //double delta_t = 0.1;
    //double G_yx = 1/pow(2*M_PI*delta_t, 3*numberOfParticles/2)*exp(-pow(y-x - 0.5*delta_t*F_x,2)/(2*delta_t));
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();

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
