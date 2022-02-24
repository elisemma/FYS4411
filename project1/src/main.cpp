#include <iostream>
#include <vector>
#include "project1/system.h"
#include "project1/particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include <fstream>


using namespace std;

struct Parameters {
    double alpha, delta_t, energy;
    int numberOfDimensions, numberOfParticles;
    bool importanceSampling, useNumerical;

    Parameters(double alpha, int numberOfDimensions, int numberOfParticles, double delta_t, bool importanceSampling, bool useNumerical) :
        alpha(alpha), delta_t(delta_t), numberOfDimensions(numberOfDimensions), numberOfParticles(numberOfParticles), importanceSampling(importanceSampling), useNumerical(useNumerical) {}
};

void gradient_alpha_search(Parameters parameters, int seed, double eta, double delta) {
    //for (auto &parameters : parametersVec) {
    double omega              = 1.0;          // Oscillator frequency.
    double stepLength         = 0.1;
    double equilibration      = 0.1;          // Amount of the total steps used
    int    numberOfSteps      = (int) 1e6;
    double alpha              = parameters.alpha;
    double alpha_change       = 1;

    while (abs(alpha_change) > delta) {
        System* system = new System(
                omega,
                alpha,
                parameters.numberOfDimensions,
                parameters.numberOfParticles,
                equilibration,
                stepLength,
                parameters.useNumerical,
                seed
            );

        system->runMetropolisSteps(numberOfSteps, parameters.delta_t, parameters.importanceSampling);

        alpha_change = eta * system->getAlphaDerivativeChange();
        alpha -= alpha_change;

        cout << "alpha: " << alpha << " alpha_change: " << alpha_change << endl;
    }

}

void do_elliptical(double omega, double alpha, double beta, double a, int numberOfDimensions, int numberOfParticles, double equilibration, double stepLength, int seed, int numberOfSteps, double delta_t, bool importanceSampling) {
    System* system = new System(
            omega,
            alpha,
            beta,
            beta,
            a,
            numberOfDimensions,
            numberOfParticles,
            equilibration,
            stepLength,
            seed
        );
    system->runMetropolisSteps(numberOfSteps, delta_t, importanceSampling);
}

int main() {
    // Seed for the random number generator
    int seed = 42;

    // int    numberOfDimensions = 1;
    // int    numberOfParticles  = 1;
    int    numberOfSteps      = (int) 1e6;
    double omega              = 1.0;          // Oscillator frequency.
    // double alpha              = 0.5;          // Variational parameter.
    double stepLength         = 0.1;          // Metropolis step length.
    double equilibration      = 0.1;          // Amount of the total steps used
    // for equilibration.
    // double delta_t            = 0.001;        // time step

    // vector<int> numberOfDimensionsVec{1, 2, 3};
    vector<int> numberOfDimensionsVec{3};
    // vector<int> numberOfParticlesVec{1, 4, 6};
    vector<int> numberOfParticlesVec{3};
    vector<double> delta_tVec{0.1, 1000};
    // int numberOfParticlesArray[] = {100};

    // double alpha_values[] = {0.01, 0.06, 0.11, 0.16, 0.21, 0.26, 0.31, 0.36, 0.41, 0.46, 0.51, 0.56, 0.61, 0.66, 0.71, 0.76, 0.81, 0.86, 0.91, 0.96};
    //vector<double> alphaVec{0.01, 0.21, 0.41, 0.61, 0.81, 0.96};
    vector<double> alphaVec{0.5};

    vector<bool> importanceSamplingVec{true, false};
    vector<bool> useNumericalVec{true, false};

    vector<Parameters> parametersVec;

    for (double alpha : alphaVec)
        for (int numberOfDimensions : numberOfDimensionsVec)
            for (int numberOfParticles : numberOfParticlesVec)
                for (double delta_t : delta_tVec)
                    for (bool importanceSampling : importanceSamplingVec)
                        for (bool useNumerical : useNumericalVec)
                            parametersVec.push_back(Parameters(alpha, numberOfDimensions, numberOfParticles, delta_t, importanceSampling, useNumerical));

    cout << "NOTE: The parameters were hijacked by an evil space pirate:((" << endl;

    // parametersVec = {Parameters(0.5, 3, 3, 0.1, true, false)};
    // parametersVec = {Parameters(0.1, 1, 1, 0.1, true, false)};
    parametersVec = {Parameters(0.1, 1, 1, 0.1, true, true)};
    Parameters p = parametersVec[0];

    double beta = 2.82843;
    double a = 0.0043*(1-2e-6);
    do_elliptical(omega, p.alpha, beta, a, p.numberOfDimensions, p.numberOfParticles, equilibration, stepLength, seed, numberOfSteps, p.delta_t, p.importanceSampling);

    // double eta   = 1e-1;
    // double delta = 1e-8;
    // gradient_alpha_search(parametersVec[0], seed, eta, delta);
    // exit(68);
    //
    // #pragma omp parallel for schedule(dynamic)
    // for (auto &parameters : parametersVec) {
    //     System* system = new System(
    //         omega,
    //         parameters.alpha,
    //         parameters.numberOfDimensions,
    //         parameters.numberOfParticles,
    //         equilibration,
    //         stepLength,
    //         parameters.useNumerical,
    //         seed
    //     );
    //
    //     system->runMetropolisSteps(numberOfSteps, parameters.delta_t, parameters.importanceSampling);
    //     double energy = system->getEnergy();
    //
    //     parameters.energy = energy;
    // }
    //
    // ofstream energyFile("output/data/energy_values.tsv");
    // energyFile << "alpha\tnumberOfDimensions\tnumberOfParticles\tdt\timportanceSampling\tenergy" << endl;
    // for (auto parameters : parametersVec)
    //     energyFile << parameters.alpha << "\t" << parameters.numberOfDimensions << "\t" << parameters.numberOfParticles << "\t" << parameters.delta_t << "\t" << parameters.importanceSampling << "\t" << parameters.energy << endl;
    // energyFile.close();

    return 0;
}
