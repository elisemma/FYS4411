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

void gradient_alpha_search(Parameters parameters, int seed){
  //for (auto &parameters : parametersVec) {
  double omega              = 1.0;          // Oscillator frequency.
  double stepLength         = 0.1;
  double equilibration      = 0.1;          // Amount of the total steps used
  int    numberOfSteps      = (int) 1e6;


  System* system = new System(
          omega,
          parameters.alpha,
          parameters.numberOfDimensions,
          parameters.numberOfParticles,
          equilibration,
          stepLength,
          seed
      );

  system->runMetropolisSteps(numberOfSteps, parameters.delta_t, parameters.importanceSampling);
  // double energy = system->getEnergy();
  // double energy = system->runMetropolisSteps(numberOfSteps, parameters.delta_t, parameters.importanceSampling);

//TODO: Verify if h=0.001 is ok as step
  double alpha_new = parameters.alpha - 0.001*system->getAlphaDerivativeChange();
  cout << "alpha_new" << alpha_new << endl;
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
    // vector<double> delta_tVec{0.1, 1000}; // TODO: Hvorfor endres ikke energy når vi endrer dt??????
    vector<double> delta_tVec{0.1, 1000}; // TODO: Hvorfor endres ikke energy når vi endrer dt??????
    // int numberOfParticlesArray[] = {100};

    // TODO: make this in a nicer way?
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

    parametersVec = {Parameters(0.3, 3, 3, 0.1, true, false)};

    // gradient_alpha_search(parametersVec[0], seed);
    // exit(69);

    #pragma omp parallel for schedule(dynamic)
    for (auto &parameters : parametersVec) {
        // TODO: Actually do something if you have use numberical
        System* system = new System(
            omega,
            parameters.alpha,
            parameters.numberOfDimensions,
            parameters.numberOfParticles,
            equilibration,
            stepLength,
            seed
        );

        system->runMetropolisSteps(numberOfSteps, parameters.delta_t, parameters.importanceSampling);
        double energy = system->getEnergy();
        // double energy = system->runMetropolisSteps(numberOfSteps, parameters.delta_t, parameters.importanceSampling);

        parameters.energy = energy;
    }

    ofstream energyFile("output/data/energy_values.tsv");
    energyFile << "alpha\tnumberOfDimensions\tnumberOfParticles\tdt\timportanceSampling\tenergy" << endl;
    for (auto parameters : parametersVec)
        energyFile << parameters.alpha << "\t" << parameters.numberOfDimensions << "\t" << parameters.numberOfParticles << "\t" << parameters.delta_t << "\t" << parameters.importanceSampling << "\t" << parameters.energy << endl;
    energyFile.close();

    return 0;
}
