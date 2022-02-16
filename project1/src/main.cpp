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
    double delta_t            = 0.001;        // time step

    // int numberOfParticlesArray[] = {1, 10, 100, 500};
    vector<int> numberOfDimensionsVec{1, 2, 3};
    vector<int> numberOfParticlesVec{1, 4, 6};
    vector<double> delta_tVec{0.1, 1000}; // TODO: Hvorfor endres ikke energy når vi endrer dt??????
    // int numberOfParticlesArray[] = {100};

    // TODO: make this in a nicer way?
    // double alpha_values[] = {0.01, 0.06, 0.11, 0.16, 0.21, 0.26, 0.31, 0.36, 0.41, 0.46, 0.51, 0.56, 0.61, 0.66, 0.71, 0.76, 0.81, 0.86, 0.91, 0.96};
    //vector<double> alphaVec{0.01, 0.21, 0.41, 0.61, 0.81, 0.96};
    vector<double> alphaVec{0.5};


    vector<vector<size_t>> config_values;

    for (size_t a = 0; a < alphaVec.size(); a++) {
        for (size_t d = 0; d < numberOfDimensionsVec.size(); d++) {
            for (size_t n = 0; n < numberOfParticlesVec.size(); n++) {
              for (size_t dt = 0; dt < delta_tVec.size(); dt++){

                  vector<size_t> config = {a, d, n, dt};
                  config_values.push_back(config);
                }
            }
        }
    }

    double energyValues[alphaVec.size()][numberOfDimensionsVec.size()][numberOfParticlesVec.size()][delta_tVec.size()];
    //
    // // TODO: Loop over dimensions from 1..3
    #pragma omp parallel for schedule(dynamic)
    for (auto config : config_values) {
        int a = config[0];
        int d = config[1];
        int n = config[2];
        int dt = config[3];

        System* system = new System(seed);
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new SimpleGaussian(system, alphaVec[a]));
        // system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setInitialState             (new RandomUniform(system, numberOfDimensionsVec[d], numberOfParticlesVec[n]));
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);
        double energy = system->runMetropolisSteps          (numberOfSteps, delta_t);
        // cout << alphaVec[a] << " " << numberOfParticlesVec[n] << " " << energy << endl;
        energyValues[a][d][n][dt] = energy;
    }

    ofstream energyFile("output/data/energy_values.tsv");
    energyFile << "alpha\tnumberOfDimensions\tnumberOfParticles\tdt\tenergy" << endl;
    for (size_t a = 0; a < alphaVec.size(); a++) {
        for (size_t d = 0; d < numberOfDimensionsVec.size(); d++) {
            for (size_t n = 0; n < numberOfParticlesVec.size(); n++) {
                for (size_t dt = 0; dt < delta_tVec.size(); dt++) {
                  energyFile << alphaVec[a] << "\t" << numberOfDimensionsVec[d] << "\t" << numberOfParticlesVec[n] << "\t" << delta_tVec[dt] << "\t" << energyValues[a][d][n][dt] << endl;
                }
            }
        }
    }
    energyFile.close();

    // for (double alpha : alpha_values) {
    //     for (int numberOfParticles : numberOfParticlesArray) {
    //         System* system = new System(seed);
    //         system->setHamiltonian              (new HarmonicOscillator(system, omega));
    //         system->setWaveFunction             (new SimpleGaussian(system, alpha));
    //         // system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    //         system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    //         system->setEquilibrationFraction    (equilibration);
    //         system->setStepLength               (stepLength);
    //         double energy = system->runMetropolisSteps          (numberOfSteps);
    //         cout << alpha << " " << numberOfParticles << " " << energy << endl;
    //     }
    // }

    return 0;
}
