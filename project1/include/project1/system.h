#pragma once
#include <vector>
#include <Math/random.h>

class System {
public:
    System(double omega, double alpha, int numberOfDimensions, int numberOfParticles, double equilibration, double stepLength);
    System(double omega, double alpha, int numberOfDimensions, int numberOfParticles, double equilibration, double stepLength, int seed);

    // System();
    // System(int seed);
    bool metropolisStep             (bool importance, double delta_t);
    void runMetropolisSteps         (int numberOfMetropolisSteps, double delta_t, bool importanceSampling);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    class Random*                   getRandomEngine()   { return m_random; }
    double getEnergy                () { return m_energy; };
    double getEnergyVariance        () { return m_energyVariance; };
    double getStepLength()              { return m_stepLength; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    double calculate_r_squared(std::vector<Particle*> particles);
    double getAlphaDerivativeChange()   {return m_alphaDerivativeChange; }
    double getDistance(int i, int j);


private:
    void initialize_system(double omega, double alpha, int numberOfDimensions, int numberOfParticles, double equilibration, double stepLength);
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
    class Random*                   m_random = nullptr;
    double                          m_alphaDerivativeChange;
    double                          m_energy;
    double                          m_energyVariance;

};
