#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    double getAlphaDerivativeChange();
    double getEnergy() { return m_energy; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_expectedTrailEnergyDerivative = 0;
    double  m_expectedTrailDerivative = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergySquared = 0;
    double  m_cululativeTrailDerivative = 0;
    double  m_cumulativeTrailEnergyDerivative = 0;
    class System* m_system = nullptr;
};
