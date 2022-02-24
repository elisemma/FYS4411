#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy() { return m_energy; }
    double getEnergyVariance() { return m_energyVariance; }
    double getAlphaDerivativeChange() { return m_alphaDerivativeChange; }
    double getAcceptedMetropolisStepRatio() { return m_acceptedMetropolisStepRatio; }

private:
    int     m_numberOfMetropolisSteps         = 0;
    int     m_stepNumber                      = 0;
    double  m_energy                          = 0;
    double  m_energyVariance                  = 0;
    double  m_cumulativeEnergy                = 0;
    double  m_cumulativeEnergySquared         = 0;
    double  m_alphaDerivativeChange           = 0;
    double  m_cumulativeTrailDerivative       = 0;
    double  m_cumulativeTrailEnergyDerivative = 0;
    double  m_numberOfAcceptedMetropolisSteps = 0;
    double  m_acceptedMetropolisStepRatio     = 0;
    class System* m_system                    = nullptr;
};
