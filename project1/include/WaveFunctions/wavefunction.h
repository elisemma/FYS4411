#pragma once

#include <vector>

class WaveFunction {
public:
    WaveFunction(class System* system);
    // int     getNumberOfParameters() { return m_numberOfParameters; }
    // std::vector<double> getParameters() { return m_parameters; }
    double getAlpha() { return alpha; }
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles, bool analytical) = 0;
protected:
    // int     m_numberOfParameters = 0;
    // std::vector<double> m_parameters = std::vector<double>();
    double alpha;
    class System* m_system = nullptr;
};
