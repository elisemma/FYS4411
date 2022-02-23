#pragma once

#include "project1/particle.h"

class WaveFunction {
public:
    WaveFunction(class System* system);
    // int     getNumberOfParameters() { return m_numberOfParameters; }
    // std::vector<double> getParameters() { return m_parameters; }
    double getAlpha() { return alpha; }
    double computeDoubleDerivativeNumerical(std::vector<Particle*> particles);
    void move_particles(std::vector<Particle*> particles, double step_length, int dim);
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    // virtual double computeDoubleDerivative(std::vector<class Particle*> particles, bool analytical) = 0;
    virtual double computeDoubleDerivativeAnalytical(std::vector<class Particle*> particles) = 0;
    // virtual double computeDoubleDerivativeNumerical(std::vector<class Particle*> particles) = 0;
    virtual std::vector<double> computeQuantumForceAnalytical(Particle* particles) = 0;
    virtual double computeDerivative(std::vector<class Particle*> particles) = 0;
protected:
    // int     m_numberOfParameters = 0;
    // std::vector<double> m_parameters = std::vector<double>();
    double alpha;
    class System* m_system = nullptr;
};
