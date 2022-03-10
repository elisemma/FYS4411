#pragma once

#include "hamiltonian.h"
#include <vector>

class Elliptical : public Hamiltonian {
public:
    // Elliptical(System* system, double omega, double gamma, double a);
    Elliptical(System* system, double omega, double gamma, double a, bool useNumerical);
    double computeLocalEnergy(std::vector<Particle*> particles);
private:
    // double getDistance(std::vector<class Particle*> particles, int i, int j);
    bool useNumerical;
    double m_omega = 0;
    double m_gamma_squared = 0;
    double m_a = 0;
};
