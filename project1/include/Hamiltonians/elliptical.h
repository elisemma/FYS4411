#pragma once

#include "hamiltonian.h"
#include <vector>

class Elliptical : public Hamiltonian {
public:
    Elliptical(System* system, double omega, double gamma, double a);
    double computeLocalEnergy(std::vector<Particle*> particles);
private:
    double getDistance(std::vector<class Particle*> particles, int i, int j);
    double m_omega = 0;
    double m_gamma_squared = 0;
    double m_a = 0;
};
