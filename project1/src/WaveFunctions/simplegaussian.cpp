#include "WaveFunctions/simplegaussian.h"
#include <cmath>
#include <cassert>
#include "WaveFunctions/wavefunction.h"
#include "project1/system.h"
#include "project1/particle.h"

using namespace std;

SimpleGaussian::SimpleGaussian(System* system, double alpha_) : WaveFunction(system) {
    assert(alpha_ >= 0);
    this->alpha = alpha_;
}

double SimpleGaussian::evaluate(vector<Particle*> particles) {
    double r_squared = m_system->calculate_r_squared(particles);
    return exp(-alpha * r_squared);
}

double SimpleGaussian::computeDoubleDerivative(vector<class Particle*> particles, bool analytical) {
    if (analytical) {
      return computeDoubleDerivativeAnalytical(particles);
    } else {
      return computeDoubleDerivativeNumerical(particles);
    }
}

double SimpleGaussian::computeDoubleDerivativeAnalytical(vector<class Particle*> particles) {
    double r_squared = m_system->calculate_r_squared(particles);
    return 2*alpha*exp(-alpha*r_squared)*(2*alpha*r_squared - 1);
}

//double SimpleGaussian::computeDoubleDerivativeNumerical(vector<class Particle*> particles, double[] r_vec_sq){
  //double double_derivative = (exp(-2*alpha*r_vec_sq(2)) - 2*exp(-2*alpha*r_vec_sq(1)) + exp(-2*alpha*r_vec_sq(0)))/pow(m_stepLength,2);
//}
double SimpleGaussian::computeDoubleDerivativeNumerical(vector<class Particle*> particles){
  //TO DO: lage h og r_squared_list
  //(exp(-2*alpha*r_squared_list(2)) -2*exp(-2*alpha*r_squared_list(1)) + exp(-2*alpha*r_squared_list(0)))/(h*h);
  return 0;
}
