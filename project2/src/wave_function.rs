use crate::boltzmann::Boltzmann;
use ndarray::{s, Array1, Array2, Array3};

impl Boltzmann {
    pub fn evaluate(&self, particles: &Array2<f64>) -> f64 {
        let q_factor = self.compute_q_factor(particles);

        let mut term1 = 0.;
        let mut term2 = 1.;

        for n in 0..self.n {
            for d in 0..self.d {
                term1 += (particles[(n, d)] - self.a[(n, d)]).powf(2.);
            }
        }

        term1 = (-0.5 * term1 / (self.sigma * self.sigma)).exp();

        for h in 0..self.h {
            term2 *= 1. + q_factor[h].exp();
        }

        term1 * term2
    }

    pub fn quantum_force(&self, particles: &Array2<f64>) -> Array2<f64> {
        let q_factor = self.compute_q_factor(particles);

        let mut sum1 = Array2::<f64>::zeros((self.n, self.d));

        for h in 0..self.h {
            sum1 += &(&(self.w.slice(s!(.., .., h))) / (1. + (-q_factor[h]).exp()));
        }

        2. * (-(particles - &self.a) / (self.sigma * self.sigma) + sum1 / (self.sigma * self.sigma))
    }

    pub fn local_energy(&self, particles: &Array2<f64>) -> f64 {
        let q_factor = self.compute_q_factor(particles);

        let mut local_energy = 0.;

        for n in 0..self.n {
            for d in 0..self.d {
                let mut sum1 = 0.;
                let mut sum2 = 0.;

                for h in 0..self.h {
                    sum1 += self.w[(n, d, h)] / (1. + (-q_factor[h]).exp());
                    sum2 += self.w[(n, d, h)].powf(2.) * q_factor[h].exp()
                        / (1. + q_factor[h].exp()).powf(2.);
                }

                let derivative_ln_psi_1 = sum1 / (self.sigma * self.sigma)
                    - (particles[(n, d)] - self.a[(n, d)]) / (self.sigma * self.sigma);
                let derivative_ln_psi_2 = -1. / (self.sigma * self.sigma)
                    + sum2 / ((self.sigma * self.sigma) * (self.sigma * self.sigma));

                local_energy += 0.5
                    * (-derivative_ln_psi_1 * derivative_ln_psi_1 - derivative_ln_psi_2
                        + (particles[(n, d)]) * (particles[(n, d)]));
            }
        }

        if self.interaction {
            for particle_1 in 0..self.n {
                for particle_2 in particle_1 + 1..self.n {
                    let mut r_norm = 0.0;
                    for d in 0..self.d {
                        r_norm +=
                            (particles[(particle_1, d)] - particles[(particle_2, d)]).powf(2.);
                    }

                    local_energy += 1. / r_norm.sqrt();
                }
            }
        }

        local_energy
    }

    pub fn derivative(&self, particles: &Array2<f64>) -> (Array3<f64>, Array2<f64>, Array1<f64>) {
        let q_factor = self.compute_q_factor(particles);

        let mut dw = Array3::<f64>::zeros((self.n, self.d, self.h));

        for h in 0..self.h {
            dw.slice_mut(s!(.., .., h)).assign(
                &(&(self.w.slice(s!(.., .., h)))
                    / (self.sigma * self.sigma * (1. + (-q_factor[h]).exp()))),
            );
        }

        let da = (particles - &self.a) / (self.sigma * self.sigma);
        let db = 1. / (1. + (-q_factor).map(|x| x.exp()));

        return (dw, da, db);
    }

    fn compute_q_factor(&self, particles: &Array2<f64>) -> Array1<f64> {
        let mut q_factor = Array1::<f64>::zeros(self.h);

        for h in 0..self.h {
            q_factor[h] = (particles * &self.w.slice(s!(.., .., h))).sum();
        }

        &self.b + q_factor
    }
}
