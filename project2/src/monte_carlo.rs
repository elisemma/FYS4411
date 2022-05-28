use crate::args::Sampler;
use crate::boltzmann::Boltzmann;
use ndarray::{Array1, Array2, Array3};

impl Boltzmann {
    pub fn monte_carlo(&mut self) -> (f64, (Array2<f64>, Array1<f64>, Array3<f64>)) {
        let mut energy = 0.;

        let mut delta_psi_a = Array2::<f64>::zeros((self.n, self.d));
        let mut delta_psi_b = Array1::<f64>::zeros(self.h);
        let mut delta_psi_w = Array3::<f64>::zeros((self.n, self.d, self.h));

        let mut derivative_psi_e_a = Array2::<f64>::zeros((self.n, self.d));
        let mut derivative_psi_e_b = Array1::<f64>::zeros(self.h);
        let mut derivative_psi_e_w = Array3::<f64>::zeros((self.n, self.d, self.h));

        let mut wfold = self.evaluate(&self.position_old);

        let mut quantum_force_old = Array2::<f64>::zeros((self.n, self.d));

        if self.sampler == Sampler::ImportanceSampling {
            quantum_force_old = self.quantum_force(&self.position_old);
        }

        for _ in 0..self.number_of_monte_carlo_cycles {
            match self.sampler {
                Sampler::BruteForce => {
                    wfold = self.sample_brute_force(wfold);
                }
                Sampler::ImportanceSampling => {
                    (wfold, quantum_force_old) = self.importance_sample(wfold, quantum_force_old);
                }
            }

            let delta_e = self.local_energy(&self.position_old);

            self.energies.push(delta_e);

            let (der_psi_w, der_psi_a, der_psi_b) = self.derivative(&self.position_old);

            delta_psi_a += &der_psi_a;
            delta_psi_b += &der_psi_b;
            delta_psi_w += &der_psi_w;

            energy += delta_e;

            derivative_psi_e_a += &(der_psi_a * delta_e);
            derivative_psi_e_b += &(der_psi_b * delta_e);
            derivative_psi_e_w += &(der_psi_w * delta_e);

            if self.use_one_body {
                for particle in 0..self.n {
                    let one_body = (self.position_old[(particle, 0)].powf(2.)
                        + self.position_old[(particle, 1)].powf(2.))
                    .sqrt();
                    self.one_body_vec.push(one_body);
                }
            }
        }

        energy /= self.number_of_monte_carlo_cycles as f64;

        derivative_psi_e_a /= self.number_of_monte_carlo_cycles as f64;
        derivative_psi_e_b /= self.number_of_monte_carlo_cycles as f64;
        derivative_psi_e_w /= self.number_of_monte_carlo_cycles as f64;

        delta_psi_a /= self.number_of_monte_carlo_cycles as f64;
        delta_psi_b /= self.number_of_monte_carlo_cycles as f64;
        delta_psi_w /= self.number_of_monte_carlo_cycles as f64;

        let energy_der_a = 2. * (derivative_psi_e_a - delta_psi_a * energy);
        let energy_der_b = 2. * (derivative_psi_e_b - delta_psi_b * energy);
        let energy_der_w = 2. * (derivative_psi_e_w - delta_psi_w * energy);

        (energy, (energy_der_a, energy_der_b, energy_der_w))
    }
}
