use crate::boltzmann::Boltzmann;
use ndarray::Array2;
use rand::Rng;

impl Boltzmann {
    pub fn importance_sample(
        &mut self,
        mut wfold: f64,
        mut quantum_force_old: Array2<f64>,
    ) -> (f64, Array2<f64>) {
        let particle = self.rng.gen_range(0..self.n);

        for d in 0..self.d {
            self.position_new[(particle, d)] = self.position_old[(particle, d)]
                + self.rand_normal(0., 1.) * self.time_step.sqrt()
                + quantum_force_old[(particle, d)] * self.time_step * self.diffusion;
        }

        let wfnew = self.evaluate(&self.position_new);
        let quantum_force_new = self.quantum_force(&self.position_new);

        let mut green_function = 0.;
        for j in 0..self.d {
            green_function += 0.5
                * (quantum_force_old[(particle, j)] + quantum_force_new[(particle, j)])
                * (self.diffusion
                    * self.time_step
                    * 0.5
                    * (quantum_force_old[(particle, j)] - quantum_force_new[(particle, j)])
                    - self.position_new[(particle, j)]
                    + self.position_old[(particle, j)]);
        }

        green_function = green_function.exp();

        let ratio = green_function * wfnew.powf(2.) / wfold.powf(2.);

        if self.rand_unif(0., 1.) <= ratio {
            for d in 0..self.d {
                self.position_old[(particle, d)] = self.position_new[(particle, d)];
                quantum_force_old[(particle, d)] = quantum_force_new[(particle, d)];
            }
            wfold = wfnew;
        }

        (wfold, quantum_force_old)
    }
}
