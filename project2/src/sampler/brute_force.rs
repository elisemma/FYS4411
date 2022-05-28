use crate::boltzmann::Boltzmann;
use rand::Rng;

impl Boltzmann {
    pub fn sample_brute_force(&mut self, mut wfold: f64) -> f64 {
        let particle = self.rng.gen_range(0..self.n);

        for d in 0..self.d {
            self.position_new[(particle, d)] =
                self.position_old[(particle, d)] + self.rand_unif(-self.delta_t, self.delta_t)
        }

        let wfnew = self.evaluate(&self.position_new);

        let ratio = wfnew.powf(2.) / wfold.powf(2.);

        if self.rand_unif(0., 1.) <= ratio {
            for j in 0..self.d {
                self.position_old[(particle, j)] = self.position_new[(particle, j)];
            }
            wfold = wfnew;
        }

        wfold
    }
}
