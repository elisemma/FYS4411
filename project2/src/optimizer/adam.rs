use crate::boltzmann::Boltzmann;
use colored::Colorize;
use ndarray::{Array1, Array2, Array3};

impl Boltzmann {
    pub fn adam(&mut self) {
        let mut first_moment_w = Array3::<f64>::zeros((self.n, self.d, self.h));
        let mut first_moment_a = Array2::<f64>::zeros((self.n, self.d));
        let mut first_moment_b = Array1::<f64>::zeros(self.h);

        let mut second_moment_w = Array3::<f64>::zeros((self.n, self.d, self.h));
        let mut second_moment_a = Array2::<f64>::zeros((self.n, self.d));
        let mut second_moment_b = Array1::<f64>::zeros(self.h);

        // Values taked from pytorch
        let epsilon = 1e-8;
        let beta_1: f64 = 0.9;
        let beta_2: f64 = 0.999;

        for i in 0..self.max_optimizer_steps {
            let (energy, (e_der_a, e_der_b, e_der_w)) = self.monte_carlo();

            if self.verbose {
                println!(
                    "{}) Energy: {:.15}",
                    (i + 1).to_string().truecolor(100, 100, 100),
                    energy.to_string().truecolor(255, 183, 195)
                );
            }

            first_moment_a = (first_moment_a * beta_1) + (1. - beta_1) * &e_der_a;
            first_moment_b = (first_moment_b * beta_1) + (1. - beta_1) * &e_der_b;
            first_moment_w = (first_moment_w * beta_1) + (1. - beta_1) * &e_der_w;

            second_moment_a = (&second_moment_a * beta_2) + (1. - beta_2) * (&e_der_a * &e_der_a);
            second_moment_b = (&second_moment_b * beta_2) + (1. - beta_2) * (&e_der_b * &e_der_b);
            second_moment_w = (&second_moment_w * beta_2) + (1. - beta_2) * (&e_der_w * &e_der_w);

            let first_moment_bias_corrected_a = &first_moment_a / (1. - beta_1.powf(i as f64 + 1.));
            let first_moment_bias_corrected_b = &first_moment_b / (1. - beta_1.powf(i as f64 + 1.));
            let first_moment_bias_corrected_w = &first_moment_w / (1. - beta_1.powf(i as f64 + 1.));

            let second_moment_bias_corrected_a =
                &second_moment_a / (1. - beta_2.powf(i as f64 + 1.));
            let second_moment_bias_corrected_b =
                &second_moment_b / (1. - beta_2.powf(i as f64 + 1.));
            let second_moment_bias_corrected_w =
                &second_moment_w / (1. - beta_2.powf(i as f64 + 1.));

            let step_a = &(self.learning_rate * first_moment_bias_corrected_a
                / (second_moment_bias_corrected_a.map(|x| x.sqrt()) + epsilon));
            let step_b = &(self.learning_rate * first_moment_bias_corrected_b
                / (second_moment_bias_corrected_b.map(|x| x.sqrt()) + epsilon));
            let step_w = &(self.learning_rate * first_moment_bias_corrected_w
                / (second_moment_bias_corrected_w.map(|x| x.sqrt()) + epsilon));

            self.optimizer_step_sizes.push(
                step_a.map(|x| x.abs()).sum() / step_a.len() as f64
                    + step_b.map(|x| x.abs()).sum() / step_b.len() as f64
                    + step_w.map(|x| x.abs()).sum() / step_w.len() as f64,
            );

            self.a -= step_a;
            self.b -= step_b;
            self.w -= step_w;
        }
    }
}
