use crate::boltzmann::Boltzmann;
use colored::Colorize;

impl Boltzmann {
    pub fn gradient_descent(&mut self) {
        for i in 0..self.max_optimizer_steps {
            let (energy, (e_der_a, e_der_b, e_der_w)) = self.monte_carlo();

            if self.verbose {
                println!(
                    "{}) Energy: {:.15}",
                    (i + 1).to_string().truecolor(100, 100, 100),
                    energy.to_string().truecolor(255, 183, 195)
                );
            }

            self.optimizer_step_sizes.push(
                e_der_a.map(|x| x.abs()).sum() / e_der_a.len() as f64
                    + e_der_b.map(|x| x.abs()).sum() / e_der_b.len() as f64
                    + e_der_w.map(|x| x.abs()).sum() / e_der_w.len() as f64,
            );

            self.a -= &(self.learning_rate * e_der_a);
            self.b -= &(self.learning_rate * e_der_b);
            self.w -= &(self.learning_rate * e_der_w);
        }
    }
}
