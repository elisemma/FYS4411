use crate::args::{Optimizer, Sampler};
use colored::*;
use ndarray::{Array, Array1, Array2, Array3};
use ndarray_rand::rand_distr::Normal;
use ndarray_rand::RandomExt;
use rand::{Rng, SeedableRng};
use std::time::Instant;
use tinymt::{TinyMT64, TinyMT64Seed};

pub struct Boltzmann {
    pub n: usize,
    pub d: usize,
    pub h: usize,
    pub omega: f64,
    pub sigma: f64,
    pub learning_rate: f64,
    pub interaction: bool,
    pub optimizer: Optimizer,
    pub sampler: Sampler,
    pub w: Array3<f64>,
    pub a: Array2<f64>,
    pub b: Array1<f64>,
    pub number_of_monte_carlo_cycles: usize,
    pub delta_t: f64,
    pub delta: f64,
    pub max_optimizer_steps: usize,
    pub rng: TinyMT64,
    pub position_new: Array2<f64>,
    pub position_old: Array2<f64>,
    pub time_step: f64,
    pub diffusion: f64,
    pub energies: Vec<f64>,
    pub use_one_body: bool,
    pub one_body_vec: Vec<f64>,
    pub verbose: bool,
    pub optimizer_step_sizes: Vec<f64>,
}

impl Boltzmann {
    pub fn new(
        num_particles: usize,
        dimensions: usize,
        number_of_monte_carlo_cycles: usize,
        interaction: bool,
        omega: f64,
        sigma: f64,
        learning_rate: f64,
        num_hidden: usize,
        sampler: Sampler,
        optimizer: Optimizer,
        seed: u64,
        delta_t: f64,
        weight_standard_deviation: f64,
        delta: f64,
        max_optimizer_steps: usize,
        use_one_body: bool,
        verbose: bool,
    ) -> Self {
        let mut rng = TinyMT64::from_seed(TinyMT64Seed::from(seed));

        let a: Array2<f64> = Array::random_using(
            (num_particles, dimensions),
            Normal::new(0., weight_standard_deviation).unwrap(),
            &mut rng,
        );
        let b: Array1<f64> = Array::random_using(
            num_hidden,
            Normal::new(0., weight_standard_deviation).unwrap(),
            &mut rng,
        );
        let w: Array3<f64> = Array::random_using(
            (num_particles, dimensions, num_hidden),
            Normal::new(0., weight_standard_deviation).unwrap(),
            &mut rng,
        );

        let time_step: f64 = 0.05;

        let position_old: Array2<f64> = time_step.sqrt()
            * Array::random_using(
                (num_particles, dimensions),
                Normal::new(0., 1.).unwrap(),
                &mut rng,
            );
        let position_new = position_old.clone();

        let one_body_vec = Vec::new();

        Self {
            n: num_particles,
            d: dimensions,
            h: num_hidden,
            omega,
            sigma,
            learning_rate,
            interaction,
            optimizer,
            sampler,
            w,
            a,
            b,
            number_of_monte_carlo_cycles,
            delta_t,
            delta,
            max_optimizer_steps,
            rng,
            position_new,
            position_old,
            time_step,
            diffusion: 0.5,
            energies: Vec::new(),
            use_one_body,
            one_body_vec,
            verbose,
            optimizer_step_sizes: Vec::new(),
        }
    }

    pub fn run(&mut self) {
        if self.verbose {
            println!(
                "Running the system with {} particles and {} dimensions with ω={}, the {} optimizer and {} and {} interactions",
                self.n.to_string().truecolor(254, 144, 0),
                self.d.to_string().truecolor(254, 144, 0),
                self.omega.to_string().truecolor(254, 144, 0),
                self.optimizer.to_string().truecolor(254, 144, 0),
                self.sampler.to_string().truecolor(254, 144, 0),
                (if self.interaction {"with"} else {"without"}).truecolor(254, 144, 0),
            );
        }

        let now = Instant::now();

        match self.optimizer {
            Optimizer::GradientDescent => {
                self.gradient_descent();
            }
            Optimizer::Adam => {
                self.adam();
            }
        }

        if self.verbose {
            println!(
                "Used {} ms",
                now.elapsed()
                    .as_millis()
                    .to_string()
                    .truecolor(98, 187, 193)
            );
        }
    }

    pub fn rand_unif(&mut self, a: f64, b: f64) -> f64 {
        self.rng.gen_range(a..b)
    }

    pub fn rand_normal(&mut self, mu: f64, sigma_square: f64) -> f64 {
        Array::random_using(1, Normal::new(mu, sigma_square).unwrap(), &mut self.rng)[0]
    }
}
