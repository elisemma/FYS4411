use clap::{ArgEnum, Parser};
use strum::EnumIter;

/// A program that applies restricted boltzmann machine to the quantum many body problem
#[derive(Parser, Debug)]
#[clap(author="Herman Brunborg <herman@brunborg.com>", version, about, long_about = None)]
pub struct Args {
    /// What to run
    #[clap(short, long, arg_enum, default_value = "normal")]
    pub run_type: RunType,

    /// The optimizer method
    #[clap(short, long, arg_enum, default_value = "gradient-descent")]
    pub optimizer: Optimizer,

    /// The sampler method
    #[clap(short, long, arg_enum, default_value = "brute-force")]
    pub sampler: Sampler,

    /// Number of dimensions in the simulation
    #[clap(short, long, default_value = "2")]
    pub dimensions: usize,

    #[clap(short, long, default_value = "2")]
    /// Number of particles in the simulation
    pub num_particles: usize,

    /// Whether to include one body density or not
    #[clap(short, long)]
    pub one_body_density: bool,

    /// Whether to include interactions or not
    #[clap(short, long)]
    pub interaction: bool,

    /// Number of monte carlo cycles to run
    #[clap(long, default_value = "10000")]
    pub monte_carlo_cycles: usize,

    /// Omega
    #[clap(long, default_value = "1")]
    pub omega: f64,

    /// Learning rate
    #[clap(long, default_value = "0.2")]
    pub learning_rate: f64,

    /// Number of hidden layers
    #[clap(long, default_value = "2")]
    pub hidden_layers: usize,

    /// Delta t
    #[clap(long, default_value = "0.5")]
    pub delta_t: f64,

    /// Standard deviation for the weight initialization
    #[clap(long, default_value = "0.25")]
    pub weight_standard_deviation: f64,

    /// Random seed to use
    #[clap(long, default_value = "42")]
    pub seed: u64,

    /// Tolerance for the optimizer
    #[clap(long, default_value = "0.000000001")]
    pub delta: f64,

    /// Maximum number of optimizer steps
    #[clap(long, default_value = "100")]
    pub max_optimizer_steps: usize,

    /// Wether to include many prints or not
    #[clap(long)]
    pub verbose: bool,
}

#[derive(ArgEnum, Debug, Clone, Copy, PartialEq, EnumIter)]
pub enum RunType {
    /// Run the system with the given config
    Normal,

    /// Run all analysis with (with custom parameters)
    Analysis,
}

#[derive(ArgEnum, Debug, Clone, Copy, PartialEq, EnumIter)]
pub enum Optimizer {
    /// Standard gradient descent
    GradientDescent,

    /// Adam
    Adam,
}

impl Optimizer {
    pub fn to_string(&self) -> &str {
        match self {
            Optimizer::GradientDescent => "gradient descent",
            Optimizer::Adam => "adam",
        }
    }
}

#[derive(ArgEnum, Debug, Clone, Copy, PartialEq, EnumIter)]
pub enum Sampler {
    /// Standard gradient descent
    BruteForce,

    /// Adam
    ImportanceSampling,
}

impl Sampler {
    pub fn to_string(&self) -> &str {
        match self {
            Sampler::BruteForce => "brute force",
            Sampler::ImportanceSampling => "importance sampling",
        }
    }
}
