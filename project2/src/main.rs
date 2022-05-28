use clap::Parser;
use colored::*;
use ndarray::{Array, Array1};
use project2::args::{Args, RunType, Optimizer, Sampler};
use project2::boltzmann::Boltzmann;
use std::fs::{create_dir_all, File};
use std::io::{Result, Write};
use strum::IntoEnumIterator;

fn get_generic_boltzmann(args: &Args) -> Boltzmann {
    let use_one_body = false;
    Boltzmann::new(
        args.num_particles,
        args.dimensions,
        args.monte_carlo_cycles,
        args.interaction,
        args.omega,
        1. / args.omega.sqrt(),
        args.learning_rate,
        args.hidden_layers,
        args.sampler,
        args.optimizer,
        args.seed,
        args.delta_t,
        args.weight_standard_deviation,
        args.delta,
        args.max_optimizer_steps,
        use_one_body,
        args.verbose,
    )
}

#[allow(dead_code)]
fn test_learning_rates() -> Result<()> {
    let args = Args::parse();
    println!(
        "{}",
        "Running system with different learning rates!".truecolor(130, 2, 99)
    );

    for interaction in [false, true] {
        let learning_rates: Array1<f64> = Array::linspace(0.00001, 0.5, 26);

        let mut f = File::create(format!(
            "output/learning_rates_{}interactions.tsv",
            if interaction { "" } else { "no_" }
        ))?;
        write!(f, "learning rate\toptimizer\tsampler\tenergy\n")?;

        println!(
            "Running {} interactions",
            (if interaction { "with" } else { "without" }).truecolor(107, 191, 89)
        );
        for optimizer in Optimizer::iter() {
            for sampler in Sampler::iter() {
                for learning_rate in &learning_rates {
                    let mut boltzmann = get_generic_boltzmann(&args);

                    boltzmann.learning_rate = *learning_rate;
                    boltzmann.optimizer = optimizer;
                    boltzmann.sampler = sampler;
                    boltzmann.interaction = interaction;

                    boltzmann.run();

                    write!(
                        f,
                        "{}\t{}\t{}\t{:.15}\n",
                        learning_rate,
                        optimizer.to_string().replace(" ", "_"),
                        sampler.to_string().replace(" ", "_"),
                        boltzmann.energies.iter().sum::<f64>() / boltzmann.energies.len() as f64
                    )?;

                    println!(
                        "Learning rate {:.7}:\tOptimizer: {}\tSampler: {}\t{}: {}",
                        learning_rate.to_string().truecolor(98, 187, 193),
                        optimizer.to_string().truecolor(98, 187, 193),
                        sampler.to_string().truecolor(98, 187, 193),
                        "〈E〉".truecolor(254, 144, 0),
                        (boltzmann.energies.iter().sum::<f64>() / boltzmann.energies.len() as f64)
                            .to_string()
                            .truecolor(205, 255, 26),
                    );
                }
            }
        }
    }

    Ok(())
}

#[allow(dead_code)]
fn test_omegas_different_sigmas() -> Result<()> {
    let mut args = Args::parse();
    println!(
        "{}",
        "Running system with different omegas and sigmas!".truecolor(130, 2, 99)
    );
    let old_omega = args.omega;

    let omegas: Array1<f64> = Array::linspace(0.25, 3., 21);
    for interaction in [true, false] {
        let mut f = File::create(format!(
            "output/omega_sigma_square_{}interactions.tsv",
            if interaction { "" } else { "no_" }
        ))?;
        write!(f, "omega\toptimizer\tsampler\tenergy\n")?;

        println!(
            "Running {} interactions",
            (if interaction { "with" } else { "without" }).truecolor(107, 191, 89)
        );

        for optimizer in Optimizer::iter() {
            for sampler in Sampler::iter() {
                for omega in &omegas {
                    args.omega = *omega;

                    let mut boltzmann = get_generic_boltzmann(&args);

                    boltzmann.interaction = interaction;
                    boltzmann.optimizer = optimizer;
                    boltzmann.sampler = sampler;

                    boltzmann.run();

                    write!(
                        f,
                        "{}\t{}\t{}\t{:.15}\n",
                        omega,
                        optimizer.to_string().replace(" ", "_"),
                        sampler.to_string().replace(" ", "_"),
                        boltzmann.energies.iter().sum::<f64>() / boltzmann.energies.len() as f64
                    )?;

                    println!(
                        "Omega {:.4}\topt: {}\tsampler: {}\t{}: {}",
                        omega.to_string().truecolor(98, 187, 193),
                        optimizer.to_string().truecolor(98, 187, 193),
                        sampler.to_string().truecolor(98, 187, 193),
                        "〈E〉".truecolor(254, 144, 0),
                        (boltzmann.energies.iter().sum::<f64>() / boltzmann.energies.len() as f64)
                            .to_string()
                            .truecolor(205, 255, 26),
                    );
                }
            }
        }
    }

    args.omega = old_omega;

    Ok(())
}

#[allow(dead_code)]
fn test_num_hidden_layers() -> Result<()> {
    let mut args = Args::parse();
    println!(
        "{}",
        "Running system with different hidden layer sizes!".truecolor(130, 2, 99)
    );

    for interaction in [true, false] {
        let mut f = File::create(format!(
            "output/hidden_layers_{}interactions.tsv",
            if interaction { "" } else { "no_" }
        ))?;
        write!(f, "hidden layers\toptimizer\tsampler\tenergy\n")?;

        println!(
            "Running {} interactions",
            (if interaction { "with" } else { "without" }).truecolor(107, 191, 89)
        );

        for optimizer in Optimizer::iter() {
            for sampler in Sampler::iter() {
                for hidden_layers in 1..=8 {
                    args.hidden_layers = hidden_layers;
                    let mut boltzmann = get_generic_boltzmann(&args);

                    boltzmann.interaction = interaction;
                    boltzmann.optimizer = optimizer;
                    boltzmann.sampler = sampler;

                    boltzmann.run();

                    write!(
                        f,
                        "{}\t{}\t{}\t{:.15}\n",
                        hidden_layers,
                        optimizer.to_string().replace(" ", "_"),
                        sampler.to_string().replace(" ", "_"),
                        boltzmann.energies.iter().sum::<f64>() / boltzmann.energies.len() as f64
                    )?;

                    println!(
                        "Hidden layers {}:\topt: {}\tsampler: {}\t{}: {}",
                        hidden_layers,
                        optimizer.to_string().truecolor(98, 187, 193),
                        sampler.to_string().truecolor(98, 187, 193),
                        "〈E〉".truecolor(254, 144, 0),
                        (boltzmann.energies.iter().sum::<f64>() / boltzmann.energies.len() as f64)
                            .to_string()
                            .truecolor(98, 187, 193),
                    );
                }
            }
        }
    }
    args.hidden_layers = 2;

    Ok(())
}

#[allow(dead_code)]
fn create_loss_functions() -> Result<()> {
    let args = Args::parse();
    println!("{}", "Creating different loss plots".truecolor(130, 2, 99));
    for interaction in [true, false] {
        println!(
            "Running {} interactions",
            (if interaction { "with" } else { "without" }).truecolor(107, 191, 89)
        );

        for optimizer in Optimizer::iter() {
            for sampler in Sampler::iter() {
                let mut boltzmann = get_generic_boltzmann(&args);

                boltzmann.interaction = interaction;
                boltzmann.optimizer = optimizer;
                boltzmann.sampler = sampler;

                boltzmann.run();

                println!(
                    "Optimizer: {}\tsampler: {}\t{}: {}",
                    optimizer.to_string().truecolor(98, 187, 193),
                    sampler.to_string().truecolor(98, 187, 193),
                    "〈E〉".truecolor(254, 144, 0),
                    (boltzmann.energies.iter().sum::<f64>() / boltzmann.energies.len() as f64)
                        .to_string()
                        .truecolor(98, 187, 193),
                );

                let mut f = File::create(format!(
                    "output/loss_function_{}_{}_{}interaction.dat",
                    optimizer.to_string().replace(" ", "_"),
                    sampler.to_string().replace(" ", "_"),
                    if interaction { "" } else { "no_" }
                ))?;
                for energy in &boltzmann.energies {
                    write!(f, "{:.15}\n", energy)?;
                }
            }
        }
    }

    Ok(())
}

#[allow(dead_code)]
fn create_tables() -> Result<()> {
    let mut args = Args::parse();
    println!("{}", "Creating tables".truecolor(130, 2, 99));
    for optimizer in Optimizer::iter() {
        for sampler in Sampler::iter() {
            args.optimizer = optimizer;
            args.sampler = sampler;

            args.num_particles = 1;
            args.dimensions = 1;
            let mut boltzmann = get_generic_boltzmann(&args);
            boltzmann.interaction = false;
            boltzmann.run();
            let mut f = File::create(format!(
                "output/table_values_1_{}_{}.dat",
                optimizer.to_string().replace(" ", "_"),
                sampler.to_string().replace(" ", "_")
            ))?;
            for energy in &boltzmann.energies {
                write!(f, "{:.15}\n", energy)?;
            }
            println!("Wrote table 1");

            args.num_particles = 2;
            args.dimensions = 2;
            let mut boltzmann = get_generic_boltzmann(&args);
            boltzmann.interaction = false;
            boltzmann.run();
            let mut f = File::create(format!(
                "output/table_values_2_{}_{}.dat",
                optimizer.to_string().replace(" ", "_"),
                sampler.to_string().replace(" ", "_")
            ))?;
            for energy in &boltzmann.energies {
                write!(f, "{:.15}\n", energy)?;
            }
            println!("Wrote table 2");


            args.omega = 0.75;
            let mut boltzmann = get_generic_boltzmann(&args);
            boltzmann.interaction = true;
            boltzmann.run();
            let mut f = File::create(format!(
                "output/table_values_3_{}_{}.dat",
                optimizer.to_string().replace(" ", "_"),
                sampler.to_string().replace(" ", "_")
            ))?;
            for energy in &boltzmann.energies {
                write!(f, "{:.15}\n", energy)?;
            }
            println!("Wrote table 3");
            args.omega = 1.;

            let mut boltzmann = get_generic_boltzmann(&args);
            boltzmann.interaction = true;
            boltzmann.run();
            let mut f = File::create(format!(
                "output/table_values_4_{}_{}.dat",
                optimizer.to_string().replace(" ", "_"),
                sampler.to_string().replace(" ", "_")
            ))?;
            for energy in &boltzmann.energies {
                write!(f, "{:.15}\n", energy)?;
            }
            println!("Wrote table 4");
        }
    }

    Ok(())
}

#[allow(dead_code)]
fn test_one_body_density() -> Result<()> {
    let mut args = Args::parse();
    println!("{}", "Running one body".truecolor(130, 2, 99));
    for interactions in [false, true] {
        for optimizer in Optimizer::iter() {
            for sampler in Sampler::iter() {
                args.optimizer = optimizer;
                args.sampler = sampler;

                let use_one_body = true;
                let mut boltzmann = Boltzmann::new(
                    args.num_particles,
                    args.dimensions,
                    args.monte_carlo_cycles,
                    interactions,
                    args.omega,
                    1. / args.omega.sqrt(),
                    args.learning_rate,
                    args.hidden_layers,
                    args.sampler,
                    args.optimizer,
                    args.seed,
                    args.delta_t,
                    args.weight_standard_deviation,
                    args.delta,
                    args.max_optimizer_steps,
                    use_one_body,
                    args.verbose,
                );

                boltzmann.run();

                println!(
                    "Optimizer: {}\tsampler: {}\t{}: {}",
                    optimizer.to_string().truecolor(98, 187, 193),
                    sampler.to_string().truecolor(98, 187, 193),
                    "〈E〉".truecolor(254, 144, 0),
                    (boltzmann.energies.iter().sum::<f64>() / boltzmann.energies.len() as f64)
                        .to_string()
                        .truecolor(98, 187, 193),
                );


                let mut f = File::create(format!(
                    "output/one_body_{}_{}_{}.dat",
                    interactions,
                    optimizer.to_string().replace(" ", "_"),
                    sampler.to_string().replace(" ", "_")
                ))?;

                for one_body in &boltzmann.one_body_vec {
                    write!(f, "{:.15}\n", one_body)?;
                }
            }
        }
    }

    Ok(())
}

#[allow(dead_code)]
fn test_optimizer_step_sizes() -> Result<()> {
    let mut args = Args::parse();
    println!("{}", "Running step sizes".truecolor(130, 2, 99));
    for interactions in [false, true] {
        for optimizer in Optimizer::iter() {
            for sampler in Sampler::iter() {
                args.optimizer = optimizer;
                args.sampler = sampler;

                let mut boltzmann = Boltzmann::new(
                    args.num_particles,
                    args.dimensions,
                    args.monte_carlo_cycles,
                    interactions,
                    args.omega,
                    1. / args.omega.sqrt(),
                    args.learning_rate,
                    args.hidden_layers,
                    args.sampler,
                    args.optimizer,
                    args.seed,
                    args.delta_t,
                    args.weight_standard_deviation,
                    args.delta,
                    args.max_optimizer_steps,
                    args.one_body_density,
                    args.verbose,
                );

                boltzmann.run();

                println!(
                    "Optimizer: {}\tsampler: {}\t{}: {}",
                    optimizer.to_string().truecolor(98, 187, 193),
                    sampler.to_string().truecolor(98, 187, 193),
                    "〈E〉".truecolor(254, 144, 0),
                    (boltzmann.energies.iter().sum::<f64>() / boltzmann.energies.len() as f64)
                        .to_string()
                        .truecolor(98, 187, 193),
                );


                let mut f = File::create(format!(
                    "output/step_size_{}_{}_{}.dat",
                    interactions,
                    optimizer.to_string().replace(" ", "_"),
                    sampler.to_string().replace(" ", "_")
                ))?;

                for one_body in &boltzmann.optimizer_step_sizes {
                    write!(f, "{:.15}\n", one_body)?;
                }
            }
        }
    }

    Ok(())
}

#[allow(dead_code)]
fn run_raw() {
    let args = Args::parse();

    let mut boltzmann = Boltzmann::new(
        args.num_particles,
        args.dimensions,
        args.monte_carlo_cycles,
        args.interaction,
        args.omega,
        1. / args.omega.sqrt(),
        args.learning_rate,
        args.hidden_layers,
        args.sampler,
        args.optimizer,
        args.seed,
        args.delta_t,
        args.weight_standard_deviation,
        args.delta,
        args.max_optimizer_steps,
        args.one_body_density,
        args.verbose,
    );

    boltzmann.run();
}

fn main() -> Result<()> {
    create_dir_all("output")?;

    let args = Args::parse();
    match args.run_type {
        RunType::Normal => {
            run_raw();
        },
        RunType::Analysis => {
            test_learning_rates()?;
            test_omegas_different_sigmas()?;
            test_num_hidden_layers()?;
            create_tables()?;
            create_loss_functions()?;
            test_one_body_density()?;
            test_optimizer_step_sizes()?;
        }
    }
    

    Ok(())
}
