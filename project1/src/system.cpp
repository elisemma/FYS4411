#include <iostream>

#include "project1/constants.h"
#include "project1/system.h"
#include <iomanip>

#include "project1/wave_functions/simple_analytical.h"
#include "project1/wave_functions/simple_numerical.h"
#include "project1/wave_functions/interactive_elliptical.h"

using namespace std;

System::System(Parameters parameters)
    :   random_engine(parameters.num_particles),
        num_particles(parameters.num_particles),
        dimensions(parameters.dimensions),
        delta_t(parameters.delta_t),
        system_data_filename(parameters.energy_filename)
{
    initialize_particles(num_particles, dimensions);

    if (parameters.system_enum == SystemEnum::SIMPLE_ANALYTICAL || parameters.system_enum == SystemEnum::SIMPLE_ANALYTICAL_IMPORTANCE) {
        wave_function = new SimpleAnalytical(num_particles, dimensions, parameters.alpha);
    } else if (parameters.system_enum == SystemEnum::SIMPLE_NUMERICAL) {
        wave_function = new SimpleNumerical(num_particles, dimensions, parameters.alpha);
    } else if (parameters.system_enum == SystemEnum::INTERACTIVE) {
        wave_function = new InteractiveElliptical(num_particles, dimensions, parameters.alpha, constants::a);
    } else if (parameters.system_enum == SystemEnum::NO_JASTROW) {
        wave_function = new InteractiveElliptical(num_particles, dimensions, parameters.alpha, 0);
    }

    if (parameters.system_enum == SystemEnum::SIMPLE_ANALYTICAL_IMPORTANCE) {
        importance_sampling = true;
    } else {
        importance_sampling = false;
    }

    if (system_data_filename != "") {
        energies = new double[constants::number_of_metropolis_steps];
        one_body_densities = new double[constants::number_of_metropolis_steps * num_particles];
    }
}


System::~System() {
    delete wave_function;
    delete[] particles;

    if (system_data_filename != "") delete[] energies;

}

void System::initialize_particles(int num_particles, int dimensions) {
    particles = new double[num_particles * dimensions];
    for (int i = 0; i < num_particles * dimensions; i++) {
        particles[i] = random_engine.next_double(-1, 1);
    }
}

void System::write_system_data_to_file() {
    if (system_data_filename == "") {
        throw runtime_error("Unexpected error: calling `write_energy_to_file` while `energy_filename` is empty");
    } 

    cout << "Writing energy to file" << endl;
    ofstream energy_file(system_data_filename + "_energy.dat");
    energy_file << std::setprecision(15);
    for (int i = 0; i < constants::number_of_metropolis_steps; i++) {
        energy_file << energies[i] << endl;
    }    
    energy_file.close();
    cout << "Finished writing energy to file" << endl;

    cout << "Writing one body to file" << endl;
    ofstream one_body_file(system_data_filename + "_one_body.dat");
    one_body_file << std::setprecision(15);
    for (int i = 0; i < constants::number_of_metropolis_steps; i++) {
        for (int particle = 0; particle < num_particles; particle++) {
            one_body_file << one_body_densities[i * num_particles + particle] << " ";
        }
        one_body_file << endl;
    }
    cout << "Finished writing one body to file" << endl;
}

void System::run_simulation() {
    int equilibrationSteps = (int) (constants::equilibration * constants::number_of_metropolis_steps);

    for (int i = 0; i < equilibrationSteps; i++) {
        if (importance_sampling) importance_metropolis_step();
        else                     metropolis_step();
    }

    double local_energy = 0;

    // Sampling
    int number_of_accepted_metropolis_steps = 0;
    double cumulative_energy = 0, cumulative_energy_squared = 0, cumulative_trail_derivative = 0, cumulative_trail_derivative_times_energy = 0;

    for (int i = 0; i < constants::number_of_metropolis_steps; i++) {
        bool accepted_step;

        if (importance_sampling) accepted_step = importance_metropolis_step();
        else                     accepted_step = metropolis_step();

        if (accepted_step) {
            number_of_accepted_metropolis_steps++;
            local_energy = wave_function->compute_local_energy(particles);
        }

        if (system_data_filename != "") {
            energies[i] = local_energy;
            for (int particle = 0; particle < num_particles; particle++) {
                one_body_densities[i * num_particles + particle] = 
                    sqrt(
                        particles[particle * dimensions + 0] * particles[particle * dimensions + 0] +
                        particles[particle * dimensions + 1] * particles[particle * dimensions + 1] +
                        particles[particle * dimensions + 2] * particles[particle * dimensions + 2]
                    );
            }
        }

        cumulative_energy  += local_energy;
        cumulative_energy_squared  += local_energy*local_energy;

        double trial_derivative = wave_function->compute_nabla_alpha(particles);

        cumulative_trail_derivative += trial_derivative;
        cumulative_trail_derivative_times_energy += trial_derivative * local_energy;
    }

    accepted_ratio = (double) number_of_accepted_metropolis_steps / constants::number_of_metropolis_steps;
    energy_expectation = cumulative_energy / constants::number_of_metropolis_steps;
    energy_variance    = (cumulative_energy_squared / constants::number_of_metropolis_steps) - (energy_expectation * energy_expectation);

    // Term 1
    double expected_trail_energy_derivative = cumulative_trail_derivative_times_energy / constants::number_of_metropolis_steps;

    // Term 2 (part 1)
    double expected_trail_derivative = cumulative_trail_derivative / constants::number_of_metropolis_steps;

    // Final
    alpha_derivative = 2 * (expected_trail_energy_derivative - expected_trail_derivative * energy_expectation);
}


bool System::metropolis_step() {
    int particle = random_engine.random_particle();
    double movement[dimensions], wave_before = wave_function->evaluate(particles);

    for (int dim = 0; dim < dimensions; dim++) {
        movement[dim] = random_engine.next_double(-delta_t, delta_t);
        particles[particle * dimensions + dim] += movement[dim];
    }

    double wave_after = wave_function->evaluate(particles);

    double ratio = (wave_after * wave_after) / (wave_before * wave_before);

    if (random_engine.next_double() < ratio) return true;

    for (int dim = 0; dim < dimensions; dim++) particles[particle * dimensions + dim] -= movement[dim];

    return false;
}

bool System::importance_metropolis_step() {
    int particle = random_engine.random_particle();
    double movement[dimensions], wave_before = wave_function->evaluate(particles), quantum_before[dimensions];

    wave_function->compute_quantum_force(particles, particle, quantum_before);

    double D = 0.5;

    for (int dim = 0; dim < dimensions; dim++) {
        movement[dim] = D * quantum_before[dim] * delta_t + random_engine.next_normal(0, 1) * sqrt(delta_t);
        
        particles[particle * dimensions + dim] += movement[dim];
    }

    double wave_after = wave_function->evaluate(particles);

    double ratio = (wave_after * wave_after) / (wave_before * wave_before);

    double quantum_after[dimensions];
    wave_function->compute_quantum_force(particles, particle, quantum_after);

    double G_xy = 0, G_yx = 0, term1 = 0, term2 = 0;

    for (int dim = 0; dim < dimensions; dim++) {
        term1 = - movement[dim] - D * delta_t * quantum_after[dim];
        G_xy += term1 * term1;
        term2 =   movement[dim] - D * delta_t * quantum_before[dim];
        G_yx += term2 * term2;
    }

    G_xy = exp(-G_xy / (4. * D * delta_t));
    G_yx = exp(-G_yx / (4. * D * delta_t));

    ratio *= G_xy / G_yx;

    if (random_engine.next_double() < ratio) return true;

    for (int dim = 0; dim < dimensions; dim++) particles[particle * dimensions + dim] -= movement[dim];

    return false;
}
