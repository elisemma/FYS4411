#include <iomanip>
#include <iostream>
#include <chrono>
#include "project1/system_enum.h"
#include "project1/system.h"
#include "project1/constants.h"
#include "project1/parameters.h"

#define RED "\033[1;31m"
#define GREEN "\033[1;32m"
#define YELLOW "\033[1;33m"
#define BLUE "\033[1;34m"
#define MAGENTA "\033[1;35m"
#define CYAN "\033[1;36m"
#define BOLD "\033[1m"
#define RESET "\033[0m"

using namespace std;

void run_systems(vector<Parameters*> parameters_vec, string filename) {
    #pragma omp parallel for schedule(dynamic)
    for (Parameters *parameters : parameters_vec) {
        System system(*parameters);
        auto start = chrono::high_resolution_clock::now();
        system.run_simulation();
        auto end = chrono::high_resolution_clock::now();
        auto duration_ms = chrono::duration_cast<chrono::milliseconds>(end - start).count();

        parameters->energy_expectation = system.energy_expectation;
        parameters->energy_variance = system.energy_variance;
        parameters->accepted_ratio = system.accepted_ratio;
        parameters->time_ms = duration_ms;

        cout << BOLD << system_enum_to_string(parameters->system_enum) << RESET
            << " with " << BOLD << parameters->num_particles << RESET
            << " particles, " << BOLD << parameters->dimensions << RESET
            << " dim and " << round(parameters->alpha*100)/100
            << " alpha " << RED << "done" << RESET
            << " in " << CYAN << duration_ms << RESET <<
            " ms" << endl;

        if (parameters->energy_filename != "") {
            system.write_system_data_to_file(); 
        }
    }

    ofstream file(filename);
    write_heading(file);
    for (auto parameters : parameters_vec) {
        parameters->add_to_file(file);
    }
}

void gradient_alpha_search(vector<Parameters*> parameters_vec, double learning_rate, double delta) {
    #pragma omp parallel for schedule(dynamic)
    for (auto parameters : parameters_vec) {
        double alpha_change = 0;
        double alpha_0 = parameters->alpha;
        double specific_learning_rate = learning_rate / parameters->num_particles ;
        vector<double> alpha_vec;
        vector<double> energy_vec;
        int i = 0;

        do {
            System system(*parameters);
            system.run_simulation();

            alpha_vec.push_back(parameters->alpha);
            energy_vec.push_back(system.energy_expectation);

            cout << "Iter " << BOLD << (i+1) << RESET << ") " << BOLD << to_string(parameters->num_particles) << RESET << 
                " particles, a_0=" << BOLD << alpha_0 << RESET <<  ", alpha: " <<
                BOLD << parameters->alpha << RESET << ", change: " << BOLD <<
                specific_learning_rate*system.alpha_derivative << RESET <<
                endl;
            alpha_change = specific_learning_rate * system.alpha_derivative;
            parameters->alpha -= alpha_change;
            i++;
        } while (abs(alpha_change) > delta && i < 50);

        ofstream file("output/gradient/alpha_" + to_string(alpha_0) + "&system_enum=" + system_enum_to_string(parameters->system_enum) + "&particles=" + to_string(parameters->num_particles) + ".tsv");
        file << setprecision(15);
        file << "alpha\tenergy" << endl;
        for (size_t i = 0; i < alpha_vec.size(); i++) {
            file << alpha_vec[i] << "\t" << energy_vec[i] << endl;
        }
        file.close();
        cout << RED << "Finished" << RESET << " with " <<
            system_enum_to_string(parameters->system_enum) << " " <<
            to_string(parameters->num_particles) << " particles, a_0=" <<
            to_string(alpha_0) << " in " << (i+1) << " iteractions" << endl;
        cout << endl;
    }
}

void clear_parameter_vec(vector<Parameters*> &parameters_vec) {
    for (auto parameter : parameters_vec) {
        delete parameter;
    }

    parameters_vec.clear();
}

int main() {
    cout << setprecision(15);

    auto start_time = chrono::high_resolution_clock::now();
    int default_num_dimensions = 3;
    double default_delta_t = 0.5;

    double alphas[] = {0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.60, 0.65, 0.7};

    // Test accepted ratio for different delta_ts
    double delta_ts[] = {2, 1.5, 1, 0.75, 0.6, 0.55, 0.5, 0.45, 0.4, 0.25, 0.1, 0.05, 0.025};
    vector<Parameters*> parameters_vec;
    for (auto system_enum : {SIMPLE_ANALYTICAL, SIMPLE_ANALYTICAL_IMPORTANCE}) {
        double alpha = 0.5;
        for (auto delta_t : delta_ts) {
            parameters_vec.push_back(new Parameters(10, default_num_dimensions, alpha, delta_t, system_enum));
        }
    }
    run_systems(parameters_vec, "output/delta_t_acceptence_rate.tsv");

    // Gradient descent on alpha
    clear_parameter_vec(parameters_vec);
    for (auto num_particles : {1, 10 ,50, 100})
        for (double alpha : alphas)
            parameters_vec.push_back(new Parameters(num_particles, default_num_dimensions, alpha, default_delta_t, SIMPLE_ANALYTICAL));
    gradient_alpha_search(parameters_vec, 0.05, 1e-8);

    clear_parameter_vec(parameters_vec);
    for (auto num_particles : {1, 10 ,50, 100})
        for (double alpha : alphas)
            parameters_vec.push_back(new Parameters(num_particles, default_num_dimensions, alpha, default_delta_t, INTERACTIVE));
    gradient_alpha_search(parameters_vec, 0.05, 1e-5);

    clear_parameter_vec(parameters_vec);
    for (SystemEnum system_enum : {SIMPLE_ANALYTICAL, SIMPLE_NUMERICAL, SIMPLE_ANALYTICAL_IMPORTANCE}) {
        double delta_t;
        delta_t = 0.5;

        for (double num_particles : {1, 10, 100, 500})
            for (double dimensions : {1, 2, 3})
                for (double alpha : alphas)
                    parameters_vec.push_back(
                        new Parameters(num_particles, dimensions, alpha, delta_t, system_enum)
                    );
    }
    run_systems(parameters_vec, "output/simple.tsv");

    // Run interactive system
    clear_parameter_vec(parameters_vec);
    double interactive_dimensions = 3;
    for (double num_particles : {10, 50, 100}) {
        for (double alpha : alphas)
            parameters_vec.push_back(
                new Parameters(num_particles, interactive_dimensions, alpha, default_delta_t, INTERACTIVE)
            );
    }
    run_systems(parameters_vec, "output/interactive.tsv");

    // Run jastrow analysis
    clear_parameter_vec(parameters_vec);
    for (int num_particles : {10, 50, 100})
        for (SystemEnum system_enum : {INTERACTIVE, NO_JASTROW})
            parameters_vec.push_back(
                new Parameters(num_particles, interactive_dimensions, 0.5, default_delta_t, system_enum, "output/data/" + system_enum_to_string(system_enum) + "&particles=" + to_string(num_particles))
            );
    run_systems(parameters_vec, "output/jastrow.tsv");

    for (auto parameter : parameters_vec) {
        delete parameter;
    }

    auto stop_time = chrono::high_resolution_clock::now();
    auto big_duration_ms = chrono::duration_cast<chrono::nanoseconds>(stop_time - start_time).count();

    cout << "Total running time: " << big_duration_ms << " ns" << endl;

    return 0;
}
