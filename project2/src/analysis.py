import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import tikzplotlib
from blocking import block

plt.style.use("ggplot")

sns.set_theme()
sns.set_style("whitegrid")
sns.set_palette("bright")


class OptimizerSampler:
    def __init__(self, raw_name, name, short_name):
        self.raw_name = raw_name
        self.name = name
        self.short_name = short_name


OPTIMIZERS = [
    OptimizerSampler("gradient_descent", "gradient descent", "GD"),
    OptimizerSampler("adam", "adam", "Adam"),
]

SAMPLERS = [
    OptimizerSampler("brute_force", "brute force", "BF"),
    OptimizerSampler("importance_sampling", "importance sampling", "IS"),
]


def plot_all_energy_per_iteration():
    for interaction in [True, False]:
        no_interaction_txt = "" if interaction else "no_"
        #  for optimizer in reversed(OPTIMIZERS):
        for sampler in SAMPLERS:
            _, ax = plt.subplots(figsize=(4, 3))
            for optimizer in reversed(OPTIMIZERS):
                raw_data = np.loadtxt(
                    f"output/loss_function_{optimizer.raw_name}_{sampler.raw_name}_{no_interaction_txt}interaction.dat"
                )
                plt.plot(
                    raw_data,
                    #  ":" if sampler == SAMPLERS[1] else "-",
                    #  label=f"{optimizer.name.capitalize()} {sampler.name}",
                    label=optimizer.name.capitalize(),
                    alpha=1 if optimizer == OPTIMIZERS[1] else 0.5,
                    color="r" if optimizer == OPTIMIZERS[1] else "b",
                    #  linewidth=0.1,
                )

            ax.set_title(
                f"Energy {'with' if interaction else 'without'} interaction per\niteration for {sampler.name}",
            )
            ax.set_ylabel("Energy")
            ax.set_xlabel("Monte carlo iteration")
            if interaction:
                ax.set_yscale("log")

            ax.legend()
            plt.tight_layout()
            plt.savefig(
                f"output/loss_functions_all_{no_interaction_txt}interaction_{sampler.raw_name}.pdf"
            )
            plt.close()


def plot_energy_per_monte_carlo():
    for interaction in [True, False]:
        no_interaction_txt = "" if interaction else "no_"
        for optimizer in OPTIMIZERS:
            for sampler in SAMPLERS:
                raw_data = np.loadtxt(
                    f"output/loss_function_{optimizer.raw_name}_{sampler.raw_name}_{no_interaction_txt}interaction.dat"
                )
                raw_data -= 3 if interaction else 2
                monte_carlo_size = 10000
                percentages = int(len(raw_data) / monte_carlo_size)
                y = np.empty(percentages)
                for i in range(percentages):
                    y[i] = abs(
                        np.mean(
                            raw_data[monte_carlo_size * i : monte_carlo_size * (i + 1)]
                        )
                    )

                #  raw_data = np.abs(raw_data - (3 if interaction else 2))
                #  plt.plot(y, label=f"{optimizer.short_name} {sampler.short_name}")
                sns.lineplot(
                    x=range(1, len(y) + 1),
                    y=y,
                    label=f"{optimizer.short_name} {sampler.short_name}",
                    linestyle="--"
                    if sampler == SAMPLERS[1] and not interaction
                    else "-",
                    linewidth=2 if interaction else 3,
                )
        plt.title(
            rf"Absolute energy error for different optimizer\\steps {'with' if interaction else 'without'} interaction"
        )
        plt.legend()
        plt.ylabel("Energy")
        plt.xlabel("Optimizer interation")
        tikzplotlib.save(
            f"output/loss_functions_smooth_{no_interaction_txt}interaction.tex",
            extra_axis_parameters=[
                "title style={align=center}",
                "xmajorticks=true",
                "ymajorticks=true",
                "mark options={mark size=2.5pt, line width=1.5pt}",
            ],
            strict=True,
        )
        plt.close()


def plot_learning_rates():
    for interaction in [True, False]:
        no_interaction_txt = "" if interaction else "no_"
        df = pd.read_csv(
            f"output/learning_rates_{no_interaction_txt}interactions.tsv", sep="\t"
        )
        for optimizer in OPTIMIZERS:
            for sampler in SAMPLERS:
                small_df = df.query(
                    f'optimizer == "{optimizer.raw_name}" & sampler == "{sampler.raw_name}"'
                )
                small_df.energy = (small_df.energy - (3 if interaction else 2)).abs()
                sns.lineplot(
                    x="learning rate",
                    y="energy",
                    data=small_df,
                    label=f"{optimizer.short_name}, {sampler.short_name}",
                    linestyle="--"
                    if sampler == SAMPLERS[1]
                    and optimizer == OPTIMIZERS[0]
                    and not interaction
                    else "-",
                    linewidth=3,
                )
        plt.legend()

        plt.ylabel("Energy")
        plt.xlabel("Learning rate")

        plt.title(
            rf"Absolute energy error for different learning\\rates {'with' if interaction else 'without'} interactions"
        )
        tikzplotlib.save(
            f"output/learning_rates_{no_interaction_txt}interactions.tex",
            extra_axis_parameters=[
                "title style={align=center}",
                "xmajorticks=true",
                "ymajorticks=true",
                #  "mark options={mark size=2.5pt, line width=20pt}",
            ],
            strict=True,
        )
        plt.close()


def omega_sigma_square_learning_rate():
    for interaction in [True, False]:
        no_interaction_txt = "" if interaction else "no_"
        df = pd.read_csv(
            f"output/omega_sigma_square_{no_interaction_txt}interactions.tsv", sep="\t"
        )
        #  for sigma_square in [False]:
        for optimizer in OPTIMIZERS:
            for sampler in SAMPLERS:
                #  small_df = df[df["sigma square"] == sigma_square]
                small_df = df.query(
                    f'optimizer == "{optimizer.raw_name}" & sampler == "{sampler.raw_name}"'
                )
                small_df.energy = (small_df.energy - (3 if interaction else 2)).abs()
                #  sigma_txt = (
                #      r"$\sigma=\frac{1}{\omega}$" if sigma_square else r"$\sigma=1$"
                #  )
                sns.lineplot(
                    x="omega",
                    y="energy",
                    data=small_df,
                    linestyle="--"
                    if sampler == SAMPLERS[1] and optimizer == OPTIMIZERS[0]
                    else "-",
                    linewidth=3,
                    label=rf"{optimizer.short_name}, {sampler.short_name}",
                )
        plt.title(
            rf"Absolute energy error for different $\omega$\\{'with' if interaction else 'without'} interactions"
        )
        plt.legend()
        plt.ylabel("Energy")
        plt.xlabel(r"$\omega$")
        plt.yscale("log")

        tikzplotlib.save(
            f"output/omega_sigma_square_{no_interaction_txt}interactions.tex",
            extra_axis_parameters=[
                "title style={align=center}",
                "xmajorticks=true",
                "ymajorticks=true",
                "mark options={mark size=2.5pt, line width=1.5pt}",
            ],
            strict=True,
        )
        plt.close()


def hidden_layer_sizes():
    for interaction in [True, False]:
        no_interaction_txt = "" if interaction else "no_"
        df = pd.read_csv(
            f"output/hidden_layers_{no_interaction_txt}interactions.tsv", sep="\t"
        )
        df.energy = (df.energy - (3 if interaction else 2)).abs()
        #  for optimizer in ["adam", "gradient_descent"]:
        for optimizer in OPTIMIZERS:
            for sampler in SAMPLERS:
                #  for sampler in ["brute_force", "importance_sampling"]:
                small_df = df.query(
                    f'optimizer == "{optimizer.raw_name}" & sampler == "{sampler.raw_name}"'
                )
                sns.lineplot(
                    x="hidden layers",
                    y="energy",
                    data=small_df,
                    label=f"{optimizer.short_name}, {sampler.short_name}",  # , {sampler}
                    #  linestyle="--" if sampler == SAMPLERS[1] else "-",
                    linewidth=3,
                )
        plt.title(
            rf"Absolute energy error for different hidden\\layers {'with' if interaction else 'without'} interactions"
        )
        plt.legend()
        plt.ylabel("Energy")
        plt.xlabel("Hidden layers")

        tikzplotlib.save(
            f"output/hidden_layers_{no_interaction_txt}interactions.tex",
            extra_axis_parameters=[
                "title style={align=center}",
                "xmajorticks=true",
                "ymajorticks=true",
                "mark options={mark size=2.5pt, line width=1.5pt}",
            ],
            strict=True,
        )
        plt.close()


def plot_optimizer_step_sizes():
    for interaction in [False, True]:
        no_interaction_txt = "" if interaction else "no_"
        for optimizer in OPTIMIZERS:
            for sampler in SAMPLERS:

                filename = f"output/step_size_{'true' if interaction else 'false'}_{optimizer.raw_name}_{sampler.raw_name}.dat"
                y = np.loadtxt(filename)

                sns.lineplot(
                    x=range(1, len(y) + 1),
                    y=y,
                    label=f"{optimizer.short_name} {sampler.short_name}",
                    linestyle="--"
                    if sampler == SAMPLERS[1] and not interaction
                    else "-",
                    linewidth=2 if interaction else 3,
                )

        plt.title(
            f"Mean optimizer step size {'with' if interaction else 'without'} interactions"
        )
        plt.legend()
        plt.ylabel("Step size")
        plt.xlabel("Optimizer interation")
        #  plt.yscale("log")

        tikzplotlib.save(
            f"output/optimizer_step_size_{no_interaction_txt}interactions.tex",
            extra_axis_parameters=[
                "title style={align=center}",
                "xmajorticks=true",
                "ymajorticks=true",
                "mark options={mark size=2.5pt, line width=1.5pt}",
            ],
            strict=True,
        )
        #  plt.savefig(f"nice_{no_interaction_txt}.pdf")
        plt.close()


def plot_one_body_density():
    for interaction in [False, True]:
        for optimizer in OPTIMIZERS:
            for sampler in SAMPLERS:
                filename = f"output/one_body_{'true' if interaction else 'false'}_{optimizer.raw_name}_{sampler.raw_name}.dat"
                data = np.sort(np.loadtxt(filename))

                steps = data.shape[0]
                data = data.reshape(-1)
                H, bins = np.histogram(data, 50, range=(0, 2))
                bins = 0.5 * (bins[:-1] + bins[1:])
                sns.lineplot(
                    x=bins,
                    y=H / (bins ** 2 * 2 * steps),
                    label=f"{optimizer.short_name}, {sampler.short_name}, {'w int' if interaction else 'w/o int'}",
                    linewidth=3,
                    linestyle="--" if interaction else "-",
                )

    plt.xscale("log")
    plt.title(f"One body density")
    plt.legend()
    plt.ylabel(r"$\rho(\vec{r})$")
    plt.xlabel(r"$|\vec{r}|$")

    tikzplotlib.save(
        f"output/one_body.tex",
        extra_axis_parameters=[
            "title style={align=center}",
            "xmajorticks=true",
            "ymajorticks=true",
            "mark options={mark size=2.5pt, line width=1.5pt}",
        ],
        strict=True,
    )
    #  plt.savefig("nice.pdf")
    plt.close()


def create_table_values():
    print(r"\hline")
    print(
        r"\textbf{N} & \textbf{D} & \textbf{Inter} & $\boldsymbol{\omega}$ & \textbf{Opti} & \textbf{Spl} & \textbf{Abs err} & $\boldsymbol{\sigma_{block}}$ \\"
    )
    print(r"\hline")
    for optimizer in OPTIMIZERS:
        for sampler in SAMPLERS:
            table_1 = np.loadtxt(
                f"output/table_values_1_{optimizer.raw_name}_{sampler.raw_name}.dat"
            )[475712:]
            block_1 = block(table_1)
            print(
                f"${1}$ & ${1}$ & w/o & 1 & {optimizer.short_name} & {sampler.short_name} & {abs(np.mean(table_1) - 0.5):.4e} & {block_1:.4e} \\\\"
            )
            print(r"\hline")

    for optimizer in OPTIMIZERS:
        for sampler in SAMPLERS:
            table_2 = np.loadtxt(
                f"output/table_values_2_{optimizer.raw_name}_{sampler.raw_name}.dat"
            )[475712:]
            block_2 = block(table_2)
            print(
                f"${2}$ & ${2}$ & w/o & 1 & {optimizer.short_name} & {sampler.short_name} & {abs(np.mean(table_2) - 2):.4e} & {block_2:.4e} \\\\"
            )
            print(r"\hline")

    for optimizer in OPTIMIZERS:
        for sampler in SAMPLERS:
            table_3 = np.loadtxt(
                f"output/table_values_3_{optimizer.raw_name}_{sampler.raw_name}.dat"
            )[475712:]
            block_3 = block(table_3)
            print(
                f"${2}$ & ${2}$ & with & 0.75 &{optimizer.short_name} & {sampler.short_name} & {abs(np.mean(table_3) - 3):.4e} & {block_3:.4e} \\\\"
            )
            print(r"\hline")

    for optimizer in OPTIMIZERS:
        for sampler in SAMPLERS:
            table_4 = np.loadtxt(
                f"output/table_values_4_{optimizer.raw_name}_{sampler.raw_name}.dat"
            )[475712:]
            block_4 = block(table_4)
            print(
                f"${2}$ & ${2}$ & with & 1 & {optimizer.short_name} & {sampler.short_name} & {abs(np.mean(table_4) - 3):.4e} & {block_4:.4e} \\\\"
            )
            print(r"\hline")

    for optimizer in OPTIMIZERS:
        for sampler in SAMPLERS:
            table_5 = np.loadtxt(
                f"output/table_values_5_{optimizer.raw_name}_{sampler.raw_name}.dat"
            )[475712:]
            block_5 = block(table_5)
            print(
                f"${2}$ & ${2}$ & with & 1.25 & {optimizer.short_name} & {sampler.short_name} & {abs(np.mean(table_5) - 3):.4e} & {block_5:.4e} \\\\"
            )
            print(r"\hline")

    #  print(
    #      f"Value 1: mu: {np.mean(table_1)}, err: {abs(0.5 - np.mean(table_1))}, sd: {block_1}"
    #  )
    #  print(
    #      f"Value 2: mu: {np.mean(table_2)}, err: {abs(2 - np.mean(table_2))}, sd: {block_2}"
    #  )
    #  print(
    #      f"Value 3: mu: {np.mean(table_3)}, err: {abs(3 - np.mean(table_3))}, sd: {block_3}"
    #  )

    #  table_2 = np.loadtxt(f"output/table_values_2.dat")[475712:]
    #  table_3 = np.loadtxt(f"output/table_values_3.dat")[475712:]


if __name__ == "__main__":
    #  plot_learning_rates()
    #  omega_sigma_square_learning_rate()
    #  hidden_layer_sizes()
    #  plot_one_body_density()
    #  plot_all_energy_per_iteration()
    #  plot_energy_per_monte_carlo()
    #  plot_optimizer_step_sizes()
    create_table_values()
