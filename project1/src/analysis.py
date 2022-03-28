# TODO: Do we need to compare with the exact results?
import argparse
import os
import re


from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import tikzplotlib

from blocking import block

plt.style.use("ggplot")

sns.set_theme()
sns.set_style("whitegrid")
sns.set_palette("bright")


class SystemEnum:
    def __init__(self, value, short_name, long_name, marker="o"):
        self._value = value
        self.short = short_name
        self.long = long_name
        self.marker = marker

    def __call__(self):
        return self._value


def particles_str(particles):
    return f"{particles} particle{'s' if particles > 1 else ''}"


def dimensions_str(dimensions):
    return f"{dimensions} dimension{'s' if dimensions > 1 else ''}"


SIMPLE_ANALYTICAL_BRUTE_FORCE = SystemEnum(
    "simple_analytical_brute_force", "Analytical", "Analytical brute force"
)
SIMPLE_NUMERICAL_BRUTE_FORCE = SystemEnum(
    "simple_numerical_brute_force", "Numerical", "Numerical brute force", marker="X"
)
SIMPLE_IMPORANCE_SAMPLING = SystemEnum(
    "simple_analytical_importance", "Importance", "Importance sampling"
)

INTERACTIVE = SystemEnum("interactive", "Interactive", "Interactive")


def plot_cpu_time():
    df = pd.read_csv("output/simple.tsv", sep="\t")

    for system_enum in [
        SIMPLE_ANALYTICAL_BRUTE_FORCE,
        SIMPLE_NUMERICAL_BRUTE_FORCE,
    ]:
        for dimensions in [1, 2, 3]:
            queried_df = df.query(
                f'system_enum=="{system_enum()}" & alpha==0.5 & dimensions=={dimensions}'
            )
            sns.lineplot(
                x="num_particles",
                y="time_ms",
                data=quescatterplotried_df,
                label=f"{system_enum.short} {dimensions_str(dimensions)}",
            )

    plt.yscale("log")
    plt.xscale("log")
    plt.legend()
    plt.xlabel(r"Number of particles")
    plt.ylabel(r"CPU time")
    plt.title(r"CPU time for analytical vs numerical brute force on a log log scale")

    tikzplotlib.save(
        "output/plots/cpu_time_analytical_vs_numerical.tex",
        extra_axis_parameters=[
            "title style={align=center}",
            "xmajorticks=true",
            "ymajorticks=true",
            "mark options={mark size=2.5pt, line width=1.5pt}",
        ],
        strict=True,
    )
    plt.show()
    plt.close()


def plot_gradient_descent():
    for system_enum in [SIMPLE_ANALYTICAL_BRUTE_FORCE, INTERACTIVE]:
        num_particles_list = [1, 10, 50, 100]
        for j, num_particles in enumerate(num_particles_list):
            colors = ["#023eff", "#ff7c00", "#1ac938", "#e8000b"]
            alphas = []
            for i, filename in enumerate(
                sorted(
                    filter(
                        lambda filename: system_enum() in filename
                        and f"particles={num_particles}." in filename,
                        os.listdir("output/gradient"),
                    )
                )
            ):
                df = pd.read_csv("output/gradient/" + filename, sep="\t")

                alphas.append(df["alpha"].iloc[-1])
                if system_enum == INTERACTIVE:
                    df = df[:10]
                else:
                    df = df[:15]

                if i == 0:
                    sns.lineplot(
                        x=df.index,
                        y="alpha",
                        data=df,
                        color=colors[j],
                    )
                else:
                    sns.lineplot(
                        x=df.index,
                        y="alpha",
                        data=df,
                        color=colors[j],
                    )

            plt.axhline(
                y=sum(alphas) / len(alphas),
                color=colors[j],
                linestyle="--",
                label=rf"{particles_str(num_particles)}, $\hat{{\alpha}}\approx{sum(alphas) / len(alphas):.4f}$",
            )
        plt.legend()

        plt.title(
            rf"Gradient search for $\alpha$ for different \\ initial $\alpha_0$ for the the {system_enum.short} case"
        )
        plt.legend()
        plt.xlabel("Iteration")
        plt.ylabel(r"$\alpha$")

        tikzplotlib.save(
            f"output/plots/gradient_descent_{system_enum()}.tex",
            extra_axis_parameters=[
                "title style={align=center}",
                "xmajorticks=true",
                "ymajorticks=true",
            ],
            strict=True,
        )

        plt.close()
        continue


def plot_simple_energy():
    df = pd.read_csv("output/simple.tsv", sep="\t")
    for dimensions in [1, 2, 3]:
        for num_particles in [1, 10, 100, 500]:
            for system_enum in [
                SIMPLE_ANALYTICAL_BRUTE_FORCE,
                SIMPLE_IMPORANCE_SAMPLING,
                SIMPLE_NUMERICAL_BRUTE_FORCE,
            ]:
                queried_df = df.query(
                    f'system_enum=="{system_enum()}" & num_particles=={num_particles} & dimensions=={dimensions}'
                )

                sns.scatterplot(
                    x="alpha",
                    y="energy_expectation",
                    data=queried_df,
                    label=f"{system_enum.long}",
                    marker=system_enum.marker,
                    #  linewidth=2,
                )

            plt.title(
                rf"Energy as a function of $\alpha$ for {dimensions_str(dimensions)} and"
                + r"\\"
                + f"{particles_str(num_particles)} in the spherical trap without interactions"
            )
            plt.xlabel(r"$\alpha$")
            plt.ylabel(r"$\langle E_L \rangle$")
            #  plt.legend()
            #  plt.grid(True)
            tikzplotlib.save(
                f"output/plots/energy_particles={num_particles}_dimensions={dimensions}.tex",
                extra_axis_parameters=[
                    "title style={align=center}",
                    "xmajorticks=true",
                    "ymajorticks=true",
                    "mark options={mark size=2.5pt, line width=1.5pt}",
                ],
                strict=True,
            )
            plt.close()


def plot_interactive_energy():
    df = pd.read_csv("output/interactive.tsv", sep="\t")
    for num_particles in [10, 50, 100]:
        queried_df = df.query(f"num_particles=={num_particles}")

        sns.scatterplot(
            x="alpha",
            y="energy_expectation",
            data=queried_df[:-1],
        )

        plt.title(
            rf"Energy as a function of $\alpha$ for {dimensions_str(3)} and"
            + r"\\"
            + f"{particles_str(num_particles)} in the elliptical trap with interactions"
        )
        plt.xlabel(r"$\alpha$")
        plt.ylabel(r"$\langle E_L \rangle$")
        plt.legend()

        tikzplotlib.save(
            f"output/plots/energy_interactive_particles={num_particles}.tex",
            extra_axis_parameters=[
                "title style={align=center}",
                "xmajorticks=true",
                "ymajorticks=true",
                "mark options={mark size=4pt}",
            ],
            strict=True,
        )
        plt.close()


def plot_energy():
    plot_simple_energy()
    plot_interactive_energy()


def plot_one_body():
    for num_particles in [10, 50, 100]:
        for i, filename in enumerate(
            [
                f"output/data/interactive&particles={num_particles}_one_body.dat",
                f"output/data/no_jastrow&particles={num_particles}_one_body.dat",
            ]
        ):
            data = np.loadtxt(filename)
            steps = data.shape[0]
            data = data.reshape(-1)
            H, bins = np.histogram(data, 30, range=(0, 2 * data.max() / 3))
            bins = 0.5 * (bins[:-1] + bins[1:])
            sns.lineplot(
                x=bins,
                y=H / (steps * num_particles * (bins ** 2)),
                label="With interactions" if i == 0 else "Without interactions",
            )

        plt.title(
            f"One body density for {particles_str(num_particles)} in the elliptical trap"
        )
        plt.legend()
        plt.ylabel(r"$\rho(\vec{r})$")
        plt.xlabel(r"$|\vec{r}|$")

        tikzplotlib.save(
            f"output/plots/one_body_{num_particles}_particles.tex",
            extra_axis_parameters=[
                "title style={align=center}",
                "xmajorticks=true",
                "ymajorticks=true",
                "mark options={mark size=2.5pt, line width=1.5pt}",
            ],
            strict=True,
        )

        plt.close()


def plot_delta_t():
    df = pd.read_csv("output/delta_t_acceptence_rate.tsv", sep="\t")
    for system_enum in [
        SIMPLE_ANALYTICAL_BRUTE_FORCE,
        SIMPLE_IMPORANCE_SAMPLING,
    ]:
        queried_df = df.query(f"system_enum=='{system_enum()}'")

        sns.scatterplot(
            x="delta_t",
            y="accepted_ratio",
            data=queried_df,
            label=f"{system_enum.long}",
            marker=system_enum.marker,
            linewidth=0,
        )

    plt.title(
        rf"Acceptence rate for different step lengths $\Delta t$ in the sperical and elliptical trap"
    )
    plt.xlabel(r"$\Delta t$")
    plt.ylabel(r"Acceptence rate")
    plt.legend()

    tikzplotlib.save(
        "output/plots/acceptence_delta_t.tex",
        extra_axis_parameters=[
            "title style={align=center}",
            "xmajorticks=true",
            "ymajorticks=true",
            "mark options={mark size=2.5pt, line width=1.5pt}",
        ],
        strict=True,
    )

    plt.close()


def run_blocking():
    blocks = []
    for use_jastrow in ["interactive", "no_jastrow"]:
        print(use_jastrow)
        for particles in [10, 50, 100]:
            stderr = block(
                f"output/data/{use_jastrow}&particles={particles}_energy.dat"
            )
            print(f"{particles} particles, std err: {stderr}")


def color_tikz_plots():
    for filename in os.listdir("output/plots"):
        if filename.endswith(".tex"):
            with open(f"output/plots/{filename}", "r") as f:
                lines = f.readlines()
            with open(f"output/plots/{filename}", "w") as f:
                for line in lines:
                    if re_search := re.search(
                        r".*\[draw=white, fill=color(\d), mark=., only marks\].*",
                        line,
                    ):
                        line = line.replace("white", f"color{re_search.group(1)}")
                    f.write(line)


def write_latex_include(filename, caption, label):
    energy_plot = rf"""\begin{{figure}}[H]
    \centering
    \input{{{filename}}}
    \caption{{{caption}}}
    \label{{{label}}}
\end{{figure}}
"""
    print(energy_plot)


def create_latex_include_energy_simple():
    for num_particles in [1, 10, 100, 500]:
        for dimensions in [1, 2, 3]:
            write_latex_include(
                f"FYS4411 Project 1/Plot/v2/energy_particles={num_particles}_dimensions={dimensions}",
                f"Energy as a function of the variational parameter, $\\alpha$, for {num_particles} particle without interactions in {dimensions} dimension, for analytical and numerical double derivative.",
                f"fig:Energy_{num_particles}_nonint_{dimensions}D",
            )


def create_latex_include_energy_interactive():
    for num_particles in [10, 50, 100]:
        write_latex_include(
            f"FYS4411 Project 1/Plot/v2/energy_interactive_particles={num_particles}",
            f"Energy as a function of the variational parameter, $\\alpha$, for ${num_particles}$ particles with interactions in 3 dimensions.",
            f"fig:Energy_{num_particles}_int",
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "A program for doing different analysis for the boson simulation"
    )
    parser.add_argument(
        "-g", "--gradient", help="Plot the gradient descent", action="store_true"
    )
    parser.add_argument("-e", "--energy", help="Plot the energy", action="store_true")
    parser.add_argument("-d", "--delta-t", help="Plot the delta t", action="store_true")
    parser.add_argument(
        "-b",
        "--block",
        help="Do statistical analyis using blocking",
        action="store_true",
    )
    parser.add_argument(
        "-c",
        "--cpu-time",
        help="Plot the cpu time",
        action="store_true",
    )
    parser.add_argument(
        "-o",
        "--one-body",
        help="Plot the one body density",
        action="store_true",
    )
    parser.add_argument("-a", "--all", help="Create all the plots", action="store_true")
    args = parser.parse_args()

    if not any(
        [
            args.gradient,
            args.energy,
            args.delta_t,
            args.block,
            args.cpu_time,
            args.one_body,
            args.all,
        ]
    ):
        print("No argument parsed")
        parser.print_help()

    if args.gradient or args.all:
        plot_gradient_descent()

    if args.energy or args.all:
        plot_energy()
        create_latex_include_energy_simple()
        create_latex_include_energy_interactive()

    if args.delta_t or args.all:
        plot_delta_t()

    if args.block or args.all:
        run_blocking()

    if args.cpu_time or args.all:
        plot_cpu_time()

    if args.one_body or args.all:
        plot_one_body()

    color_tikz_plots()
