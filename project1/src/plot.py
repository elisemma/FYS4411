import csv
from matplotlib import pyplot as plt


def plot_energy():
    particles1 = []
    particles2 = []
    particles3 = []
    with open("output/data/energy_values.tsv") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        for row in reader:
            print(row[2])
            if row[2] == "1":
                particles1.append(float(row[3]))
            elif row[2] == "4":
                particles2.append(float(row[3]))
            elif row[2] == "6":
                particles3.append(float(row[3]))

    plt.plot(particles1[1:], label="1")
    plt.plot(particles2[1:], label="4")
    plt.plot(particles3[1:], label="6")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    plot_energy()
