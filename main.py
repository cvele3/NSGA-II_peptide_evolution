from PeptideEvolutionNSGAII import NSGA_II
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

def visualize_pareto_fronts(pareto_fronts):

    ax = plt.figure().gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    plt.title("Prikaz pareto fronti")
    plt.xlabel("Pareto fronta")
    plt.ylabel("Vjerojatnost postojanja AMP svojstva (ff_amp_probability)")

    colors = ["#" + ''.join([np.random.choice(list('0123456789ABCDEF')) for _ in range(6)])
              for _ in range(len(pareto_fronts))]

    for front_index, front in enumerate(pareto_fronts):
        for peptide in front:
            _, _, ff_amp_probability, ff_length = peptide
            plt.scatter(front_index, ff_amp_probability, c=colors[front_index])

    plt.show()


GA = NSGA_II(
    lowerRange=8,
    upperRange=15,
    population_size=100,
    offspring_size=20,
    num_generations=30,
    num_solutions_tournament=2,
    mutation_probability=0.1,
    penalty_function_reducer=0.5
)

pareto_fronts = GA.calculate()

visualize_pareto_fronts(pareto_fronts)

for solution in pareto_fronts[0]:
    print(solution)