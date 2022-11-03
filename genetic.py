from random import choices, randint, randrange, random
from typing import List, Optional, Callable, Tuple
from PileLoads import PileGroup
from Optym_GA import Optimize
from functools import partial
import time

start_time = time.time()

Genome = List[int]
Population = List[Genome]
PopulateFunc = Callable[[], Population]
FitnessFunc = Callable[[Genome], int]
SelectionFunc = Callable[[Population, FitnessFunc], Tuple[Genome, Genome]]
CrossoverFunc = Callable[[Genome, Genome], Tuple[Genome, Genome]]
MutationFunc = Callable[[Genome], Genome]
PrinterFunc = Callable[[Population, int, FitnessFunc], None]


def generate_genome(length: int) -> Genome:
    genom = []
    for i in range(int(length//2)):
        genom.append(choices([0, 1])[0])
        genom.append(choices([0, 1,2,3])[0])
    return genom


def generate_population(size: int, genome_length: int) -> Population:
    return [generate_genome(genome_length) for _ in range(size)]





def single_point_crossover(a: Genome, b: Genome) -> Tuple[Genome, Genome]:
    if len(a) != len(b):
        raise ValueError("Genomes a and b must be of same length")

    length = len(a)
    if length < 2:
        return a, b

    p = randint(1, length - 1)
    return a[0:p] + b[p:], b[0:p] + a[p:]


def mutation(genome: Genome, num: int = 1, probability: float = 0.5) -> Genome:
    for _ in range(num):
        index = randrange(len(genome))
        if random() > probability:
            genome[index] = genome[index]
        else:
            if index%2==0:
                genome[index] = abs(genome[index] - 1)
            else:
                genome[index] = choices([0, 1, 2, 3])[0]
    return genome


def population_fitness(population: Population, fitness_func: FitnessFunc) -> int:
    return sum([fitness_func(genome) for genome in population])


def selection_pair(population: Population, fitness_func: FitnessFunc) -> Population:
    return choices(
        population=population,
        weights=[1/fitness_func(gene) for gene in population],
        k=2
    )


def sort_population(population: Population, fitness_func: FitnessFunc) -> Population:
    return sorted(population, key=fitness_func, reverse=False)


def genome_to_string(genome: Genome) -> str:
    return "".join(map(str, genome))


def print_stats(population: Population, generation_id: int, fitness_func: FitnessFunc):
    print("GENERATION %02d" % generation_id)
    print("=============")
    print("Population: [%s]" % ", ".join([genome_to_string(gene) for gene in population]))
    print("Avg. Fitness: %f" % (population_fitness(population, fitness_func) / len(population)))
    sorted_population = sort_population(population, fitness_func)
    print(
        "Best: %s (%f)" % (genome_to_string(sorted_population[0]), fitness_func(sorted_population[0])))
    print("Worst: %s (%f)" % (genome_to_string(sorted_population[-1]),
                              fitness_func(sorted_population[-1])))
    print("")

    return sorted_population[0]


def fitness(genom):
    global Piles_grp, Result
    Piles_grp.Generate_piles(genom)

    Results.import_piles(Piles_grp.piles)
    Results.piles_matrix()
    Results.forces()
    normal_forces = Results.single_normal_forces()
    deformation = Results.single_deformation()



    sum_Pmax =0
    sum_Pmin = 0
    n_min=0
    n_max=0
    for force in normal_forces:
        if force >0:
            if force <Piles_grp.P_max:
                sum_Pmax += Piles_grp.P_max-force
            else:
                sum_Pmax += 2*abs(Piles_grp.P_max - force)
            n_max +=1
        else:
            if force>Piles_grp.P_min:
                sum_Pmin += abs(Piles_grp.P_min-force)
            else:
                sum_Pmin += 2*abs(Piles_grp.P_min - force)
            n_min +=1

    #print(f'{sum_Pmax}....{sum_Pmin}...{(sum_Pmax/P_max+sum_Pmin/abs(P_min)) /(n_max+n_min)}')
    value =  (sum_Pmax/Piles_grp.P_max+sum_Pmin/abs(Piles_grp.P_min)) /(n_max+n_min)+max((abs(deformation)-Piles_grp.def_lim)/Piles_grp.def_lim,0)
    return value




def run_evolution(
        populate_func: PopulateFunc,
        fitness_func: FitnessFunc,
        selection_func: SelectionFunc = selection_pair,
        crossover_func: CrossoverFunc = single_point_crossover,
        mutation_func: MutationFunc = mutation,
        generation_limit: int = 100,
        printer: Optional[PrinterFunc] = None) \
        -> Tuple[Population, int]:
    population = populate_func()

    for i in range(generation_limit):
        population = sorted(population, key=lambda genome: fitness_func(genome), reverse=False)

        if printer is not None:
            printer(population, i, fitness_func)

        if fitness_func(population[0]) <= 0:
            break

        next_generation = population[0:2]

        for j in range(max(int((len(population) / 2 - 1)*.7),10)):
            parents = selection_func(population, fitness_func)
            offspring_a, offspring_b = crossover_func(parents[0], parents[1])
            offspring_a = mutation_func(offspring_a)
            offspring_b = mutation_func(offspring_b)
            next_generation += [offspring_a, offspring_b]

#        for j in range(int(len(population) / 2) - 1):
#            parents = selection_func(population, fitness_func)
#            offspring_a, offspring_b = crossover_func(parents[0], parents[1])
#            offspring_a = mutation_func(offspring_a)
#            offspring_b = mutation_func(offspring_b)
#            next_generation += [offspring_a, offspring_b]

        population = next_generation

    return population, i




Piles_grp = Optimize(P_max=2000,P_min=-1500,def_lim=32)
Piles_grp.Generate_allowed_piles(B_min=-6, B_max=10, L_min=-5, L_max=5, dx=2, dy=2, d_edge=0.6,N=4,L=30)
Results = PileGroup()
Results.pile_parameter(Ep=33 * 10 ** 9, Gp=12 * 10 ** 9, Ap=109378 * 10 ** -6, Jp=58791 * 10 ** (4 - 4 * 3), m=0,
                       soil_type="cohesive", kd=4615 * 1000)
Results.single_load('Loads.csv')


genom = generate_genome(2*len(Piles_grp.all_piles))
Piles_grp.Generate_piles(genom)



Results.import_piles(Piles_grp.all_piles)

Results.piles_matrix()
Results.forces()
print(genom)
print(fitness(genom))


population, generations = run_evolution(
    populate_func=partial(generate_population, size=50, genome_length=2*len(Piles_grp.all_piles)),
    fitness_func=partial(fitness),
    generation_limit=60
)
print(population[0])
print(f'fitness {fitness(population[0])}')
Piles_grp.Generate_piles(population[0])

Results.import_piles(Piles_grp.piles)
Results.piles.to_csv('final_piles.csv')
print( Results.single_normal_forces())
print(Results.single_deformation())
print("--- %s seconds ---" % (time.time() - start_time))