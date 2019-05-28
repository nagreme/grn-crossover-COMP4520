from config import Config
from grn import Grn

def main():
    pop = []
    for i in range(Config.pop_size):
        pop.append(Grn())

    simulate(pop)

#runs the regulatory simulation on the population
def simulate(pop):
    #we perform the simulation individually on each grn
    for i in range(Config.pop_size):
        print('Running simulation on GRN {}:'.format(i))
        grn = pop[i]
        #We run a fixed number of simulation timesteps, which can set in the config file
        for j in range(Config.sim_steps):
            print(' Step {}:'.format(j))
            #before we being the first timestep, insert the initial proteins in the Grn
            if j == 0:
                grn.push_initial_proteins()

            #The simulation consists of four steps:

            #first, we allow proteins to bind to genes
            bind_events = grn.bind()
            print_bind_events(bind_events)

            #second, active genes produce their product proteins
            produce_events = grn.produce()
            print_produce_events(produce_events)

            #third, the proteins diffuse (spread out) across the neighbouring genes
            grn.diffuse()

            #fourth, proteins "decay" (their concentrations are reduced by a small amount)
            grn.decay()

            print() #print a blank line to separate steps

    #evaluate the fitness of the individuals in the population
    eval_fitness(pop)

#update the fitness values of the networks in the population
def eval_fitness(pop):
    for i in range(Config.pop_size):
        grn = pop[i]
        #for now, we'll use the number of output proteins to determine fitness
        #(i.e. more output proteins = "better")
        #this is not very realistic, but it's good enough for testing purposes for now...
        #fitness value will be between 0 and 10
        #note: by convention, fitness is minimized in evolutionary algs (0 = best fitness, 10 = worst fitness)
        grn.fitness = max(10 - len(grn.output_proteins), 0)

def std_crossover(pop):
    # do we always do crossovers or do we have some probability of doing them? this parameter should be placed in the config file
    # this should probably be handled by the caller that's building the next generation population

    # roulette wheel selection of parents form pop
    select_parents(pop)

    # perform crossover, contruct two children
    # pick random indexm split parent gene arrays and initial proteins at that point
    # create two new grns and give them half and half


# use roulette wheel selection (slice size proportional to fitness)
def select_parents(pop):
    # this creates a new list, alternatively use pop.sort() to sort in place
    sorted_pop = sorted(pop, key=lambda x: x.fitness, reverse=True)
    # find total (sum) fitness
    # get random val (spin the wheel)
    # got through individuals from most to least fit, normalize value
    # see if the wheel landed in their chunk

    return 0


def print_bind_events(events):
    for (protein, gene) in events:
        print('  Protein {} bound to Gene {}'.format(protein.seq.to01(), gene.index))

    if not events:
        print('  No bind events')

def print_produce_events(events):
    for (gene, protein) in events:
        print('  Gene {} produced Protein {}'.format(gene.index, protein.seq.to01()))

    if not events:
        print('  No produce events')

main()
