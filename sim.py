from config import Config
from grn import Grn
from numpy import random
from collections import OrderedDict
from utils import Utils
from protein import *

def main():
    pop = []
    for i in range(Config.pop_size):
        pop.append(Grn())

    simulate(pop)

    # print("Population")
    # print(pop[0])
    # print()
    # print(pop[1])
    # print()
    # print(pop[2])
    # print()
    #
    # print("Crossover")
    # print()
    # c1,c2 = std_crossover(pop)
    #
    # print("Children")
    # print(c1)
    # print(c1.initial_proteins)
    # print()
    # print(c2)
    # print(c2.initial_proteins)
    # print()



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
    p1, p2 = select_parents(pop)
    # p1 and p2 are Grn objects

    # print("Selected parents:")
    # print(p1)
    # print(p1.initial_proteins)
    # print()
    # print(p2)
    # print(p2.initial_proteins)
    # print()

    # perform crossover, contruct two children

    # create two new grns and give them half and half
    c1 = Grn()
    c2 = Grn()
    # add stuff to children (genes and initial proteins)

    # we need to reset the initial proteins since they're randomly generated in the constructor
    # we could pass a flag to the constructor but I'd rather not change it for now
    c1.initial_proteins = OrderedDict()
    c2.initial_proteins = OrderedDict()

    # pick random indexm split parent gene arrays and initial proteins at that point
    split_index = random.randint(0, Config.num_genes)

    # Note: split index of zero means no crossing over occurs

    # print("gene split index:", split_index)
    # print()

    #TODO: Currently genes in the children keep their bound proteins:
    # I could null these if the parent doent' need them anymore (depends on the simulation process)
    # Or I could create new gene objects and manually set the binding and product sequences (and bind thresh and production rate)

    c1.genes = p1.genes[0:split_index] + p2.genes[split_index:]
    c2.genes = p2.genes[0:split_index] + p1.genes[split_index:]

    split_index = random.randint(0, Config.num_initial_proteins)

    # print("initial prot split index:", split_index)
    # print()

    # pull the keys out of the ordered dict
    # can't use .keys() because it returns an odict_keys object that doesn't support indexing
    p1_prot_list = list(p1.initial_proteins)
    p2_prot_list = list(p2.initial_proteins)

    for i in range(0, Config.num_initial_proteins):
        # pull out next initial protein
        p1_prot = p1.initial_proteins[p1_prot_list[i]]
        p2_prot = p2.initial_proteins[p2_prot_list[i]]


        if i < split_index:
            c1.initial_proteins[p1_prot_list[i]] = p1_prot
            c2.initial_proteins[p2_prot_list[i]] = p2_prot
        else:
            c1.initial_proteins[p2_prot_list[i]] = p2_prot
            c2.initial_proteins[p1_prot_list[i]] = p1_prot



    # check if both children have enough initial_proteins
    check_init_prots(c1)
    check_init_prots(c2)

    return c1,c2


# use roulette wheel selection (slice size proportional to fitness) (implemented with replacement)
# returns the Grn objects that are the selected parents
def select_parents(pop):
    # this creates a new list, alternatively use pop.sort() to sort in place
    sorted_pop = sorted(pop, key=lambda x: x.fitness, reverse=False)

    # find total (sum) fitness
    fitness_sum = 0

    for individual in sorted_pop:
        # Recall best fitness = 0, worst fitness = 10
        fitness_sum += (Config.worst_fitness - individual.fitness)

    # get random vals for each parent (spin the wheel)
    p1 = random.random()
    p2 = random.random()

    # keep track with flags
    p1_found = False
    p2_found = False

    p1_index = -1
    p2_index = -1

    # go through individuals from most to least fit grn, normalize value (this value is the size of their wheel chunk)
    index = -1 # index of current grn
    cutoff = 0 # where we're at in the wheel

    # print("selecting parents")
    # print(p1,p2)
    # print(0.0)

    # keep going until we find both parents
    while not (p1_found and p2_found):
        index += 1
        cutoff += (Config.worst_fitness - sorted_pop[index].fitness) / fitness_sum
        # print(cutoff)

        if  not p1_found and p1 < cutoff:
            p1_index = index
            p1_found = True

        if  not p2_found and p2 < cutoff:
            p2_index = index
            p2_found = True

    # print()
    return sorted_pop[p1_index], sorted_pop[p2_index]


def check_init_prots(child):
    # Code copied from grn constructor but without the tries because we really need this to work
    while len(child.initial_proteins) < Config.num_initial_proteins:
        print("smaller")
        print(child.initial_proteins)
        seq = Utils.rand_bitvec(Config.num_protein_bits - 1)
        seq.insert(0, 0) # make this an internal protein so that it will bind

        if seq.to01() not in child.initial_proteins:
            protein = Protein(-1, seq=seq) # there is no src gene that produced this protein. Use a -1 to indicate this
            child.initial_proteins[protein.seq.to01()] = protein


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
