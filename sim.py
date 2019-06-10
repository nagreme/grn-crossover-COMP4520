from config import Config
from grn import Grn
import networkx as nx
import matplotlib.pyplot as plt


IPROT_COLOUR = '#678fe0'
PROT_COLOUR = '#67cce0'
OPROT_COLOUR = '#67e0b8'


def main():
    pop = []
    for i in range(Config.pop_size):
        pop.append(Grn())

    graphs = []

    simulate(pop, graphs)




#runs the regulatory simulation on the population
def simulate(pop, graphs):
    #we perform the simulation individually on each grn
    for i in range(Config.pop_size):
        print('Running simulation on GRN {}:'.format(i))
        grn = pop[i]
        DG = nx.DiGraph()
        graphs.append(DG)

        #We run a fixed number of simulation timesteps, which can set in the config file
        for j in range(Config.sim_steps):
            print(' Step {}:'.format(j))
            #before we being the first timestep, insert the initial proteins in the Grn
            if j == 0:
                grn.push_initial_proteins()

                # add the initial proteins to the graph so we can colour them differently
                for k in grn.internal_proteins.keys():
                    DG.add_node('P:'+ k)


            #The simulation consists of four steps:

            #first, we allow proteins to bind to genes
            bind_events = grn.bind()
            print_bind_events(bind_events, DG)

            #second, active genes produce their product proteins
            produce_events = grn.produce()
            print_produce_events(produce_events, DG)

            #third, the proteins diffuse (spread out) across the neighbouring genes
            grn.diffuse()

            #fourth, proteins "decay" (their concentrations are reduced by a small amount)
            grn.decay()

            print() #print a blank line to separate steps

        # set up colours for the different types of protein
        # recall all the initial proteins are at the beginning of the node list
        # and the output proteins' sequences start with '1'
        colors = [IPROT_COLOUR for i in range(Config.num_initial_proteins)] + [PROT_COLOUR for i in range(len(DG.nodes)-Config.num_initial_proteins)]

        node_index = 0
        for node in DG.nodes():
            if 'P:1' in node:
                colors[node_index] = OPROT_COLOUR
            node_index += 1

        nx.draw(DG, with_labels=True, node_color=colors, arrowsize=15, alpha=0.7, edge_color='#555555')
        plt.show()
        plt.clf()

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

def print_bind_events(events, DG):
    for (protein, gene) in events:
        print('  Protein {} bound to Gene {} ({})'.format(protein.seq.to01(), gene.index, gene.binding_seq.to01()))
        DG.add_edge('P:'+protein.seq.to01(), 'G{}:{}'.format(gene.index, gene.binding_seq.to01()))

    if not events:
        print('  No bind events')

def print_produce_events(events, DG):
    for (gene, protein) in events:
        print('  Gene {} ({}) produced Protein {}'.format(gene.index, gene.binding_seq.to01(), protein.seq.to01()))
        DG.add_edge('G{}:{}'.format(gene.index, gene.binding_seq.to01()), 'P:'+protein.seq.to01())

    if not events:
        print('  No produce events')

main()
