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
        grn = pop[i]
        #We run a fixed number of simulation timesteps, which can set in the config file
        for j in range(Config.sim_steps):
            #before we being the first timestep, insert the initial proteins in the Grn
            if j == 0:
                grn.push_initial_proteins()

            #The simulation consists of four steps:
            
            #first, we allow proteins to bind to genes
            grn.bind()

            #second, active genes produce their product proteins
            grn.produce()

            #third, the proteins diffuse (spread out) across the neighbouring genes
            grn.diffuse()

            #fourth, proteins "decay" (their concentrations are reduced by a small amount)
            grn.decay()
    
main()
