from config import Config
from grn import Grn

def main():
    pop = []
    for i in range(Config.pop_size):
        pop.append(Grn())

    simulate(pop)

def simulate(pop):
    for i in range(Config.pop_size):
        grn = pop[i]
        for j in range(Config.sim_steps):
            if j == 0:
                grn.push_initial_proteins()

            grn.bind()
            grn.produce()
            grn.diffuse()
            grn.decay()
    
main()
