#This class holds all of the adjustable parameters that the simulation uses.
#These values remain fixed throughout the simulation.
class Config():
    #self-explanetory parameters
    pop_size = 3
    num_genes = 8
    num_initial_proteins = 4
    sim_steps = 2
    worst_fitness = 10 # inclusive

    #number of bits in a protein sequence, binding sequence, or output sequence
    num_protein_bits = 4

    #when binding, a protein sequence and binding sequence may differ by a maximum of this number of bits
    binding_seq_play = 2

    #proteins are removed from the simulation if their concentration level dips below this threshold (across all genes)
    min_protein_conc = 0.01

    #proteins decay at the end of each timestep. This value (between 0 and 1) indicates the proportion protein that is lost in this decay.
    #(i.e. new_protein_conc = curr_conc - curr_conc * decay_factor)
    decay_factor = 0.2
