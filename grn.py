from config import Config
from bitarray import bitarray
from utils import Utils
from protein import *
from gene import Gene
from collections import OrderedDict
import numpy as np
import copy

#This class represents a single gene regulatory network
class Grn():
    def __init__(self):
        #current fitness of the network
        #right now, this is a value between 0 and 10
        #note: by convention, fitness is minimized in evolutionary algs (0 = best fitness, 10 = worst fitness)
        #I'm initializing to the illegal value -1 here to indicate that the fitness hasn't been evaluated yet
        #it will be updated after the simulation is finished (see sim.eval_fitness())
        self.fitness = -1
        
        #these proteins bind to other genes (and do not affect program output)
        self.internal_proteins = {}
        
        #these proteins affect program output (and do not bind to other genes)
        self.output_proteins = {}
        
        #these proteins are injected into the GRN at the start of the regulatory simulation
        #I'm using an ordered dictionary here - it works just like a regular dictionary, except you can iterate
        #through it in a fixed order (like an array). This should make it a little easier to split the container in half
        #when you're doing crossover.
        self.initial_proteins = OrderedDict()
        
        #randomly initialize the initial proteins
        #note: we want to make sure we have no duplicate sequences here
        #purpose of the "tries" var below is just in case we don't have enough protein bits to generate unique combinations within a reasonable timeframe...
        tries = 0
        while len(self.initial_proteins) < Config.num_initial_proteins and tries < 100 * Config.num_initial_proteins: 
            seq = Utils.rand_bitvec(Config.num_protein_bits - 1)
            seq.insert(0, 0) #make this an internal protein so that it will bind

            if seq.to01() not in self.initial_proteins:
                protein = Protein(-1, seq=seq) #there is no src gene that produced this protein. Use a -1 to indicate this
                self.initial_proteins[protein.seq.to01()] = protein

            else:
                tries += 1

        if len(self.initial_proteins) < Config.num_initial_proteins:
            print('Unable to generate unique random initial protein sequences.')
            print('Try increasing the Config.num_protein_bits, or decreasing the number of initial proteins.')
            exit(1)

        #Create the genes and initialze them randomly
        self.genes = []
        for i in range(Config.num_genes):
            self.genes.append(Gene(i))

    #inserts the initial proteins into the network
    #this must be done before the simulation starts
    def push_initial_proteins(self):
        for key in self.initial_proteins:
            protein = self.initial_proteins[key]
            #insert a copy, so that if the protein decays so much that it gets removed, we don't lose it for good
            clone = copy.deepcopy(protein)
            self.internal_proteins[clone.seq.to01()] = clone

    #binds internal proteins to genes
    #returns a list of bind events, where each element is a 2-tuple of the form (protein that bound, gene that it bound to)
    def bind(self):
        bind_events = []
        
        #go through every gene in the network and look at the proteins that exist above it
        for pos in range(Config.num_genes):
            gene = self.genes[pos]
            concs_at_pos = [] #list of tuples of the form (protein, concentration value) for all proteins above the current position
            conc_sum_at_pos = 0.0
            #go through all internal proteins and record the concentration that is present over the current gene
            for seq in self.internal_proteins:
                protein = self.internal_proteins[seq]
                #we only need to look at proteins that have a sequence that is a close enough match to the binding site sequence
                if gene.can_bind(protein):
                    concs_at_pos.append( (protein, protein.concs[pos]) )
                    conc_sum_at_pos += protein.concs[pos]

            #if there were some proteins above this gene that are capable of binding, pick one
            if len(concs_at_pos) > 0:
                #To select a protein:
                #Create a "roulette wheel" (i.e. a probability distribution), and assign each protein a chunk of it.
                #The chunk is proportional in size to the amount of protein present over the current gene.
                #Next, "spin" the wheel to select a protein (i.e. generate a random number and see which chunk's
                #range it falls into)
                
                index = 0 #index of current protein from concs_at_pos list
                #the first protein's chunk will go from 0.0 to the following value (note: we're normalizing on the fly here so that the sum of all of the chunks will be 1 - that's what the division's for)
                cutoff = concs_at_pos[index][1] / conc_sum_at_pos
                #"spin" our wheel
                r = np.random.random()
                #now, go through the chunks, computing their start and end points,
                #until you find the chunk who's range contains the random value r
                #eg. Suppose our chunks are:
                # 0.0 - 0.25, 0.26 - 0.5, 0.51 - 0.75, 0.76 - 1.0
                #and suppose our random value r = 0.4
                #Then, on the first iteration, cutoff will be 0.25 (end of first chunk). Loop will continue.
                #On second iteration, cutoff will be 0.5 (end of second chunk). Loop will stop. Index will be 1.
                #So we will select chunk 1 (i.e. protein 1 from the concs_at_pos list).
                while r > cutoff:
                    index += 1
                    cutoff += concs_at_pos[index][1] / conc_sum_at_pos

                #grab the selected protein and bind it to the gene
                selected = concs_at_pos[index][0]
                gene.bind(selected)
                bind_events.append( (selected, gene) )
                
                self._add_product_protein(gene)

            #if there were no viable proteins above the current gene, clear any existing binding from previous iteration
            else:
                gene.clear_binding()

        return bind_events

    #goes through the genes and checks to see if any of them have activated
    def _add_product_protein(self, gene):
        if gene.product_protein is not None:
            #check if protein already exists in this grn - if so, just add the current gene as a src
            key = gene.product_seq.to01()
            if key in self.internal_proteins:
                self.internal_proteins[key].add_src(gene.index)
                gene.product_protein = self.internal_proteins[key]

            elif key in self.output_proteins:
                self.output_proteins[key].add_src(gene.index)
                gene.product_protein = self.output_proteins[key]

            #otherwise protein has not yet been inserted into this Grn, so do that
            else:
                if gene.product_protein.type == ProteinTypes.INTERNAL:
                    self.internal_proteins[gene.product_seq.to01()] = gene.product_protein

                else:
                    self.output_proteins[gene.product_seq.to01()] = gene.product_protein
            

    #causes active genes to product their output protein
    #returns a list of production events, where each element is a 2-tuple of the form (gene that produced, protein that was produced)
    def produce(self):
        produce_events = []
        
        for pos in range(Config.num_genes):
            gene = self.genes[pos]
            if gene.product_protein is not None:
                protein = gene.product_protein
                #increment the concentration in the position directly above this gene by the gene's production rate
                #We need to clamp to make sure the conc stays between 0 and 1.
                #Note: the protein will diffuse into the neighbouring positions in the following timestamps (see diffuse())
                protein.concs[pos] = Utils.clamp(protein.concs[pos] + gene.production_rate, 0.0, 1.0)
                produce_events.append( (gene, protein) )

        return produce_events

    #causes the proteins to "spread outward" from the src gene that produced them
    def diffuse(self):
        self._diffuse(self.internal_proteins)
        self._diffuse(self.output_proteins)

    #This method implements a discritized differential equation to do the actual diffusion
    #This approach was adapted from the "Forward Euler scheme" described here:
    #http://hplgit.github.io/num-methods-for-PDEs/doc/pub/diffu/sphinx/._main_diffu001.html
    def _diffuse(self, proteins):
        F = 0.3 #a "physical" constant used to control the behaviour of the liquid that's diffusing
        #these parameters control the left/right balance of the diffusion. If one is higher, more protein diffuses
        #in that direction (like a wave).
        right = 1.0
        left = 1.0

        #see the diffusion equation paper for more details
        for key in proteins:
            protein = proteins[key]
            prev_concs = protein.concs

            protein.concs[0] = prev_concs[0] + F * (-2.0 * prev_concs[0] + left * prev_concs[1])

            for j in range(1, Config.num_genes - 1):
                new_conc = prev_concs[j] + F * (right * prev_concs[j - 1] - 2.0 * prev_concs[j] + left * prev_concs[j + 1])
                protein.concs[j] = min(new_conc, 1.0)

            protein.concs[Config.num_genes - 1] = prev_concs[Config.num_genes - 1] + F * (right * prev_concs[Config.num_genes - 1] - 2.0 * prev_concs[Config.num_genes - 1])

    #causes the concentration levels of all proteins to decrease slightly ("decay")
    #Using decay means that a protein can only occupy a binding site as long as it it is continually reinforced.
    def decay(self):
        self._decay(self.internal_proteins)
        self._decay(self.output_proteins)

    #do the actual decay
    #After decay, this method also checks to see if a protein has dropped below the minimum allowable concentration level.
    #if so, it is removed from the Grn.
    def _decay(self, proteins):
        keys_to_remove = []
        
        for key in proteins:
            protein = proteins[key]
            above_threshold = False
            for pos in range(len(protein.concs)):
                #compute how much the protein should decay (based on the Config.decay_factor, then subtract that amount from the current position.
                #note: we need to make sure we don't accidentally cause the conc to dip below zero - that's what the max() is for
                protein.concs[pos] = max(0.0, protein.concs[pos] - protein.concs[pos] * Config.decay_factor)
                #check if it's still above the threshold
                #note: the entire protein (all positions) must be below the threshold in order for the protein to be removed from the Grn - hence the "or" operation here
                above_threshold = above_threshold or protein.concs[pos] >= Config.min_protein_conc

            if not above_threshold:
                #look through the gene(s) that produced this protein
                for src in protein.srcs:
                    #if a gene is still producing the protein, stop production
                    if self.genes[src].product_protein == protein:
                        self.genes[src].product_protein = None

                #add to the remove list (note: we can't remove it right here because changing the dictionary size will mess up the outer loop's iteration)
                keys_to_remove.append(key)

        #remove the decayed proteins from this grn
        for key in keys_to_remove:
            del proteins[key]
