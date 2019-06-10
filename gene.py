from config import Config
from utils import Utils
from protein import *
import numpy as np

class Gene():
    def __init__(self, index):
        #the position of this gene within the enclosing Grn's gene array
        self.index = index

        #a binary sequence that proteins can bind to in order to activate this gene
        self.binding_seq = Utils.rand_bitvec(Config.num_protein_bits)

        #a binary sequence that will be produced (in the form of a protein) when this gene is activated
        self.product_seq = Utils.rand_bitvec(Config.num_protein_bits)

        #the currently bound protein (if any)
        self.bound_protein = None

        #the protein currently being produced by this gene (if any)
        #note: this points to a protein object, not a sequence
        self.product_protein = None

        #the amount of protein that is produced by this gene (if it's active) at each simulation timestep
        #this is a floating point value between 0 and 1
        self.production_rate = np.random.random()

        #the concentration of protein that must be present above this gene in order for the protein to be
        #considered for binding. This is a floating point value between 0 and 1
        self.bind_threshold = np.random.random()

    #binds the given protein to this gene
    #if this gene is not already active, it begins to produce its product protein
    def bind(self, protein):
        self.bound_protein = protein
        if self.product_protein is None:
            #note: we initialize the protein with zero concentrations.
            #The concentration will be incremented when this gene produce the product (see grn.produce(), which is called after the binding has completed)
            self.product_protein = Protein(self.index, self.product_seq, zero_concs=True)

    #decouples a bound protein from this gene
    #if this gene was producing product protein, it stops
    def clear_binding(self):
        self.bound_protein = None
        self.product_protein = None

    #checks to see if the given protein would meet the criteria for binding to this gene
    #returns true or false
    def can_bind(self, protein):
        #find the minimum number of bits that must match in order for binding to occur
        min_affinity = Config.num_protein_bits - Config.binding_seq_play

        #Binding Criteria 1: check to see if the protein sequence matches this gene's binding sequence closely enough
        above_affinity = Utils.same_count(self.binding_seq, protein.seq) > min_affinity

        #Binding Criteria 2: check to see if the protein has a high enough concentration to bind to this gene (using this gene's binding threshold)
        above_threshold = protein.concs[self.index] >= self.bind_threshold

        #in order for the protein to bind, it must meet both of the above criteria
        return above_affinity and above_threshold


    def __repr__(self):
        return "{}.{}(index={}, binding_seq={}, product_seq={}, bound_protein={}, production_rate={}, bind_threshold={})".format(self.__module__, type(self).__name__, self.index, self.binding_seq.to01(), self.product_seq.to01(), self.bound_protein, self.production_rate, self.bind_threshold)

    def __str__(self):
        return "Gene: {}||{} (bound protein: {})".format(self.binding_seq.to01(), self.product_seq.to01(), self.bound_protein)
        # return "G:{}".format(self.binding_seq.to01())
