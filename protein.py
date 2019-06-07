from config import Config
from utils import Utils
from enum import Enum
import numpy as np

#there are 2 types of proteins in the simulation:
class ProteinTypes(Enum):
    #these ones can only bind other genes (they cannot affect the generated program code)
    INTERNAL = 0

    #these ones can only be used to affect the generated program code (they cannot bind)
    OUTPUT = 1

class Protein():
    #Parameters:
    #src [int] is the index of the gene that produced this protein, or -1 if this is an initial protein
    #seq [bitarray] is the binary string to use for this protein.
    #           If this argument is not passed, the protein sequence is randomly initialized
    #zero_concs [bool] may be set to True to create a protein with a concentration of zero above each gene.
    #               Otherwise, the concentration levels are randomly initialized
    def __init__(self, src, seq=None, zero_concs=False):
        if seq is None:
            self.seq = Utils.rand_bitvec(Config.num_protein_bits)
        else:
            self.seq = seq

        if zero_concs:
            self.concs = np.zeros(Config.num_genes)
        else:
            self.concs = np.random.random(Config.num_genes)

        self.type = ProteinTypes(self.seq[0]) #the top bit determines the protein's type

        #note: a protein may have more than one source gene (i.e. two genes may produce the same protein)
        #use a set to prevent duplicates and provide lookup in constant time
        self.srcs = set()
        self.srcs.add(src)

    def add_src(self, src):
        self.srcs.add(src)

    def __repr__(self):
        return "{}.{}(seq={}, srcs={})".format(self.__module__, type(self).__name__, self.seq.to01(), self.srcs)

    def __str__(self):
        return "P:{}".format(self.seq.to01())
