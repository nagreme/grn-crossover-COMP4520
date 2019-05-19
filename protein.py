from config import Config
from utils import Utils
from enum import Enum
import numpy as np

class ProteinTypes(Enum):
    INTERNAL = 0
    OUTPUT = 1

class Protein():
    def __init__(self, src, seq=None, zero_concs=False):
        if seq is None:
            self.seq = Utils.rand_bitvec(Config.num_protein_bits)
        else:
            self.seq = seq

        if zero_concs:
            self.concs = np.zeros(Config.num_genes)
        else:
            self.concs = np.random.random(Config.num_genes)
            
        self.type = ProteinTypes(self.seq[0]) #msb determines type
        self.srcs = set( (src,) ) #use a set to prevent duplicates and provide lookup in constant time

    def add_src(self, src):
        self.srcs.add(src)
