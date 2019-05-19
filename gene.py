from config import Config
from utils import Utils
from protein import *
import numpy as np

class Gene():
    def __init__(self, index):
        self.index = index
        self.binding_seq = Utils.rand_bitvec(Config.num_protein_bits)
        self.product_seq = Utils.rand_bitvec(Config.num_protein_bits)
        self.bound_protein = None
        self.product_protein = None
        self.production_rate = np.random.random()
        self.bind_threshold = np.random.random()

    def bind(self, protein):
        self.bound_protein = protein
        if self.product_protein is None:
            self.product_protein = Protein(self.index, self.product_seq)

    def clear_binding(self):
        self.bound_protein = None
        self.product_protein = None

    def can_bind(self, protein):
        min_affinity = Config.num_protein_bits - Config.binding_seq_play
        above_affinity = Utils.same_count(self.binding_seq, protein.seq) > min_affinity
        above_threshold = protein.concs[self.index] >= self.bind_threshold

        return above_affinity and above_threshold
