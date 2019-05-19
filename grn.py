from config import Config
from bitarray import bitarray
from utils import Utils
from protein import *
from gene import Gene
import numpy as np
import copy

class Grn():
    def __init__(self):
        self.internal_proteins = {}
        self.output_proteins = {}
        self.initial_proteins = {}
        for i in range(Config.num_initial_proteins):
            seq = Utils.rand_bitvec(Config.num_protein_bits - 1)
            seq.insert(0, 0) #make this an internal protein
            protein = Protein(-1, seq=seq)
            #note: if we have duplicates, hashing them like this will eliminate them
            #the grn will then have one fewer initial protein...
            self.initial_proteins[protein.seq.to01()] = protein
        
        self.genes = []
        for i in range(Config.num_genes):
            self.genes.append(Gene(i))

    def push_initial_proteins(self):
        for key in self.initial_proteins:
            protein = self.initial_proteins[key]
            clone = copy.deepcopy(protein)
            self.internal_proteins[clone.seq.to01()] = clone
            
    def bind(self):
        for pos in range(Config.num_genes):
            gene = self.genes[pos]
            already_bound = set()
            concs_at_pos = []
            conc_sum_at_pos = 0.0
            for seq in self.internal_proteins:
                protein = self.internal_proteins[seq]
                if protein not in already_bound and gene.can_bind(protein):
                    concs_at_pos.append( (protein, protein.concs[pos]) )
                    conc_sum_at_pos += protein.concs[pos]

            if len(concs_at_pos) > 0:
                r = np.random.random()

                index = 0
                cutoff = concs_at_pos[index][1] / conc_sum_at_pos
                while r > cutoff:
                    index += 1
                    cutoff += concs_at_pos[index][1] / conc_sum_at_pos

                selected = concs_at_pos[index][0]
                gene.bind(selected)
                already_bound.add(selected)

            else:
                gene.clear_binding()

    def regulate(self):
        for pos in range(Config.num_genes):
            gene = self.genes[pos]
            if gene.product_protein is None:
                #check if protein already exists - if so, just add the current gene as a src
                key = gene.product_seq.to01()
                if key in self.internal_proteins:
                    self.internal_proteins[key].add_src(pos)
                    gene.product_protein = self.internal_proteins[key]

                elif key in self.output_proteins:
                    self.output_proteins[key].add_src(pos)
                    gene.product_protein = self.output_proteins[key]

                #otherwise protein does not exist yet - create it
                else:
                    gene.product_protein = Protein(pos, seq=gene.product_seq, zero_concs=True)
                    if gene.product_protein.type == ProteinTypes.INTERNAL:
                        self.internal_proteins[gene.product_seq.to01()] = gene.product_protein

                    else:
                        self.output_proteins[gene.product_seq.to01()] = gene.product_protein
            

    def produce(self):
        for pos in range(Config.num_genes):
            gene = self.genes[pos]
            if gene.product_protein is not None:
                protein = gene.product_protein
                protein.concs[pos] = Utils.clamp(protein.concs[pos] + gene.production_rate, 0.0, 1.0)

    def diffuse(self):
        self._diffuse(self.internal_proteins)
        self._diffuse(self.output_proteins)

    def _diffuse(self, proteins):
        F = 0.3
        right = 1.0
        left = 1.0

        to_remove = []
        for key in proteins:
            protein = proteins[key]
            prev_concs = protein.concs

            protein.concs[0] = prev_concs[0] + F * (-2.0 * prev_concs[0] + left * prev_concs[1])

            for j in range(1, Config.num_genes - 1):
                new_conc = prev_concs[j] + F * (right * prev_concs[j - 1] - 2.0 * prev_concs[j] + left * prev_concs[j + 1])
                protein.concs[j] = min(new_conc, 1.0)

            protein.concs[Config.num_genes - 1] = prev_concs[Config.num_genes - 1] + F * (right * prev_concs[Config.num_genes - 1] - 2.0 * prev_concs[Config.num_genes - 1])

    def decay(self):
        self._decay(self.internal_proteins)
        self._decay(self.output_proteins)
                
    def _decay(self, proteins):
        keys_to_remove = []
        
        for key in proteins:
            protein = proteins[key]
            above_threshold = False
            for pos in range(len(protein.concs)):
                protein.concs[pos] = max(0.0, protein.concs[pos] - protein.concs[pos] * Config.decay_factor)
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
