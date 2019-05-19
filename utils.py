from bitarray import bitarray
import numpy as np
import math

class Utils():
    @staticmethod
    def rand_bitvec(size):
        bv = bitarray()
        num_bytes = int(math.ceil(size / 8))
        bv.frombytes(np.random.bytes(num_bytes))
        if size % 8:
            bv = bv[:size]
        
        return bv

    @staticmethod
    def clamp(val, low, high):
        return min(max(val, low), high)

    @staticmethod
    def same_count(x, y):
        return len(x) - (x ^ y).count()
