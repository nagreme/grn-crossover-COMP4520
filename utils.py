from bitarray import bitarray
import numpy as np
import math

#some handy method for manipulating bitarrays and floating point values
class Utils():
    #creates a random bitvec of the given size
    @staticmethod
    def rand_bitvec(size):
        bv = bitarray()
        #this is awkward, but it's the way the library works...
        num_bytes = int(math.ceil(size / 8))
        bv.frombytes(np.random.bytes(num_bytes))
        if size % 8:
            bv = bv[:size]
        
        return bv

    #Checks to see if val is between low and high. If so, returns val.
    #If val < low, returns low. If val > high, returns high.
    @staticmethod
    def clamp(val, low, high):
        return min(max(val, low), high)

    #Compares two bitarrays x and y, and returns the number of bits that are the same
    #(have the same value in the same position)
    @staticmethod
    def same_count(x, y):
        return len(x) - (x ^ y).count()
