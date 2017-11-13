import numpy as np
import random

def randomW(parray):
    """
    given an array of weights, returns an integer from the range(0, len(parray)), with
    probability given by the normalized weight in that position of the array
    """
    s = float(np.sum(parray))
    parray /= s

    r = random.random()
    P = 0.0
    for i in range(0, len(parray)):
        P += parray[i]
        if r <= P:
            return i
