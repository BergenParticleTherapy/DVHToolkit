from math import *
import numba as nb
import numpy as np

@nb.njit
def HPM(x):
    """An accurate but quick approximation of the gaussian cumulative density function. 
        See Vazquez-Leal et al. Math Prob Eng 124029 (2012)"""   
    return 1. / (exp(-358 * x / 23 + 111 * atan(37 * x / 294)) + 1)