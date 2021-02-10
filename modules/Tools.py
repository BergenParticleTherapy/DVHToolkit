from math import *


def HPM(x):
    """An accurate but quick approximation of the gaussian cumulative density function. 
        See Vazquez-Leal et al. Math Prob Eng 124029 (2012)"""
    try:
        return 1 / (exp(-358 * x / 23 + 111 * atan(37 * x / 294)) + 1)
    except OverflowError as e:
        print(f"OVERFLOW error {e} with x={x}")
        return 0
