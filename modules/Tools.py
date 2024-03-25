from math import *
import numba as nb

@nb.njit
def HPM(x):
    """An accurate but quick approximation of the gaussian cumulative density function. 
        See Vazquez-Leal et al. Math Prob Eng 124029 (2012)"""   
    return 1. / (exp(-358 * x / 23 + 111 * atan(37 * x / 294)) + 1)

def sanitizeName(raw):
    out = raw.replace("/", "-").replace("?","-qst-")
    out = out.replace(":",".").replace("\"","'")
    out = out.replace("\\","-").replace("<","lt")
    out = out.replace(">","gt").replace("*","-st-")
    return out