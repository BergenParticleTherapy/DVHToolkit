from math import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.optimize import basinhopping

from tkinter import *

from ..Tools import *

class Result:
    def __init__(self, fun, x):
        self.fun = fun
        self.x = x

    def items(self):
        d = {'x': self.x, 'fun': self.fun}
        return d.items()

def gradient_respecting_bounds(bounds, fun, eps=1e-8):
    """bounds: list of tuples (lower, upper)
        FROM https://stackoverflow.com/questions/52208363/scipy-minimize-violates-given-bounds"""

    def gradient(x):
        fx = fun(x)
        grad = np.zeros(len(x))
        for k in range(len(x)):
            d = np.zeros(len(x))
            d[k] = eps if x[k] + eps <= bounds[k][1] else -eps
            grad[k] = (fun(x + d) - fx) / d[k]
        return grad
    return gradient


def doMatrixMinimization(self, progress):
    """Original method from Piero Fossati (rewritten from Matlab)."""

    N = self.options.matrixSize.get()
    aList = np.linspace(self.options.aFrom.get(), self.options.aTo.get(), num=N)
    bList = np.linspace(self.options.bFrom.get(), self.options.bTo.get(), num=N)
    nList = np.linspace(self.options.nFrom.get(), self.options.nTo.get(), num=N)
    mList = np.linspace(self.options.mFrom.get(), self.options.mTo.get(), num=N)
    TD50List = np.linspace(self.options.TD50From.get(), self.options.TD50To.get(), num=N)

    if self.options.fixA.get():
        aList = aList[0:1]
    if self.options.fixB.get():
        bList = bList[0:1]
    if self.options.fixN.get():
        nList = nList[0:1]
    if self.options.fixM.get():
        mList = mList[0:1]
    if self.options.fixTD50.get():
        TD50List = TD50List[0:1]

    if self.options.NTCPcalculation.get() == "LKB":
        params = [nList, mList, TD50List]
        if progress:
            progress['maximum'] = len(nList)
    else:
        params = [aList, bList, [1]]
        if progress:
            progress['maximum'] = len(aList)

    toxArray = np.zeros([len(k) for k in params], dtype=np.longdouble)
    argTuple = ()
    for name, patient in list(self.patients.items()):
        argTuple += ((patient.getTox() >= self.options.toxLimit.get(), patient.getGEUD),)

    for idx1, val1 in enumerate(params[0]):
        if progress:
            progress.step(1)
            progress.update_idletasks()

        for idx2, val2 in enumerate(params[1]):
            for idx3, val3 in enumerate(params[2]):  # [1] if Logit is used
                error = 0

                for tox, dose in argTuple:
                    if self.options.NTCPcalculation.get() == "LKB":
                        n = val1
                        m = val2
                        TD50 = val3
                        gEUD = dose(n)
                        NTCP = HPM((gEUD - TD50) / (m * TD50))
                    else:
                        a = val1
                        b = val2
                        Dpercent = dose
                        NTCP = 1 - 1 / (1 + exp(a + b * Dpercent))

                    if self.options.optimizationMetric.get() == "LS":
                        error += (tox - NTCP) ** 2
                    else:
                        if tox:
                            if NTCP > 0:
                                error -= log(NTCP)
                            else:
                                error += 2500  # assume minimum probability of ~10^1000
                        else:
                            if NTCP < 1:
                                error -= log(1 - NTCP)
                            else:
                                error += 2500

                toxArray[idx1, idx2, idx3] = error

    minIdx = np.unravel_index(np.argmin(toxArray), np.shape(toxArray))
    self.bestParameters = [params[0][minIdx[0]], params[1][minIdx[1]], params[2][minIdx[2]]]
    if self.options.NTCPcalculation.get() == "Logit":
        self.bestParameters.pop()

    self.calculateNTCP()

    if progress:
        progress['value'] = 0

    return Result(toxArray[minIdx], self.bestParameters)


def doGradientOptimization(self, progress):
    class MyTakeStep(object):
        def __init__(self, options, idx, stepsize=1):
            self.stepsize = stepsize
            self.options = options
            self.idx = idx

        def __call__(self, x):
            s = self.stepsize

            if self.options.NTCPcalculation.get() == "Logit":
                if not self.options.fixA.get():
                    x[self.idx['a']] += np.random.uniform(-s * self.options.basinHoppingAsize.get(), s * self.options.basinHoppingAsize.get())
                if not self.options.fixB.get():
                    x[self.idx['b']] += np.random.uniform(-s * self.options.basinHoppingBsize.get(), s * self.options.basinHoppingBsize.get())
            elif self.options.NTCPcalculation.get() == "LKB":
                if not self.options.fixN.get():
                    x[self.idx['n']] += np.random.uniform(-s * self.options.basinHoppingNsize.get(), s * self.options.basinHoppingNsize.get())
                if not self.options.fixM.get():
                    x[self.idx['m']] += np.random.uniform(-s * self.options.basinHoppingMsize.get(), s * self.options.basinHoppingMsize.get())
                if not self.options.fixTD50.get():
                    x[self.idx['TD50']] += np.random.uniform(-s * self.options.basinHoppingTD50size.get(), s * self.options.basinHoppingTD50size.get())

            eps = 1e-6
            for idx in range(len(x)):  # Don't push x beyond bounds
                x[idx] = max(x[idx], bounds[idx][0]) - eps
                x[idx] = min(x[idx], bounds[idx][1]) + eps

            return x

    def funLogitLS(x, *args):
        error = 0

        a = self.options.fixN.get() and p['n'] or x[self.idx['a']]
        b = self.options.fixM.get() and p['m'] or x[self.idx['b']]
        if self.NTCPTimeDict:
            Lambda = self.options.fixLambda.get() and p['lambda'] or x[self.idx['lambda']]
            gamma = self.options.fixGamma.get() and p['gamma'] or x[self.idx['gamma']]
        else:
            Lambda = gamma = None

        for tox, Dpercent, time in args:
            NTCP = 1 - 1 / (1 + exp(a + b * Dpercent))
            if time and not tox:
                NTCP *= (1 - exp(-(Lambda * time)**gamma))
            error += (tox - NTCP) ** 2

        # Lower and Upper bounds to optimization ("all cases below X Gy should be negative")
        if self.options.NTCPBoundWeight.get() and self.options.NTCPUseBound.get():
            NTCP = 1 - 1 / (1 + exp(a + b * self.options.NTCPBoundLower.get()))
            error += len(args) * self.options.NTCPBoundWeight.get() * (NTCP ** 2)

            NTCP = 1 - 1 / (1 + exp(a + b * self.options.NTCPBoundHigher.get()))
            error += len(args) * self.options.NTCPBoundWeight.get() * (1 - NTCP)**2

        return error

    def funLKBLS(x, *args):
        error = 0

        n = self.options.fixN.get() and p['n'] or x[self.idx['n']]
        m = self.options.fixM.get() and p['m'] or x[self.idx['m']]
        TD50 = self.options.fixTD50.get() and p['TD50'] or x[self.idx['TD50']]
        if self.NTCPTimeDict:
            Lambda = self.options.fixLambda.get() and p['lambda'] or x[self.idx['lambda']]
            gamma = self.options.fixGamma.get() and p['gamma'] or x[self.idx['gamma']]
        else:
            Lamda = gamma = None

        for tox, GEUDspline, time in args:
            gEUD = GEUDspline(x[0])
            if gEUD < 0:
                NTCP = 0
            else:
                NTCP = HPM((GEUDspline(n) - TD50) / (m * TD50))
            if time and not tox:
                NTCP *= (1 - exp(-(Lambda * time)**gamma))
            error += (tox - NTCP) ** 2

        # Lower and Upper bounds to optimization ("all cases below X Gy should be negative")
        if self.options.NTCPBoundWeight.get() and self.options.NTCPUseBound.get():
            NTCP = HPM((self.options.NTCPBoundLower.get() - TD50) / (m * TD50))
            error += len(args) * self.options.NTCPBoundWeight.get() * (NTCP ** 2)

            NTCP = HPM((self.options.NTCPBoundUpper.get() - TD50) / (m * TD50))
            error += len(args) * self.options.NTCPBoundWeight.get() * (1 - NTCP)**2

        return error

    def funLogitLLH(x, *args):
        error = 0

        a = self.options.fixN.get() and p['n'] or x[self.idx['a']]
        b = self.options.fixM.get() and p['m'] or x[self.idx['b']]
        eps = 1e-20

        for tox, Dpercent in args:
            NTCP = 1 - 1 / (1 + exp(a + b * Dpercent))

            if tox:
                error -= log(NTCP + eps)
            else:
                error -= log(1 - NTCP + eps)

        # Lower and Upper bounds to optimization ("all cases below X Gy should be negative")
        if self.options.NTCPBoundWeight.get() and self.options.NTCPUseBound.get():
            NTCP = 1 - 1 / (1 + exp(a + b * self.options.NTCPBoundLower.get()))
            error -= len(args) * self.options.NTCPBoundWeight.get() * log(max(1 - NTCP, 1e-323))

            NTCP = 1 - 1 / (1 + exp(a + b * self.options.NTCPBoundHigher.get()))
            error -= len(args) * self.options.NTCPBoundWeight.get() * log(max(NTCP, 1e-323))

        return error

    def funLogitLLHtime(x, *args):
        error = 0

        a = self.options.fixN.get() and p['n'] or x[self.idx['a']]
        b = self.options.fixM.get() and p['m'] or x[self.idx['b']]
        Lambda = self.options.fixLambda.get() and p['lambda'] or x[self.idx['lambda']]
        gamma = self.options.fixGamma.get() and p['gamma'] or x[self.idx['gamma']]
        eps = 1e-20

        for tox, Dpercent, time in args:
            NTCP = 1 - 1 / (1 + exp(a + b * Dpercent))

            if tox:
                error -= log(max(NTCP + eps))
            else:
                if time:
                    NTCP *= (1 - exp(-(Lambda * time)**gamma))

                error -= log(max(1 - NTCP + eps))

        # Lower and Upper bounds to optimization ("all cases below X Gy should be negative")
        if self.options.NTCPBoundWeight.get() and self.options.NTCPUseBound.get():
            NTCP = 1 - 1 / (1 + exp(a + b * self.options.NTCPBoundLower.get()))
            error -= len(args) * self.options.NTCPBoundWeight.get() * log(max(1 - NTCP, 1e-323))

            NTCP = 1 - 1 / (1 + exp(a + b * self.options.NTCPBoundHigher.get()))
            error -= len(args) * self.options.NTCPBoundWeight.get() * log(max(NTCP, 1e-323))

        return error

    def funLKBLLH(x, *args):
        error = 0

        n = self.options.fixN.get() and p['n'] or x[self.idx['n']]
        m = self.options.fixM.get() and p['m'] or x[self.idx['m']]
        TD50 = self.options.fixTD50.get() and p['TD50'] or x[self.idx['TD50']]
        eps = 1e-20

        for tox, GEUDspline in args:
            gEUD = GEUDspline(n)
            
            if gEUD <= 0:
                NTCP = 0
            else:
                NTCP = HPM((gEUD - TD50) / (m * TD50))

            if tox:
                error -= log(NTCP + eps)
            else:
                error -= log(1-NTCP + eps)
            

        # Lower and Upper bounds to optimization ("all cases below X Gy should be negative")
        if self.options.NTCPBoundWeight.get() and self.options.NTCPUseBound.get():
            NTCP = HPM((self.options.NTCPBoundLower.get() - TD50) / (m * TD50))
            error -= len(args) * self.options.NTCPBoundWeight.get() * log(1 - NTCP + eps)

            NTCP = HPM((self.options.NTCPBoundUpper.get() - TD50) / (m * TD50))
            error -= len(args) * self.options.NTCPBoundWeight.get() * log(NTCP + eps)

        return error

    def funLKBLLHtime(x, *args):
        error = 0

        n = self.options.fixN.get() and p['n'] or x[self.idx['n']]
        m = self.options.fixM.get() and p['m'] or x[self.idx['m']]
        TD50 = self.options.fixTD50.get() and p['TD50'] or x[self.idx['TD50']]

        Lambda = self.options.fixLambda.get() and p['lambda'] or x[self.idx['lambda']]
        gamma = self.options.fixGamma.get() and p['gamma'] or x[self.idx['gamma']]
        eps = 1e-20

        for tox, GEUDspline, time in args:
            gEUD = GEUDspline(n)
            if gEUD <= 0:
                NTCP = 0
            else:
                NTCP = HPM((gEUD - TD50) / (m * TD50))

            if tox:
                error -= np.log(NTCP + eps)
            else:
                if time: 
                    NTCP *= (1 - exp(-(Lambda * time)**gamma))
                error -= np.log(1 - NTCP + eps)

        # Lower and Upper bounds to optimization ("all cases below X Gy should be negative")
        if self.options.NTCPBoundWeight.get() and self.options.NTCPUseBound.get():
            NTCP = HPM((self.options.NTCPBoundLower.get() - TD50) / (m * TD50))
            error -= len(args) * self.options.NTCPBoundWeight.get() * log(max(1 - NTCP, 1e-323))

            NTCP = HPM((self.options.NTCPBoundUpper.get() - TD50) / (m * TD50))
            error -= len(args) * self.options.NTCPBoundWeight.get() * log(max(NTCP, 1e-323))

        return error

    def print_fun(x, f, accepted):
        if progress:
            progress.step(1)
            progress.update_idletasks()

    if progress:
        progress['maximum'] = self.options.basinHoppingIterations.get()

    argTuple = ()
    if self.options.NTCPcalculation.get() == "Logit":
        for name, patient in list(self.patients.items()):
            if self.NTCPTimeDict:
                NTCPTime = self.NTCPTimeDict and self.NTCPTimeDict[name.split("_")[0].split("tox")[0]] / 12 or None
                argTuple += ((patient.getTox() >= self.options.toxLimit.get(), patient.getDpercent(), NTCPTime),)
            else:
                argTuple += ((patient.getTox() >= self.options.toxLimit.get(), patient.getDpercent()),)

    elif self.options.NTCPcalculation.get() == "LKB":
        for name, patient in list(self.patients.items()):
            if self.NTCPTimeDict:
                NTCPTime = self.NTCPTimeDict and self.NTCPTimeDict[name.split("_")[0].split("tox")[0]] / 12 or None
                argTuple += ((patient.getTox() >= self.options.toxLimit.get(), patient.fastGetGEUD, NTCPTime),)
            else:
                argTuple += ((patient.getTox() >= self.options.toxLimit.get(), patient.fastGetGEUD),)

    mytakestep = MyTakeStep(self.options, self.idx)

    if (self.options.NTCPcalculation.get() == "LKB" and self.options.fixN.get() and self.options.fixM.get() and self.options.fixTD50.get()) \
        or (self.options.NTCPcalculation.get() == "Logit" and self.options.fixA.get() and self.options.fixB.get()):
        # Fixed parameters, so only return LLH
        # E.g. for external model validation on dataset

        funs = {'LKB': {'LLH': funLKBLLH, 'LS': funLKBLS}, 'Logit': {'LLH': funLogitLLH, 'LS': funLogitLS}}
        fun = funs[self.options.NTCPcalculation.get()][self.options.optimizationMetric.get()]

        if self.NTCPTimeDict and self.options.NTCPcalculation.get() == "LKB":
            fun = funLKBLLHtime

        p = {'a': self.options.aFrom.get(), 'b': self.options.bFrom.get(),
             'n': self.options.nFrom.get(), 'm': self.options.mFrom.get(), 'TD50': self.options.TD50From.get()}

        params = self.options.NTCPcalculation.get() == "LKB" and [p['n'], p['m'], p['TD50']] or [p['a'], p['b']]

        error = fun(params, *argTuple)
        return Result(error, params)

    if self.options.NTCPcalculation.get() == "Logit":
        boundsA = (self.options.aFrom.get(), self.options.aTo.get())
        boundsB = (self.options.bFrom.get(), self.options.bTo.get())
        boundsLambda = (self.options.lambdaFrom.get(), self.options.lambdaTo.get())
        boundsGamma = (self.options.gammaFrom.get(), self.options.gammaTo.get())
        bounds = list()

        if not self.options.fixA.get():
            bounds.append(boundsA)
        if not self.options.fixB.get():
            bounds.append(boundsB)

        if self.NTCPTimeDict:
            if not self.options.fixLambda.get():
                bounds.append(boundsLambda)
            if not self.options.fixGamma.get():
                bounds.append(boundsGamma)

        p = dict()
        p['a'] = self.options.aFrom.get()
        p['b'] = self.options.bFrom.get()
        p['lambda'] = self.options.lambdaFrom.get()
        p['gamma'] = self.options.gammaFrom.get()

        x0 = np.mean(bounds, axis=1)

    else:  # FOR LKB
        boundsN = (self.options.nFrom.get(), self.options.nTo.get())
        boundsM = (self.options.mFrom.get(), self.options.mTo.get())
        boundsTD50 = (self.options.TD50From.get(), self.options.TD50To.get())
        boundsLambda = (self.options.lambdaFrom.get(), self.options.lambdaTo.get())
        boundsGamma = (self.options.gammaFrom.get(), self.options.gammaTo.get())
        bounds = list()

        if not self.options.fixN.get():
            bounds.append(boundsN)
        if not self.options.fixM.get():
            bounds.append(boundsM)
        if not self.options.fixTD50.get():
            bounds.append(boundsTD50)

        if self.NTCPTimeDict:
            if not self.options.fixLambda.get():
                bounds.append(boundsLambda)
            if not self.options.fixGamma.get():
                bounds.append(boundsGamma)

        x0 = list()
        p = dict()

        p['n'] = self.options.nFrom.get()
        p['m'] = self.options.mFrom.get()
        p['TD50'] = self.options.TD50From.get()
        p['lambda'] = self.options.lambdaFrom.get()
        p['gamma'] = self.options.gammaFrom.get()

        x0 = np.mean(bounds, axis=1)

    funs = {'LKB': {'LLH': funLKBLLH, 'LS': funLKBLS}, 'Logit': {'LLH': funLogitLLH, 'LS': funLogitLS}}
    fun = funs[self.options.NTCPcalculation.get()][self.options.optimizationMetric.get()]

    if self.NTCPTimeDict and self.options.NTCPcalculation.get() == "LKB":
        fun = funLKBLLHtime
    elif self.NTCPTimeDict and self.options.NTCPcalculation.get() == "Logit":
        fun = funLogitLLHtime

    res = basinhopping(fun, x0, niter=self.options.basinHoppingIterations.get(), T=self.options.basinHoppingTemperature.get(),
                       minimizer_kwargs={'args': argTuple, 'method': 'Powell', 'bounds': bounds},
                       take_step=mytakestep, callback=print_fun)

    self.bestParameters = res.x

    self.calculateNTCP()

    if progress:
        progress['value'] = 0

    if self.options.NTCPcalculation.get() == "Logit":
        res["TD5%"] = self.calculateTDxFromLogit(5)
        res["TD50%"] = self.calculateTDxFromLogit(50)
    else:
        res["TD50%"] = None
        res["TD5%"] = None

    return res


def profileLikelihood(self):
    """For some reason this article gets cited for this method:
        Stryhn, H, and J Christensen. “Confidence Intervals by the Profile Likelihood Method, 
        with Applications in Veterinary Epidemiology,” 2003. http://www.sciquest.org.nz/elibrary/edition/5008."""

    q = 1 - self.options.confidenceIntervalPercent.get() / 100
    gamma = st.chi2.isf(q, df=1) / 2

    if self.options.NTCPcalculation.get() == "Dpercent":
        aInit = self.bestParameters[0]
        aValues = np.arange(aInit * 0.5, aInit * 2, aInit / 100)
        bInit = self.bestParameters[1]
        bValues = np.arange(bInit * 0.5, bInit * 2, bInit / 100)

        LLH_a = np.zeros(len(aValues))
        LLH_b = np.zeros(len(bValues))

        for idx, a in enumerate(aValues):
            error = 0
            for patient in list(self.patients.values()):
                tox = patient.getTox() >= self.options.toxLimit.get()
                NTCP = 1 - 1 / (1 + exp(a + bInit * patient.getDpercent()))
                if tox:
                    if NTCP > 0:
                        error -= log(NTCP)
                    else:
                        error += 2500  # assume minimum probability of 10^-1000
                elif not tox:
                    if NTCP < 1:
                        error -= log(1 - NTCP)
                    else:
                        error += 2500  # assume minimum probability of 10^-1000

            LLH_a[idx] = -error

        aCIvalue = np.max(LLH_a) - gamma
        idx = np.argwhere(np.diff(np.sign(LLH_a - aCIvalue))).flatten()
        aLow = min(aValues[idx[0]], aValues[idx[1]])
        aHigh = max(aValues[idx[0]], aValues[idx[1]])

        for idx, b in enumerate(bValues):
            error = 0
            for patient in list(self.patients.values()):
                tox = patient.getTox() >= self.options.toxLimit.get()
                NTCP = 1 - 1 / (1 + exp(aInit + b * patient.getDpercent()))
                if tox:
                    if NTCP > 0:
                        error -= log(NTCP)
                    else:
                        error += 2500  # assume minimum probability of 10^-1000
                elif not tox:
                    if NTCP < 1:
                        error -= log(1 - NTCP)
                    else:
                        error += 2500  # assume minimum probability of 10^-1000
            LLH_b[idx] = -error

        bCIvalue = np.max(LLH_b) - gamma
        idx = np.argwhere(np.diff(np.sign(LLH_b - bCIvalue))).flatten()
        bLow = min(bValues[idx[0]], bValues[idx[1]])
        bHigh = max(bValues[idx[0]], bValues[idx[1]])

        fig, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(20, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        plt.title(f"{self.options.confidenceIntervalPercent.get()}% profile likelihood scan over NTCP = 1 - 1 / (1 + exp(a + b * D{self.options.NTCPcalculationDpercent.get()}%))")
        ax1.plot(aValues, LLH_a, label="Log Likelihood")
        ax1.plot([aValues[0], aValues[-1]], [aCIvalue, aCIvalue], "-", label=f"chi2(1)/2 = {gamma:.2f} reduction in LLH")
        ax1.set_xlabel("a values")
        ax1.set_ylabel("log likelihood")
        # ax1.set_xscale('log')
        ax1.legend(loc='lower right')
        ax2.plot(bValues, LLH_b, label="Log likelihood")
        ax2.plot([bValues[0], bValues[-1]], [bCIvalue, bCIvalue], "-", label=f"chi2(1)/2 = {gamma:.2f} reduction in LLH")
        ax2.set_xlabel("b values")
        ax2.set_ylabel("log likelihood")
        # ax2.set_xscale('log')
        ax2.legend(loc='lower right')

        plt.savefig("Output/profileLikelihood.png")

        return [[aLow, aHigh], [bLow, bHigh]]

    if self.options.NTCPcalculation.get() == "LKB":
        nInit = self.bestParameters[0]
        nValues = np.arange(self.options.nFrom.get(), self.options.nTo.get(), nInit / 50)
        mInit = self.bestParameters[1]
        mValues = np.arange(mInit * 0.2, mInit * 3, mInit / 100)
        TD50Init = self.bestParameters[2]
        TD50Values = np.arange(TD50Init * 0.5, TD50Init * 3, TD50Init / 100)

        LLH_n = np.zeros(len(nValues))
        LLH_m = np.zeros(len(mValues))
        LLH_TD50 = np.zeros(len(TD50Values))

        for idx, n in enumerate(nValues):  # scan over n
            error = 0
            for patient in list(self.patients.values()):
                tox = patient.getTox() >= self.options.toxLimit.get()
                NTCP = HPM((patient.getGEUD(n) - TD50Init) / (mInit * TD50Init))
                if tox:
                    if NTCP > 0:
                        error -= log(NTCP)
                    else:
                        error += 2500  # assume minimum probability of 10^-1000
                elif not tox:
                    if NTCP < 1:
                        error -= log(1 - NTCP)
                    else:
                        error += 2500  # assume minimum probability of 10^-1000

            LLH_n[idx] = -error

        nCIvalue = np.max(LLH_n) - gamma
        idx = np.argwhere(np.diff(np.sign(LLH_n - nCIvalue))).flatten()
        nLow = nValues[idx[0]]
        nHigh = nValues[idx[1]]

        for idx, m in enumerate(mValues):  # Scan over m
            error = 0
            for patient in list(self.patients.values()):
                tox = patient.getTox() >= self.options.toxLimit.get()
                NTCP = HPM((patient.getGEUD(nInit) - TD50Init) / (m * TD50Init))
                if tox:
                    if NTCP > 0:
                        error -= log(NTCP)
                    else:
                        error += 2500  # assume minimum probability of 10^-1000
                elif not tox:
                    if NTCP < 1:
                        error -= log(1 - NTCP)
                    else:
                        error += 2500  # assume minimum probability of 10^-1000

            LLH_m[idx] = -error

        mCIvalue = np.max(LLH_m) - gamma
        idx = np.argwhere(np.diff(np.sign(LLH_m - mCIvalue))).flatten()
        mLow = mValues[idx[0]]
        mHigh = mValues[idx[1]]

        for idx, TD50 in enumerate(TD50Values):  # Scan over TD50
            error = 0
            for patient in list(self.patients.values()):
                tox = patient.getTox() >= self.options.toxLimit.get()
                NTCP = HPM((patient.getGEUD(nInit) - TD50) / (mInit * TD50))
                if tox:
                    if NTCP > 0:
                        error -= log(NTCP)
                    else:
                        error += 2500  # assume minimum probability of 10^-1000
                elif not tox:
                    if NTCP < 1:
                        error -= log(1 - NTCP)
                    else:
                        error += 2500  # assume minimum probability of 10^-1000

            LLH_TD50[idx] = -error

        TD50CIvalue = np.max(LLH_TD50) - gamma
        idx = np.argwhere(np.diff(np.sign(LLH_TD50 - TD50CIvalue))).flatten()
        TD50Low = TD50Values[idx[0]]
        TD50High = TD50Values[idx[1]]

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=False, sharey=False, figsize=(20, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        plt.title(f"{self.options.confidenceIntervalPercent.get()}% profile likelihood scan over NTCP = LKB(n,m,TD50)")
        ax1.plot(nValues, LLH_n, label="Log Likelihood")
        ax1.plot([nValues[0], nValues[-1]], [nCIvalue, nCIvalue], "-", label=f"chi2(1)/2 = {gamma:.2f} reduction in LLH")
        ax1.set_xlabel("n values")
        ax1.set_ylabel("log likelihood")
        ax1.legend(loc='lower right')
        ax2.plot(mValues, LLH_m, label="Log likelihood")
        ax2.plot([mValues[0], mValues[-1]], [mCIvalue, mCIvalue], "-", label=f"chi2(1)/2 = {gamma:.2f} reduction in LLH")
        ax2.set_xlabel("m values")
        ax2.set_ylabel("log likelihood")
        ax2.legend(loc='lower right')
        ax3.plot(TD50Values, LLH_TD50, label="Log likelihood")
        ax3.plot([TD50Values[0], TD50Values[-1]], [TD50CIvalue, TD50CIvalue], "-", label=f"chi2(1)/2 = {gamma:.2f} reduction in LLH")
        ax3.set_xlabel("TD50 values")
        ax3.set_ylabel("log likelihood")
        ax3.legend(loc='lower right')

        plt.savefig("Output/profileLikelihood.png")

        return [[nLow, nHigh], [mLow, mHigh], [TD50Low, TD50High]]
