from math import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.optimize import basinhopping

from tkinter import *
from tkinter import ttk

from ..Tools import *


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

    class Result:
        def __init__(self, fun, x):
            self.fun = fun
            self.x = x

        def items(self):
            d = {'x': self.x, 'fun': self.fun}
            return d.items()

    return Result(toxArray[minIdx], self.bestParameters)


def doGradientOptimization(self, progress):
    class MyTakeStep(object):
        def __init__(self, options, stepsize=1):
            self.stepsize = stepsize
            self.options = options

        def __call__(self, x):
            s = self.stepsize
            if self.options.NTCPcalculation.get() == "Logit":
                x[0] += np.random.uniform(-s * self.options.basinHoppingAsize.get(), s * self.options.basinHoppingAsize.get())
                x[1] += np.random.uniform(-s * self.options.basinHoppingBsize.get(), s * self.options.basinHoppingBsize.get())
            elif self.options.NTCPcalculation.get() == "LKB":
                x[0] += np.random.uniform(-s * self.options.basinHoppingNsize.get(), s * self.options.basinHoppingNsize.get())
                x[1] += np.random.uniform(-s * self.options.basinHoppingMsize.get(), s * self.options.basinHoppingMsize.get())
                x[2] += np.random.uniform(-s * self.options.basinHoppingTD50size.get(), s * self.options.basinHoppingTD50size.get())

            eps = 1e-6
            for idx in range(len(x)):  # Don't push x beyond bounds
                x[idx] = max(x[idx], bounds[idx][0]) + eps
                x[idx] = min(x[idx], bounds[idx][1]) - eps

            return x

    def funLogitLS(x, *args):
        error = 0
        for tox, Dpercent, time in args:
            NTCP = 1 - 1 / (1 + exp(x[0] + x[1] * Dpercent))
            if time and not tox:
                NTCP *= (1 - exp(-(x[2] * time)**x[3]))
            error += (tox - NTCP) ** 2
        return error

    def funLKBLS(x, *args):
        error = 0
        for tox, GEUDspline, time in args:
            NTCP = HPM((GEUDspline(x[0]) - x[2]) / (x[1] * x[2]))
            if time and not tox:
                NTCP *= (1 - exp(-(x[3] * time)**x[4]))
            error += (tox - NTCP) ** 2
        return error

    def funLogitLLH(x, *args):
        error = 0
        for tox, Dpercent, time in args:
            NTCP = 1 - 1 / (1 + exp(x[0] + x[1] * Dpercent))
            if time and not tox:
                NTCP *= (1 - exp(-(x[2] * time)**x[3]))

            if tox:
                error -= log(max(NTCP, 1e-323))
            else:
                error -= log(max(1 - NTCP, 1e-323))

        return error

    def funLKBLLH(x, *args):
        error = 0
        for tox, GEUDspline, time in args:
            NTCP = HPM((GEUDspline(x[0]) - x[2]) / (x[1] * x[2]))
            if time and not tox:
                NTCP *= (1 - exp(-(x[3] * time)**x[4]))

            if tox:
                error -= log(max(NTCP, 1e-323))
            else:
                error -= log(max(1 - NTCP, 1e-323))

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
            NTCPTime = self.NTCPTimeDict and self.NTCPTimeDict[name.split("_")[0].split("tox")[0]] / 12 or None
            argTuple += ((patient.getTox() >= self.options.toxLimit.get(), patient.getDpercent(), NTCPTime),)

    elif self.options.NTCPcalculation.get() == "LKB":
        for name, patient in list(self.patients.items()):
            NTCPTime = self.NTCPTimeDict and self.NTCPTimeDict[name.split("_")[0].split("tox")[0]] / 12 or None
            argTuple += ((patient.getTox() >= self.options.toxLimit.get(), patient.getGEUD, NTCPTime),)

    mytakestep = MyTakeStep(self.options)
    if self.options.NTCPcalculation.get() == "Logit":
        if not self.NTCPTimeDict:
            bounds = ((self.options.aFrom.get(), self.options.fixA.get() and self.options.aFrom.get() or self.options.aTo.get()),
                      (self.options.bFrom.get(), self.options.fixB.get() and self.options.bFrom.get() or self.options.bTo.get()))
        else:
            bounds = ((self.options.aFrom.get(), self.options.fixA.get() and self.options.aFrom.get() or self.options.aTo.get()),
                      (self.options.bFrom.get(), self.options.fixB.get() and self.options.bFrom.get() or self.options.bTo.get()),
                      (self.options.lambdaFrom.get(), self.options.fixLambda.get() and self.options.lambdaFrom.get() or self.options.lambdaTo.get()),
                      (self.options.gammaFrom.get(), self.options.fixGamma.get() and self.options.gammaFrom.get() or self.options.gammaTo.get()))

        if len(self.bestParameters):
            x0 = np.array(self.bestParameters[:2])
        else:
            x0 = np.mean(bounds, axis=1)

    else:
        if not self.NTCPTimeDict:
            bounds = ((self.options.nFrom.get(), self.options.fixN.get() and self.options.nFrom.get() or self.options.nTo.get()),
                      (self.options.mFrom.get(), self.options.fixM.get() and self.options.mFrom.get() or self.options.mTo.get()),
                      (self.options.TD50From.get(), self.options.fixTD50.get() and self.options.TD50From.get() or self.options.TD50To.get()))
        else:
            bounds = ((self.options.nFrom.get(), self.options.fixN.get() and self.options.nFrom.get() or self.options.nTo.get()),
                      (self.options.mFrom.get(), self.options.fixM.get() and self.options.mFrom.get() or self.options.mTo.get()),
                      (self.options.TD50From.get(), self.options.fixTD50.get() and self.options.TD50From.get() or self.options.TD50To.get()),
                      (self.options.lambdaFrom.get(), self.options.fixLambda.get() and self.options.lambdaFrom.get() or self.options.lambdaTo.get()),
                      (self.options.gammaFrom.get(), self.options.fixGamma.get() and self.options.gammaFrom.get() or self.options.gammaTo.get()))

        if len(self.bestParameters):
            x0 = np.array(self.bestParameters)
        else:
            x0 = np.mean(bounds, axis=1)

    if self.options.optimizationMetric.get() == "LLH" and self.options.NTCPcalculation.get() == "LKB":
        fun = funLKBLLH
    elif self.options.optimizationMetric.get() == "LLH" and self.options.NTCPcalculation.get() == "Logit":
        fun = funLogitLLH
    elif self.options.optimizationMetric.get() == "LS" and self.options.NTCPcalculation.get() == "LKB":
        fun = funLKBLS
    elif self.options.optimizationMetric.get() == "LS" and self.options.NTCPcalculation.get() == "Logit":
        fun = funLogitLS
    else:
        raise Exception(f"No valid optimization function for optimization metric {self.options.optimizationMetric.get()} and NTCP={self.options.NTCPcalculation.get()}")
        fun = None

    res = basinhopping(fun, x0, niter=self.options.basinHoppingIterations.get(), T=self.options.basinHoppingTemperature.get(),
                       minimizer_kwargs={'args': argTuple, 'method': 'TNC', 'bounds': bounds}, take_step=mytakestep, callback=print_fun)

    self.bestParameters = res.x

    self.calculateNTCP()

    if progress:
        progress['value'] = 0

    res["TD5%"] = self.calculateTDxFromLogit(5)
    res["TD50%"] = self.calculateTDxFromLogit(50)

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
