from math import *
import numpy as np
import matplotlib.pyplot as plt

from ..Tools import *


class ParameterSpace:
    def __init__(self, nIterations, model, log, cohort):
        self.model = model

        if self.model == "Logit":
            self.parameters = {'a': None, 'b': None}

        else:
            self.parameters = {'n': None, 'm': None, 'TD50': None}

        self.parameterSpace = dict()
        for p in self.parameters.keys():
            self.parameterSpace[p] = np.zeros(nIterations)

        self.LLH = np.zeros(nIterations)

        self.CI = dict()

        self.idx = 0
        self.idxLLH = 0
        self.bootstrapCorrectionMethod = None
        self.log = log
        self.percentile = None
        self.upperPercent = None
        self.lowerPercent = None
        self.cohort = cohort

    def addPoint(self, x1, x2, x3=None):
        if self.model == "Logit":
            self.parameterSpace['a'][self.idx] = x1
            self.parameterSpace['b'][self.idx] = x2
        else:
            self.parameterSpace['n'][self.idx] = x1
            self.parameterSpace['m'][self.idx] = x2
            self.parameterSpace['TD50'][self.idx] = x3

        self.idx += 1

    def print2(self, string):
        self.log(string)
        print(string)

    def addPointLLH(self, LLH):
        self.LLH[self.idxLLH] = LLH
        self.idxLLH += 1

    def setParameters(self, paramDict):
        for k, v in paramDict.items():
            self.parameters[k] = v

    def getParameters(self):
        if self.model == "Logit":
            return [self.parameters[k] for k in ['a', 'b']]
        else:
            return [self.parameters[k] for k in ['n', 'm', 'TD50']]

    def calculateCI(self):
        for k, v in self.parameters.items():
            if v:
                self.CI[k] = [np.percentile(self.parameterSpace[k], perc) for perc in [self.lowerPercent, self.upperPercent]]

        if self.model == "Logit":
            self.parameterSpace['TD50'] = -self.parameterSpace['a'] / self.parameterSpace['b']
            self.CI["TD50"] = [np.percentile(self.parameterSpace["TD50"], perc) for perc in [self.lowerPercent, self.upperPercent]]

            self.parameterSpace['TD5'] = -(log(19) + self.parameterSpace['a']) / self.parameterSpace['b']
            self.CI["TD5"] = [np.percentile(self.parameterSpace["TD5"], perc) for perc in [self.lowerPercent, self.upperPercent]]

    def printCI(self):
        for k, v in self.parameters.items():
            if v:
                self.print2(f"{k} = {v} ({self.CI[k][0]:.3f} - {self.CI[k][1]:.3f})")

        if self.model == "Logit":
            parameters = {**self.parameters, **{"TD50": -self.parameters['a'] / self.parameters['b']}}
            parameters['TD5'] = -(log(19) + self.parameters['a']) / self.parameters['b']
            for k in ["TD5", "TD50"]:
                self.print2(f"{k} = {parameters[k]} ({self.CI[k][0]:.3f} - {self.CI[k][1]:.3f}))")

    def trim(self):
        for p in self.parameters.keys():
            self.parameterSpace[p] = np.trim_zeros(self.parameterSpace[p])

        self.LLH = np.trim_zeros(self.LLH)

    def setPercentile(self, percentile):
        self.percentile = percentile
        self.lowerPercent = (100 - self.percentile) / 2
        self.upperPercent = 100 - self.lowerPercent

    def setBootstrapCorrectionMethod(self, method):
        self.bootstrapCorrectionMethod = method

    def applyPivot(self):
        plist = list(self.parameters.keys()) + ["TD50"]
        if self.model == "Logit":
            plist += ["TD5"]

        for p in self.parameters.keys():
            self.CI[p] = [2 * self.parameters[p] - self.CI[p][1], 2 * self.parameters[p] - self.CI[p][0]]
            lim = [np.min(self.parameterSpace[p]), np.max(self.parameterSpace[p])]

            if self.CI[p][0] < lim[0]:
                self.CI[p] = [k + lim[0] - self.CI[p][0] for k in self.CI[p]]

            if self.CI[p][1] > lim[1]:
                self.CI[p] = [k - (self.CI[p][1] - lim[1]) for k in self.CI[p]]

        correctionStatistics = {'mean': {p: np.mean(self.parameterSpace[p]) for p in self.parameters.keys()},
                                'median': {p: np.median(self.parameterSpace[p]) for p in self.parameters.keys()},
                                'none': {k: v for k, v in self.parameters.items()}}

        self.parameters = {k: 2 * v - correctionStatistics[self.bootstrapCorrectionMethod][k] for k, v in self.parameters.items()}

    def printResults(self, log):
        pString = ", ".join([f"{k} = {v:.3f}" for k, v in self.parameters.items()])
        self.print2(f"Using the {self.bootstrapCorrectionMethod} pivot bias correction method, the corrected best fits were {pString}.")

        self.print2(f"The bias corrected bootstrapped confidence intervals were:")
        for p in self.parameters.keys():
            self.print2(f"{p}\t= [{self.CI[p][0]:.2f}, {self.CI[p][1]:.2f}]")

    def plotResults(self):
        nPlots = 6  # self.model == "LKB" and 6 or 5
        fig, axs = plt.subplots(1, nPlots, sharex=False, sharey=False, figsize=(4 * nPlots, 5))

        plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        plt.ylabel("Frequency")

        ci_95 = dict()
        ci_68 = dict()

        ps_95 = dict()
        ps_68 = dict()

        filter_95 = np.ones(len(self.parameterSpace["TD50"]), dtype=bool)  # p is last valid parameter
        filter_68 = np.ones(len(self.parameterSpace["TD50"]), dtype=bool)  # p is last valid parameter

        if self.model == "LKB":
            parameters = self.parameters
        else:
            # We want to plot TD50 as well
            parameters = {**self.parameters, **{"TD50": -self.parameters['a'] / self.parameters['b']}}
            parameters['TD5'] = -(log(19) + self.parameters['a']) / self.parameters['b']

        for idx, p in enumerate(parameters.keys()):
            axs[idx].hist(x=self.parameterSpace[p], bins=50)
            axs[idx].set_xlabel(f"{p} values")
            axs[idx].set_title(f"Parameter space for {p} ({self.cohort})")

            axs[idx].plot([self.CI[p][0]] * 2, axs[idx].get_ybound(), 'k-',
                          [self.CI[p][1]] * 2, axs[idx].get_ybound(), 'k-')

            ci_95[p] = [2 * parameters[p] - np.percentile(self.parameterSpace[p], k) for k in [97.5, 2.5]]
            ci_68[p] = [2 * parameters[p] - np.percentile(self.parameterSpace[p], k) for k in [84, 16]]

            filter_95 *= (self.parameterSpace[p] >= ci_95[p][0]) & (self.parameterSpace[p] <= ci_95[p][1])
            filter_68 *= (self.parameterSpace[p] >= ci_68[p][0]) & (self.parameterSpace[p] <= ci_68[p][1])

            print(p, ci_95[p])

        if self.model == "LKB":
            for secondPar in ['n', 'm']:
                idx += 1

                axs[idx].plot(self.parameterSpace["TD50"], self.parameterSpace[secondPar], 'o', c="red", zorder=0, label="Outside 95% LLH")
                axs[idx].plot(self.parameterSpace["TD50"][filter_95], self.parameterSpace[secondPar][filter_95], 'o', c="yellow", zorder=5, label="Within 95% LLH")
                axs[idx].plot(self.parameterSpace["TD50"][filter_68], self.parameterSpace[secondPar][filter_68], 'o', c="green", zorder=10, label="Within 68% LLH")
                axs[idx].legend()
                axs[idx].set_xlabel("TD50 parameter")
                axs[idx].set_ylabel(f"{secondPar} parameter")
                axs[idx].set_title(f"TD50 vs {secondPar} for {self.cohort}")

        else:
            idx += 1

            axs[idx].plot(self.parameterSpace["a"], self.parameterSpace["b"], 'o', c="red", zorder=0, label="Outside 95% LLH")
            axs[idx].plot(self.parameterSpace["a"][filter_95], self.parameterSpace["b"][filter_95], 'o', c="yellow", zorder=5, label="Within 95% LLH")
            axs[idx].plot(self.parameterSpace["a"][filter_68], self.parameterSpace["b"][filter_68], 'o', c="green", zorder=10, label="Within 68% LLH")
            axs[idx].legend()
            axs[idx].set_xlabel("a parameter")
            axs[idx].set_ylabel(f"b parameter")
            axs[idx].set_title(f"a vs b for {self.cohort}")

        idx += 1
        axs[idx].hist(x=self.LLH, bins=50)
        axs[idx].set_xlabel("Log Likelihood values")
        axs[idx].set_title(f"Log Likelihood for {self.cohort}")

    def writeToFile(self):
        with open(f"Output/bootstrapParameterSpace.csv", "w") as out:
            if self.model == "LKB":
                out.write("cohort,n,m,TD50\n")
                for idx in range(len(self.parameterSpace['n'])):
                    n = self.parameterSpace['n'][idx]
                    m = self.parameterSpace['m'][idx]
                    TD50 = self.parameterSpace['TD50'][idx]
                    out.write(f"{self.cohort},{n},{m},{TD50}\n")
            else:
                out.write("cohort,a,b,TD50\n")
                for idx in range(len(self.parameterSpace['a'])):
                    a = self.parameterSpace['a'][idx]
                    b = self.parameterSpace['b'][idx]
                    TD50 = self.parameterSpace['TD50'][idx]
                    out.write(f"{self.cohort},{a},{b},{TD50}\n")
