from math import *
import numpy as np
import matplotlib.pyplot as plt

from ..Tools import *


class ParameterSpace:
    def __init__(self, nIterations, log, cohort, options, idxParam):
        self.options = options

        self.model = self.options.NTCPcalculation.get()
        self.idxParam = idxParam

        if self.model == "Logit":
            self.parameters = {'a': None, 'b': None}

        else:
            self.parameters = {'n': None, 'm': None, 'TD50': None}

        self.parameterSpace = dict()
        for p in self.parameters.keys():
            self.parameterSpace[p] = np.zeros(nIterations)

        self.originalLLH = np.zeros(nIterations)
        self.originalNagelkerke = np.zeros(nIterations)
        self.originalMcFadden = np.zeros(nIterations)

        self.trainingLLH = np.zeros(nIterations)
        self.trainingNagelkerke = np.zeros(nIterations)
        self.trainingMcFadden = np.zeros(nIterations)

        self.testLLH = np.zeros(nIterations)
        self.testNagelkerke = np.zeros(nIterations)
        self.testMcFadden = np.zeros(nIterations)

        self.CI = dict()

        self.idx = 0
        self.idxOriginalLLH = 0
        self.idxTrainingLLH = 0
        self.idxTestLLH = 0
        self.bootstrapCorrectionMethod = None
        self.log = log
        self.percentile = None
        self.upperPercent = None
        self.lowerPercent = None
        self.cohort = cohort
        self.nIsLinear = self.options.nIsLinear.get()

    def addPoint(self, x):
        if self.model == "Logit":
            self.parameterSpace['a'][self.idx] = self.options.fixA.get() and self.options.aFrom.get() or x[self.idxParam['a']]
            self.parameterSpace['b'][self.idx] = self.options.fixB.get() and self.options.bFrom.get() or x[self.idxParam['b']]
        else:
            self.parameterSpace['n'][self.idx] = self.options.fixN.get() and self.options.nFrom.get() or x[self.idxParam['n']]
            self.parameterSpace['m'][self.idx] = self.options.fixM.get() and self.options.mFrom.get() or x[self.idxParam['m']]
            self.parameterSpace['TD50'][self.idx] = self.options.fixTD50.get() and self.options.TD50From.get() or x[self.idxParam['TD50']]

        self.idx += 1

    def print2(self, string):
        self.log(string)
        print(string)

    def addPointOriginalLLH(self, LLH):
        self.originalLLH[self.idxOriginalLLH] = LLH
        self.idxOriginalLLH += 1

    def addPointOriginalNagelkerke(self, nagelkerke): # ASSUME THIS ALWAYS HAPPENS AFTER LLH
        self.originalNagelkerke[self.idxOriginalLLH-1] = nagelkerke

    def addPointOriginalMcFadden(self, mcfadden):
        self.originalMcFadden[self.idxOriginalLLH-1] = mcfadden

    def addPointTrainingLLH(self, LLH):
        self.trainingLLH[self.idxTrainingLLH] = LLH
        self.idxTrainingLLH += 1

    def addPointTrainingNagelkerke(self, nagelkerke): # ASSUME THIS ALWAYS HAPPENS AFTER LLH
        self.trainingNagelkerke[self.idxTrainingLLH-1] = nagelkerke

    def addPointTrainingMcFadden(self, mcfadden):
        self.trainingMcFadden[self.idxTrainingLLH-1] = mcfadden

    def addPointTestLLH(self, LLH):
        self.testLLH[self.idxTestLLH] = LLH
        self.idxTestLLH += 1

    def addPointTestNagelkerke(self, nagelkerke): # ASSUME THIS ALWAYS HAPPENS AFTER LLH
        self.testNagelkerke[self.idxTestLLH-1] = nagelkerke

    def addPointTestMcFadden(self, mcfadden):
        self.testMcFadden[self.idxTestLLH-1] = mcfadden

    def setParameters(self, paramDict):
        for k, v in paramDict.items():
            self.parameters[k] = v

    def getParameters(self):
        params = list()
        if self.model == "Logit":
            if not self.options.fixA.get():
                params.append(self.parameters['a'])
            if not self.options.fixB.get():
                params.append(self.parameters['b'])
        else:
            if not self.options.fixN.get():
                params.append(self.parameters['n'])
            if not self.options.fixM.get():
                params.append(self.parameters['m'])
            if not self.options.fixTD50.get():
                params.append(self.parameters['TD50'])

        return params

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

        self.originalLLH = np.trim_zeros(self.originalLLH)
        self.trainingLLH = np.trim_zeros(self.trainingLLH)
        self.testLLH = np.trim_zeros(self.testLLH)

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

    def printResults(self, log):
        pString = ", ".join([f"{k} = {v:.3f}" for k, v in self.parameters.items()])
        self.print2(f"The original parameters from the NTCP fit were {pString}.")

        correctionStatistics = {'mean': {p: np.mean(self.parameterSpace[p]) for p in self.parameters.keys()},
                                'median': {p: np.median(self.parameterSpace[p]) for p in self.parameters.keys()},
                                'none': {k: v for k, v in self.parameters.items()}}

        self.parameters = {k: 2 * v - correctionStatistics[self.bootstrapCorrectionMethod][k] for k, v in self.parameters.items()}

        # Fix any parameters outside physical limits

        pString = ", ".join([f"{k} = {v:.3f}" for k, v in self.parameters.items()])
        self.print2(f"Using the {self.bootstrapCorrectionMethod} pivot bias correction method, the corrected best fits were {pString}.")

        self.print2(f"The bias corrected bootstrapped confidence intervals were:")
        for p in self.parameters.keys():
            self.print2(f"{p}\t= [{self.CI[p][0]:.2f}, {self.CI[p][1]:.2f}]")

    def plotResults(self):
        nPlots = self.model == "LKB" and 6 or 6
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
            if p == "n" and not self.nIsLinear:
                low = np.min(self.parameterSpace[p])
                high = np.max(self.parameterSpace[p])
                bins = np.logspace(np.log10(low * 0.8), np.log10(high * 1.2), 50)
                axs[idx].hist(x=self.parameterSpace[p], bins=bins)
                axs[idx].set_xscale('log')
            else:
                axs[idx].hist(x=self.parameterSpace[p], bins=50)

            axs[idx].set_xlabel(f"{p} values")
            axs[idx].set_title(f"Parameter space for {p} ({self.cohort})")

            uncorrectedCI = [np.percentile(self.parameterSpace[p], perc) for perc in [self.lowerPercent, self.upperPercent]]

            # Uncorrected
            axs[idx].plot([uncorrectedCI[0]] * 2, axs[idx].get_ybound(), "g-",
                          [uncorrectedCI[1]] * 2, axs[idx].get_ybound(), "g-")

            axs[idx].plot([self.CI[p][0]] * 2, axs[idx].get_ybound(), 'k-',
                          [self.CI[p][1]] * 2, axs[idx].get_ybound(), 'k-')

            if self.bootstrapCorrectionMethod == "none":
                ci_95[p] = [np.percentile(self.parameterSpace[p], k) for k in [97.5, 2.5]]
                ci_68[p] = [np.percentile(self.parameterSpace[p], k) for k in [84, 16]]
            else:
                ci_95[p] = [2 * parameters[p] - np.percentile(self.parameterSpace[p], k) for k in [97.5, 2.5]]
                ci_68[p] = [2 * parameters[p] - np.percentile(self.parameterSpace[p], k) for k in [84, 16]]

            filter_95 *= (self.parameterSpace[p] >= ci_95[p][0]) & (self.parameterSpace[p] <= ci_95[p][1])
            filter_68 *= (self.parameterSpace[p] >= ci_68[p][0]) & (self.parameterSpace[p] <= ci_68[p][1])

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
                if secondPar == 'n' and not self.nIsLinear:
                    axs[idx].set_yscale('log')

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
        axs[idx].hist(x=self.originalLLH, bins=50)
        axs[idx].set_xlabel("Log Likelihood values")
        axs[idx].set_title(f"Log Likelihood for {self.cohort}")

    def writeToFile(self, patients=None):
        with open(f"Output/bootstrapParameterSpace.csv", "w") as out:
            if self.model == "LKB":
                out.write("cohort,n,m,TD50,Original Nagelkerke R2,Original McFadden R2,Test Nagelkerke R2,Test McFadden R2,Training Nagelkerke R2,Training McFadden R2,")
                for name, patient in list(patients.items()):
                    out.write(f"NTCP {name} [%],")
                for name, patient in list(patients.items()):
                    out.write(f"gEUD {name} [%],")
                out.write("\n")

                assert self.idx == self.idxTestLLH == self.idxTrainingLLH
                assert self.idxOriginalLLH == 1

                originalPseudoR2Nagelkerke = self.testNagelkerke[0]
                originalPseudoR2McFadden = self.testMcFadden[0]

                for idx in range(len(self.parameterSpace['n'])):
                    n = self.parameterSpace['n'][idx]
                    m = self.parameterSpace['m'][idx]
                    TD50 = self.parameterSpace['TD50'][idx]
                    
                    trainingPseudoR2Nagelkerke = self.trainingNagelkerke[idx]
                    trainingPseudoR2McFadden = self.trainingMcFadden[idx]
                    testPseudoR2Nagelkerke = self.testNagelkerke[idx]
                    testPseudoR2McFadden = self.testMcFadden[idx]

                    out.write(f"{self.cohort},{n:.4f},{m:.4f},{TD50:.3f},")
                    out.write(f"{originalPseudoR2Nagelkerke:.4f}, {originalPseudoR2McFadden:.4f},")
                    out.write(f"{testPseudoR2Nagelkerke:.4f}, {testPseudoR2McFadden:.4f},")
                    out.write(f"{trainingPseudoR2Nagelkerke:.4f}, {trainingPseudoR2McFadden:.4f},")

                    # Write out separate NTCP values + NagelKerke's R2
                    for name, patient in list(patients.items()):
                        gEUD = patient.getGEUD(n)
                        NTCP = HPM((gEUD - TD50) / (m * TD50))
                        NTCPpercent = NTCP * 100
                        out.write(f"{NTCPpercent:.3f},")
                    for name, patient in list(patients.items()):
                        gEUD = patient.getGEUD(n)
                        out.write(f"{gEUD:.2f},")
                    out.write("\n")

            else:
                out.write("cohort,a,b,TD50,Original Nagelkerke R2,Original McFadden R2,Test Nagelkerke R2,Test McFadden R2,Training Nagelkerke R2,Training McFadden R2\n")
                for idx in range(len(self.parameterSpace['a'])):
                    a = self.parameterSpace['a'][idx]
                    b = self.parameterSpace['b'][idx]
                    TD50 = self.parameterSpace['TD50'][idx]

                    trainingPseudoR2Nagelkerke = self.trainingNagelkerke[idx]
                    trainingPseudoR2McFadden = self.trainingMcFadden[idx]
                    testPseudoR2Nagelkerke = self.testNagelkerke[idx]
                    testPseudoR2McFadden = self.testMcFadden[idx]

                    out.write(f"{self.cohort},{a:.4f},{b:.4f},{TD50:3f},")
                    out.write(f"{originalPseudoR2Nagelkerke:.4f}, {originalPseudoR2McFadden:.4f},")
                    out.write(f"{testPseudoR2Nagelkerke:.4f}, {testPseudoR2McFadden:.4f},")
                    out.write(f"{trainingPseudoR2Nagelkerke:.4f}, {trainingPseudoR2McFadden:.4f}\n")