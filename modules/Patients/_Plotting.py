import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import *

from tkinter import *
from tkinter import ttk

from ..Tools import *


def drawSigmoid(self, log, style1, style2):
    """Plot the patient outcomes vs sigmoid probability from the best parameter set."""

    D5_lower = None
    D50_lower = None
    D5_upper = None
    D50_upper = None

    if self.options.NTCPcalculation.get() == "LKB":
        n = self.bestParameters[0]
        m = self.bestParameters[1]
        TD50 = self.bestParameters[2]

    else:
        a = self.bestParameters[0]
        b = self.bestParameters[1]

    x = np.arange(0, 120, 0.1)
    yExtra = []

    if self.options.NTCPcalculation.get() == "LKB":
        y = np.fromiter([HPM((k - TD50) / (m * TD50)) for k in x], float)
        if self.pSpace:
            for n_, m_, TD50_ in zip(self.pSpace.parameterSpace['n'], self.pSpace.parameterSpace['m'], self.pSpace.parameterSpace['TD50']):
                if (self.options.fixN.get() or self.pSpace.CI['n'][0] < n_ < self.pSpace.CI['n'][1]) and \
                    (self.options.fixM.get() or self.pSpace.CI['m'][0] < m_ < self.pSpace.CI['m'][1]) and \
                        (self.options.fixTD50.get() or self.pSpace.CI['TD50'][0] < TD50_ < self.pSpace.CI['TD50'][1]):
                    yExtra.append(np.fromiter([HPM((k - TD50_) / (m_ * TD50_)) for k in x], float))

            yMaxEmpirical = np.zeros(len(x))
            yMinEmpirical = np.zeros(len(x))

            for xIdx, xVal in enumerate(x):
                minY = 2
                maxY = -1
                for yIdx in range(len(yExtra)):
                    minY = min(minY, yExtra[yIdx][xIdx])
                    maxY = max(maxY, yExtra[yIdx][xIdx])
                yMaxEmpirical[xIdx] = maxY
                yMinEmpirical[xIdx] = minY

            D5_lower = x[np.argmin(abs(yMaxEmpirical - 0.05))]
            D5_upper = x[np.argmin(abs(yMinEmpirical - 0.05))]
            D50_lower = x[np.argmin(abs(yMaxEmpirical - 0.5))]
            D50_upper = x[np.argmin(abs(yMinEmpirical - 0.5))]

    else:
        y = np.fromiter([1 - 1 / (1 + exp(a + b * k)) for k in x], float)
        if self.pSpace:
            for a_, b_, TD50_ in zip(self.pSpace.parameterSpace['a'], self.pSpace.parameterSpace['b'], self.pSpace.parameterSpace['TD50']):
                if self.pSpace.CI['TD50'][0] < TD50_ < self.pSpace.CI['TD50'][1]:
                    yExtra.append(np.fromiter([1 - 1 / (1 + exp(a_ + b_ * k)) for k in x], float))

            yMaxEmpirical = np.zeros(len(x))
            yMinEmpirical = np.zeros(len(x))

            for xIdx, xVal in enumerate(x):
                minY = 2
                maxY = -1
                for yIdx in range(len(yExtra)):
                    minY = min(minY, yExtra[yIdx][xIdx])
                    maxY = max(maxY, yExtra[yIdx][xIdx])
                yMaxEmpirical[xIdx] = maxY
                yMinEmpirical[xIdx] = minY

            D5_lower = x[np.argmin(abs(yMaxEmpirical - 0.05))]
            D5_upper = x[np.argmin(abs(yMinEmpirical - 0.05))]
            D50_lower = x[np.argmin(abs(yMaxEmpirical - 0.5))]
            D50_upper = x[np.argmin(abs(yMinEmpirical - 0.5))]

    size = len(self.patients)

    px = np.zeros(size)
    py = np.zeros(size)
    ps = np.zeros(size)

    idx = 0

    for name, patient in self.patients.items():
        if self.options.NTCPcalculation.get() == "LKB":
            px[idx] = patient.getGEUD(n)
        else:
            px[idx] = patient.getDpercent()
        py[idx] = patient.getTox() >= self.options.toxLimit.get()
        nn = name.split("_")[0].split("tox")[0]
        if self.NTCPTimeDict:
            py[idx] *= np.exp(-(0.38 * self.NTCPTimeDict[nn] / 12)**1.37)
        idx += 1

    style1 = style1.pop(0)
    style2 = style2.pop(0)

    plt.plot(x, y, "-", color=style1, zorder=15, linewidth=3, label=f"{self.cohort}")
    plt.plot(px, py, "o", color=style2, zorder=0)

    if self.options.confidenceIntervalShowModels.get():
        for each in yExtra:
            plt.plot(x, each, '-', color="black", linewidth=0.25, zorder=0)

    if self.pSpace:
        plt.fill_between(x, yMinEmpirical, yMaxEmpirical, color=style2,
                         alpha=0.3, zorder=10)  # label=f"{self.options.confidenceIntervalPercent.get():.0f}% confidence interval")

    plt.ylim([-0.1, 1.1])
    Dlabel = self.options.useNTCPcc.get() and "cc" or "%"
    CIstr = self.pSpace and f" with {self.options.confidenceIntervalPercent.get()}% CI" or ""
    if self.options.NTCPcalculation.get() == "LKB":
        plt.xlabel("gEUD [Gy]")
        if D50_lower:
            plt.title(f"LKB{CIstr}; n = {n:.3f} ({self.pSpace.CI['n'][0]:.3f}-{self.pSpace.CI['n'][1]:.3f}), m = {m:.3f} ({self.pSpace.CI['m'][0]:.3f}-{self.pSpace.CI['m'][1]:.3f}), "
                      f"TD50 = {TD50:.2f} ({D50_lower:.2f}-{D50_upper:.2f}) Gy.")
        else:
            plt.title(f"LKB{CIstr}; n = {n:.3f}, m = {m:.3f}, TD50 = {TD50:.2f} Gy.")

    else:
        plt.xlabel(f"D{self.options.NTCPcalculationDpercent.get()}{Dlabel}[Gy]")
        if not D50_lower:
            plt.title(f"Logit{CIstr} using D{self.options.NTCPcalculationDpercent.get():.0f}{Dlabel}. "
                      f"TD5 = {self.calculateTDxFromLogit(5):.1f} Gy, TD50 = {self.calculateTDxFromLogit(50):.1f} Gy.")
        else:
            plt.title(f"Logit{CIstr} using D{self.options.NTCPcalculationDpercent.get():.0f}{Dlabel}. "
                      f"TD5 = {self.calculateTDxFromLogit(5):.1f} ({D5_lower:.1f}-{D5_upper:.1f}) Gy, "
                      f"TD50 = {self.calculateTDxFromLogit(50):.1f} ({D50_lower:.1f}-{D50_upper:.1f}) Gy.")
            log(f"TD5 = {self.calculateTDxFromLogit(5):.1f} ({D5_lower:.1f}-{D5_upper:.1f}) Gy "
                f"({self.options.confidenceIntervalPercent.get()}% CI)")
            log(f"TD50 = {self.calculateTDxFromLogit(50):.1f} ({D50_lower:.1f}-{D50_upper:.1f}) Gy "
                f"({self.options.confidenceIntervalPercent.get()}% CI)")
    plt.ylabel("Normal Tissue Complication Probability")
    plt.legend()

    # plt.savefig(f"Output/LKBgraph_{self.cohort}_{self.structure}.png")

    if self.options.NTCPcalculation.get() == "LKB":
        of = open(f"Output/LKBmodel_{self.cohort}_{self.structure}.csv", "w")
        of.write(f"n = {n:.3f}, m = {m:.3f}, TD50 = {TD50:.2f} Gy.\n\n")
        of.write("gEUD tox from patients\n")
        of.write("\n".join([" ".join([str(k[0]), str(k[1])]) for k in zip(px, py)]))
        of.write("\n\ngEUD LKB from model\n")
        of.write("\n".join([" ".join([str(k[0]), str(k[1])]) for k in zip(x, y)]))
        of.close()
    else:
        of = open(f"Output/LogitModel_{self.cohort}_{self.structure}.csv", "w")
        of.write(f"a = {a:.3f}, b = {b:.3f}.\n\n")
        of.write("gEUD tox from patients\n")
        of.write("\n".join([" ".join([str(k[0]), str(k[1])]) for k in zip(px, py)]))
        of.write("\n\ngEUD LKB from model\n")
        of.write("\n".join([" ".join([str(k[0]), str(k[1])]) for k in zip(x, y)]))
        of.close()


def drawAUROC(self, extraPatients, progress):
    """Plot the AUROC curve for the chosen patient cohorts."""

    def auroc_error(theta, nplus, nminus):
        """From Boule et al, Acta Oncologica 48 (2009)

        https://www.tandfonline.com/doi/full/10.1080/02841860903078513"""
        Q1 = (theta) / (2 - theta)
        Q2 = (2 * theta**2) / (1 + theta)
        nom = theta * (1 - theta) + (nplus - 1) * (Q1 - theta**2) + (nminus - 1) * (Q2 - theta**2)
        denom = nplus * nminus
        root = sqrt(nom / denom)
        return root

    size = (len(self.patients))
    for cohort in extraPatients:
        size += len(cohort.patients)

    gEUDlist = np.zeros(size)

    nList = np.arange(self.options.nFrom.get(), self.options.nTo.get(), 0.01)

    aurocX = np.zeros(len(nList))
    aurocY = np.zeros(len(nList))
    aurocError = np.zeros(len(nList))
    aurocIDX = 0
    progress['maximum'] = len(nList)

    for n in nList:
        progress.step(1)
        progress.update_idletasks()
        idx = 0
        for patient in list(self.patients.values()):
            gEUDlist[idx] = patient.getGEUD(n)
            idx += 1
        for cohort in extraPatients:
            for patient in list(cohort.patients.values()):
                gEUDlist[idx] = patient.getGEUD(n)
                idx += 1

        gEUDlist.sort()

        doseThresholds = [-1]
        doseThresholds += [(gEUDlist[k] + gEUDlist[k + 1]) / 2 for k in range(len(gEUDlist) - 1)]
        doseThresholds += [1000]
        y = np.zeros(len(doseThresholds))
        x = np.zeros(len(doseThresholds))

        for idx, dose in enumerate(doseThresholds):
            tp = fn = nP = nN = 0
            for patient in list(self.patients.values()):
                outcome = patient.getTox() >= self.options.toxLimit.get()
                prediction = patient.getGEUD(n) > dose
                tp += outcome and prediction
                fn += (not outcome) and prediction
                nP += outcome
                nN += not outcome

            for cohort in extraPatients:
                for patient in list(cohort.patients.values()):
                    outcome = patient.getTox() >= self.options.toxLimit.get()
                    prediction = patient.getGEUD(n) > dose
                    tp += outcome and prediction
                    fn += (not outcome) and prediction
                    nP += outcome
                    nN += not outcome

            if (nN):
                x[idx] = fn / nN
            if (nP):
                y[idx] = tp / nP

        area = np.sum([(x[k] - x[k + 1]) * y[k + 1] for k in range(len(doseThresholds) - 1)])
        aurocX[aurocIDX] = n
        aurocY[aurocIDX] = area
        aurocError[aurocIDX] = auroc_error(area, nP, nN)
        aurocIDX += 1

    maxval = np.amax(aurocY)
    maxind = np.where(aurocY == maxval)
    xlimits = [aurocX[maxind[0][0]], aurocX[maxind[0][-1]]]
    xmean = (xlimits[0] + xlimits[1]) / 2
    xwidth = (xlimits[1] - xlimits[0]) / 2

    plt.errorbar(aurocX, aurocY, xerr=None, yerr=aurocError)
    plt.plot([xmean, xmean], [0, maxval], "g--")
    plt.ylim(0, 1.1)
    plt.title(f"AUROC for {self.cohort}. Max is {maxval:.2f} at n={xmean:.2f} +- {xwidth:.2f}.")
    plt.ylabel('Area Under ROC')
    plt.xlabel('gEUD n-value')
    plt.show()
    progress['value'] = 0
