import random
import time
import numpy as np
import matplotlib.pyplot as plt
from math import *
from tkinter import *
from ..Patients import *


def calculateLKBuncert(self):
    """Based on http://stacks.iop.org/1742-6596/489/i=1/a=012087?key=crossref.181b59106e0d253de74e704220e16c36.

     Choose between profile likelihood method, parametric and non-parametric bootstrapping (time demanding!!!!!)
     Profile likelihood method: Vary each parameter individually until the LLH is decreased by an amount equal to half the
       the critical value of chi2(1) distribution at the desired significance level
       Critical value: http://www3.med.unipmn.it/~magnani/pdf/Tavole_chi-quadrato.pdf. 1.92?

     Non-parametric: Randomly choose patients, with replacement, from sample 1000-2000 times and calculate n,m,TD50 from sample
       The confidence interval is chosen from the distribution of n,m,TD50 values
       Assumptions: The cohort covers the different parameters, a robust method.
       My assumption: We choose the same number of patients.

     Parametric bootstrapping: A synthesised population is derived from the found n,m,TD50
       For each patient, calculate the NTCP using the optimized n,m,TD50 values.
       Generate a random number rn [0,1] for each patient; rn < NTCP -> tox, rn > NTCP -> no tox
       Find n,m,TD50 from the synthesized population using these values, repeat 1000-2000 times
       Assumption: The cohort describes the real population. """

    ## REDO THIS
    if not len(self.bestParameters):
        cohortList = list(self.patients.values())
        for cohort in cohortList:
            cohort.options = self.options
            cohort.results = Results(cohort.options)

        if self.options.NTCPcalculation.get() == "Logit":
            for patients in cohortList:
                patients.calculateDpercent(self.options.NTCPcalculationDpercent.get())
                patients.bestParameters = self.bestParameters

        primaryCohort = cohortList[0]
        secondaryCohorts = len(cohortList) > 1 and cohortList[1:] or {}

        if self.options.optimizationScheme.get() == "GradientDescent":
            res = primaryCohort.doGradientOptimization(secondaryCohorts, self.progress)

        elif self.options.optimizationScheme.get() == "MatrixMinimization":
            res = primaryCohort.doMatrixMinimization(secondaryCohorts, self.progress)

        self.bestParameters = res.x
        self.results = Results(self.options)
        self.results.setParameters(res.x)

    nIterations = self.options.confidenceIntervalIterations.get()
    origIt = self.options.basinHoppingIterations.get()
    self.options.basinHoppingIterations.set(2)
    self.progress['maximum'] = nIterations

    """
    nHist = np.zeros(nIterations)
    mHist = np.zeros(nIterations)
    TD50Hist = np.zeros(nIterations)
    LLHhist = np.zeros(nIterations)

    nInit = self.bestParameters[0]
    mInit = self.bestParameters[1]
    if self.options.NTCPcalculation.get() == "LKB":
        TD50Init = self.bestParameters[2]
    else:
        TD50Init = -nInit / mInit
    """

    grade = self.options.toxLimit.get()

    for patientsInCohort in self.patients.values():
        patientsInCohort.calculateNTCP()

    cohortList = list(self.patients.values())
    primaryCohort = cohortList[0]
    secondaryCohorts = len(cohortList) > 1 and cohortList[1:] or {}
    time1 = time.time()
    idx = 0
    if self.options.confidenceIntervalMethod.get() == "ParametricBootstrapping":
        for patientsInCohort in self.patients.values():
            patientsInCohort.saveTox()

        for k in range(nIterations):
            self.progress.step(1)
            self.progress.update_idletasks()

            for patientsInCohort in self.patients.values():
                for patient in patientsInCohort.patients.values():
                    rn = random.random()
                    ntcp = patient.getNTCP()
                    patient.setTox(rn < ntcp and grade or 0)

            if self.options.optimizationScheme.get() == "GradientDescent":
                res = primaryCohort.doGradientOptimization(secondaryCohorts, None)
            elif self.options.optimizationScheme.get() == "MatrixMinimization":
                res = primaryCohort.doMatrixMinimization(secondaryCohorts, None)

            if res.fun < -self.options.confidenceIntervalLikelihoodLimit.get():
                continue

            primaryCohort.results.addHistogram(res.x)
#            nHist[idx] = res.x[0]
#            mHist[idx] = res.x[1]
#            if self.options.NTCPcalculation.get() == "LKB":
#                TD50Hist[idx] = res.x[2]
            LLHhist[idx] = -res.fun
            idx += 1

            for patientsInCohort in list(self.patients.values()):
                patientsInCohort.restoreTox()

    elif self.options.confidenceIntervalMethod.get() == "NonParametricBootstrapping":
        patientZip = []
        for cohort, patientsInCohort in list(self.patients.items()):
            for patientName in list(patientsInCohort.patients.keys()):
                patientZip.append((cohort, patientName))
        nPatients = len(patientZip)

        for k in range(nIterations):
            self.progress.step(1)
            self.progress.update_idletasks()

            newPatientCohort = Patients(self.options)
            newPatientCohort.bestParameters = list(self.bestParameters)

            nTox = 0
            for n in range(nPatients):
                thisPatient = Patient(None)
                randomPatientID = random.randint(0, nPatients - 1)

                thisCohort = patientZip[randomPatientID][0]
                thisName = patientZip[randomPatientID][1]
                thisPatient.setTox(self.patients[thisCohort].patients[thisName].getTox())
                nTox += thisPatient.getTox() >= self.options.toxLimit.get()
                if self.options.NTCPcalculation.get() == "LKB":
                    thisPatient.nList = self.patients[thisCohort].patients[thisName].nList
                    thisPatient.GEUDlist = self.patients[thisCohort].patients[thisName].GEUDlist
                else:
                    thisPatient.Dpercent = self.patients[thisCohort].patients[thisName].Dpercent
                thisPatient.setID(thisName)
                while thisName in newPatientCohort.patients:
                    thisName += "_"
                newPatientCohort.patients[thisName] = thisPatient

            if nTox == 0:
                continue

            if self.options.optimizationScheme.get() == "GradientDescent":
                res = newPatientCohort.doGradientOptimization({}, None)
            elif self.options.optimizationScheme.get() == "MatrixMinimization":
                res = newPatientCohort.doMatrixMinimization({}, None)

            del newPatientCohort

            if res.fun < -self.options.confidenceIntervalLikelihoodLimit.get():
                continue

            nHist[idx] = res.x[0]
            mHist[idx] = res.x[1]
            if self.options.NTCPcalculation.get() == "LKB":
                TD50Hist[idx] = res.x[2]
            LLHhist[idx] = -res.fun
            idx += 1

    elif self.options.confidenceIntervalMethod.get() == "ProfileLikelihood":
        res = primaryCohort.profileLikelihood(secondaryCohorts)
        if self.options.NTCPcalculation.get() == "Logit":
            self.log(f"a = {self.bestParameters[0]:.3f} ({res[0][0]:.3f} - {res[0][1]:.3f})")
            self.log(f"b = {self.bestParameters[1]:.3f} ({res[1][0]:.3f} - {res[1][1]:.3f})")
        else:
            self.log(f"n = {self.bestParameters[0]:.2f} ({res[0][0]:.2f}-{res[0][1]:.2f})")
            self.log(f"m = {self.bestParameters[1]:.2f} ({res[1][0]:.2f}-{res[1][1]:.2f})")
            self.log(f"TD50 = {self.bestParameters[2]:.2f} ({res[2][0]:.2f}-{res[2][1]:.2f})")

        self.confidenceInterval = res
        return
    time2 = time.time()
    self.options.basinHoppingIterations.set(origIt)

    nHist = np.trim_zeros(nHist)
    mHist = np.trim_zeros(mHist)
    TD50Hist = np.trim_zeros(TD50Hist)
    LLHhist = np.trim_zeros(LLHhist)

    if self.options.NTCPcalculation.get() == "LKB":
        with open(f"Output/LKBuncert{self.options.confidenceIntervalMethod.get()}.csv", "w") as out:
            out.write("n,m,TD50\n")
            for k in range(len(nHist)):
                out.write(f"{nHist[k]},{mHist[k]},{TD50Hist[k]}\n")

        lowerPercent = (100 - self.options.confidenceIntervalPercent.get()) / 2
        upperPercent = 100 - lowerPercent        

        nMedian = np.median(nHist)
        mMedian = np.median(mHist)
        TD50Median = np.median(TD50Hist)

        nMean = np.mean(nHist)
        mMean = np.mean(mHist)
        TD50Mean = np.mean(TD50Hist)

        print(f"nMedian = {nMedian:.2f}, mMedian = {mMedian:.2f}, TD50Median = {TD50Median:.2f}")
        print(f"nMean = {nMean:.2f}, mMean = {mMean:.2f}, TD50Mean = {TD50Mean:.2f}")

        nCI = [np.percentile(nHist, lowerPercent), np.percentile(nHist, upperPercent)]
        mCI = [np.percentile(mHist, lowerPercent), np.percentile(mHist, upperPercent)]
        TD50CI = [np.percentile(TD50Hist, lowerPercent), np.percentile(TD50Hist, upperPercent)]

        # pivot theorem
        nCIPivot = [2*nInit - nCI[1], 2*nInit-nCI[0]]
        mCIPivot = [2*mInit - mCI[1], 2*mInit-mCI[0]]
        TD50CIPivot = [2*TD50Init- TD50CI[1], 2*TD50Init-TD50CI[0]]

        # Don't push the CI outside the preset limits...
        nLim = [np.min(nHist), np.max(nHist)]
        mLim = [np.min(mHist), np.max(mHist)]
        TD50Lim = [np.min(TD50Hist), np.max(TD50Hist)]

        if nCIPivot[0] < nLim[0]:
            nCIPivot = [k + nLim[0] - nCIPivot[0] for k in nCIPivot]
        if nCIPivot[1] > nLim[1]:
            nCIPivot = [k - (nCIPivot[1] - nLim[1]) for k in nCIPivot]
        if mCIPivot[0] < mLim[0]:
            mCIPivot = [k + mLim[0] - mCIPivot[0] for k in mCIPivot]
        if mCIPivot[1] > mLim[1]:
            mCIPivot = [k - (mCIPivot[1] - mLim[1]) for k in mCIPivot]
        if TD50CIPivot[0] < TD50Lim[0]:
            TD50CIPivot = [k + TD50Lim[0] - TD50CIPivot[0] for k in TD50CIPivot]
        if TD50CIPivot[1] > TD50Lim[1]:
            TD50CIPivot = [k - (TD50CIPivot[1] - TD50Lim[1]) for k in TD50CIPivot]

        self.log(f"\nFinished Confidence Interval tests ({(time2-time1)/60:.1f} minutes).")
        self.log(f"{self.options.confidenceIntervalPercent.get()}% CI calculated as "
                 f"the percentiles of a {self.options.confidenceIntervalMethod.get()} procedure.")
        self.log(f"The original parameters were n = {nInit:.2f}, m = {mInit:.2f} and TD50 = {TD50Init:.1f} Gy.")
        if self.options.bootstrapCorrectionMethod.get() == 'mean':
            self.log("Using the MEAN bias correction method, the corrected best fits are:")
            self.log(f"n = {2*nInit-nMean:.2f};\tm = {2*mInit-mMean:.2f};TD50 = {2*TD50Init-TD50Mean:.2f}\t")
        elif self.options.bootstrapCorrectionMethod.get() == 'median':
            self.log("Using the MEDIAN bias correction method, the corrected best fits are:")
            self.log(f"n = {2*nInit-nMedian:.2f};\tm = {2*mInit-mMedian:.2f};TD50 = {2*TD50Init-TD50Median:.2f}\t")
        if self.options.bootstrapCorrectionMethod.get() == 'none':
            self.log("Bootstrapped confidence intervals:")
            self.log(f"n\t= [{nCI[0]:.2f},  {nCI[1]:.2f}]")
            self.log(f"m\t= [{mCI[0]:.2f},  {mCI[1]:.2f}]")
            self.log(f"TD50\t= [{TD50CI[0]:.1f},  {TD50CI[1]:.1f}]")
        else:
            self.log("\"Pivot\" bootstrapped confidence intervals:")
            self.log(f"n\t= [{nCIPivot[0]:.2f},  {nCIPivot[1]:.2f}]")
            self.log(f"m\t= [{mCIPivot[0]:.2f},  {mCIPivot[1]:.2f}]")
            self.log(f"TD50\t= [{TD50CIPivot[0]:.1f},  {TD50CIPivot[1]:.1f}]")

        self.confidenceIntervalNone = [nCI, mCI, TD50CI]
        self.confidenceIntervalMean = [nCIPivot, mCIPivot, TD50CIPivot]
        self.confidenceIntervalMedian = [nCIPivot, mCIPivot, TD50CIPivot]
        self.bestParametersNone = [nInit, mInit, TD50Init]
        self.bestParametersMean = [2*nInit-nMean, 2*mInit-mMean, 2*TD50Init-TD50Mean]
        self.bestParametersMedian = [2*nInit-nMedian, 2*mInit-mMedian, 2*TD50Init-TD50Median]

        if self.options.bootstrapCorrectionMethod.get() == "none":
            self.bestParameters = self.bestParametersNone
            self.confidenceInterval = self.confidenceIntervalNone

        elif self.options.bootstrapCorrectionMethod.get() == "median":
            self.bestParameters = self.bestParametersMedian
            self.confidenceInterval = self.confidenceIntervalMedian

        elif self.options.bootstrapCorrectionMethod.get() == "mean":
            self.bestParameters = self.bestParametersMean
            self.confidenceInterval = self.confidenceIntervalMean

        self.aHist = nHist
        self.bHist = mHist
        self.TD50Hist = TD50Hist

        fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(1, 6, sharex=False, sharey=False, figsize=(30, 5))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        plt.ylabel("Frequency")
        ax1.hist(x=nHist, bins=50)
        ax1.set_xlabel("n value distribution")
        ax2.hist(x=mHist, bins=50)
        ax2.set_xlabel("m value distribution")
        ax3.hist(x=TD50Hist, bins=50)
        ax3.set_xlabel("TD50 value distribution")
        
        y11, y12 = ax1.get_ybound()
        y21, y22 = ax2.get_ybound()
        y31, y32 = ax3.get_ybound()
        
        if self.options.bootstrapCorrectionMethod.get() in ["mean", "median"]:
            ax1.plot([nCIPivot[0], nCIPivot[0]], [y11, y12], "k-")            
            ax1.plot([nCIPivot[1], nCIPivot[1]], [y11, y12], "k-")            
            ax2.plot([mCIPivot[0], mCIPivot[0]], [y21, y22], "k-")            
            ax2.plot([mCIPivot[1], mCIPivot[1]], [y21, y22], "k-")
            ax3.plot([TD50CIPivot[0], TD50CIPivot[0]], [y31, y32], "k-")
            ax3.plot([TD50CIPivot[1], TD50CIPivot[1]], [y31, y32], "k-")
        else:
            ax1.plot([nCI[0], nCI[0]], [y11, y12], "k-")            
            ax1.plot([nCI[1], nCI[1]], [y11, y12], "k-")            
            ax2.plot([mCI[0], mCI[0]], [y21, y22], "k-")            
            ax2.plot([mCI[1], mCI[1]], [y21, y22], "k-")
            ax3.plot([TD50CI[0], TD50CI[0]], [y31, y32], "k-")
            ax3.plot([TD50CI[1], TD50CI[1]], [y31, y32], "k-")

        n95 = [2*nInit - np.percentile(nHist, 97.5), 2*nInit - np.percentile(nHist,  2.5)]
        n68 = [2*nInit - np.percentile(nHist, 84.0), 2*nInit - np.percentile(nHist, 16.0)]
        m95 = [2*mInit - np.percentile(mHist, 97.5), 2*mInit - np.percentile(mHist,  2.5)]
        m68 = [2*mInit - np.percentile(mHist, 84.0), 2*mInit - np.percentile(mHist, 16.0)]
        TD5095 = [2*TD50Init - np.percentile(TD50Hist, 97.5), 2*TD50Init - np.percentile(TD50Hist,  2.5)]
        TD5068 = [2*TD50Init - np.percentile(TD50Hist, 84.0), 2*TD50Init - np.percentile(TD50Hist, 16.0)]
        
        if self.options.fixN.get():
            n95 = [self.options.nFrom.get(), self.options.nTo.get()]
            n68 = [self.options.nFrom.get(), self.options.nTo.get()]
        if self.options.fixM.get():
            m95 = [self.options.mFrom.get(), self.options.mTo.get()]
            m68 = [self.options.mFrom.get(), self.options.mTo.get()]
        if self.options.fixTD50.get():
            TD5095 = [self.options.TD50From.get(), self.options.TD50To.get()]
            TD5068 = [self.options.TD50From.get(), self.options.TD50To.get()]

        if not self.options.fixN.get():
            nHist95 = nHist[(n95[0] < nHist) & (nHist < n95[1]) & (m95[0] < mHist) & (mHist < m95[1]) & (TD5095[0] < TD50Hist) & (TD50Hist < TD5095[1])]
            nHist68 = nHist[(n68[0] < nHist) & (nHist < n68[1]) & (m68[0] < mHist) & (mHist < m68[1]) & (TD5068[0] < TD50Hist) & (TD50Hist < TD5068[1])]
        else:
            nHist95 = nHist
            nHist68 = nHist
            
        if not self.options.fixM.get():
            mHist95 = mHist[(n95[0] < nHist) & (nHist < n95[1]) & (m95[0] < mHist) & (mHist < m95[1]) & (TD5095[0] < TD50Hist) & (TD50Hist < TD5095[1])]
            mHist68 = mHist[(n68[0] < nHist) & (nHist < n68[1]) & (m68[0] < mHist) & (mHist < m68[1]) & (TD5068[0] < TD50Hist) & (TD50Hist < TD5068[1])]
        else:
            mHist95 = mHist
            mHist68 = mHist
        
        if not self.options.fixTD50.get():
            TD50Hist95 = TD50Hist[(n95[0] < nHist) & (nHist < n95[1]) & (m95[0] < mHist) & (mHist < m95[1]) & (TD5095[0] < TD50Hist) & (TD50Hist < TD5095[1])]
            TD50Hist68 = TD50Hist[(n68[0] < nHist) & (nHist < n68[1]) & (m68[0] < mHist) & (mHist < m68[1]) & (TD5068[0] < TD50Hist) & (TD50Hist < TD5068[1])]
        else:
            TD50Hist95 = TD50Hist
            TD50Hist68 = TD50Hist

        ax4.plot(TD50Hist, nHist, "o", c=self.options.fixN.get() and "green" or "red", zorder=0, label="Outside 95% LLH")
        if not self.options.fixN.get() and not self.options.fixTD50.get():
            ax4.plot(TD50Hist95, nHist95, "o", c="yellow", zorder=5, label="Within 95% LLH")
            ax4.plot(TD50Hist68, nHist68, "o", c="green", zorder=10, label="Within 68% LLH")
        ax4.legend()
        ax4.set_xlabel("TD50 parameter")
        ax4.set_ylabel("n parameter")
        ax5.plot(TD50Hist, mHist, "o", c=self.options.fixM.get() and "green" or "red", zorder=0, label="Outside 95% LLH")
        if not self.options.fixM.get() and not self.options.fixTD50.get():
            ax5.plot(TD50Hist95, mHist95, "o", c="yellow", zorder=5, label="Within 95% LLH")
            ax5.plot(TD50Hist68, mHist68, "o", c="green", zorder=10, label="Within 68% LLH")
        ax5.legend()
        ax5.set_xlabel("TD50 parameter")
        ax5.set_ylabel("m parameter")
        ax6.hist(x=LLHhist, bins=50)
        ax6.set_xlabel("Log Likelihood values")            
        plt.show()
        
    elif self.options.NTCPcalculation.get() == "Logit":
        lowerPercent = (100 - self.options.confidenceIntervalPercent.get()) / 2
        upperPercent = 100 - lowerPercent

        # 50% = 1 - 1 / exp(a + b TD50)
        # exp(a + b TD50) = 1
        # a + b TD50 = 0
        # TD50 = - a / b

        # LIKEWISE
        # 5% = 1 - 1 / exp(a + b TD5)
        # exp(a + b TD5) = 19
        # a + b TD5 = ln(19)
        # TD5 = (ln(19) - a) / b

        TD50Hist = [ -a/b for a,b in zip(nHist, mHist) ]

        with open(f"Output/LKBuncert{self.options.confidenceIntervalMethod.get()}.csv", "w") as out:
            out.write("a,b,TD50\n")
            for k in range(len(nHist)):
                out.write(f"{nHist[k]},{mHist[k]},{TD50Hist[k]}\n")
        
        if self.options.fixA.get():
            nCI = [nInit, nInit]
        else:
            nCI = [np.percentile(nHist, lowerPercent), np.percentile(nHist, upperPercent)]
        if self.options.fixB.get():
            mCI = [mInit, mInit]
        else:
            mCI = [np.percentile(mHist, lowerPercent), np.percentile(mHist, upperPercent)]

        TD50CI = [np.percentile(TD50Hist, lowerPercent), np.percentile(TD50Hist, upperPercent)]
        
        nMedian = np.median(nHist)
        mMedian = np.median(mHist)
        TD50Median = np.median(TD50Hist)
        nMean = np.mean(nHist)
        mMean = np.mean(mHist)
        TD50Mean = np.mean(TD50Hist)
        
        nCIPivot = [2*nInit - nCI[1], 2*nInit-nCI[0]]
        mCIPivot = [2*mInit - mCI[1], 2*mInit-mCI[0]]
        TD50CIPivot = [2*TD50Init - TD50CI[1], 2*TD50Init-TD50CI[0]]

        # Don't push the CI outside the preset limits...
        nLim = [np.min(nHist), np.max(nHist)]
        mLim = [np.min(mHist), np.max(mHist)]
        TD50Lim = [np.min(TD50Hist), np.max(TD50Hist)]
        
        if nCIPivot[0] < nLim[0]:
            nCIPivot = [k + nLim[0] - nCIPivot[0] for k in nCIPivot]
        if nCIPivot[1] > nLim[1]:
            nCIPivot = [k - (nCIPivot[1] - nLim[1]) for k in nCIPivot]
        if mCIPivot[0] < mLim[0]:
            mCIPivot = [k + mLim[0] - mCIPivot[0] for k in mCIPivot]
        if mCIPivot[1] > mLim[1]:
            mCIPivot = [k - (mCIPivot[1] - mLim[1]) for k in mCIPivot]
        if TD50CIPivot[0] < TD50Lim[0]:
            TD50CIPivot = [k + TD50Lim[0] - TD50CIPivot[0] for k in TD50CIPivot]
        if TD50CIPivot[1] > TD50Lim[1]:
            TD50CIPivot = [k - (TD50CIPivot[1] - TD50Lim[1]) for k in TD50CIPivot]
        
        self.log(f"\nFinished Confidence Interval tests ({(time2-time1)/60:.1f} minutes).")
        self.log(f"{self.options.confidenceIntervalPercent.get()}% CI calculated as "
                 f"the percentiles of a {self.options.confidenceIntervalMethod.get()} procedure.")
        self.log(f"The original parameters were a = {nInit:.2f} and b = {mInit:.2f} (-> TD50 = {TD50Init:.2f} Gy).")
        if self.options.bootstrapCorrectionMethod.get() == 'mean':
            self.log("Using the MEAN bias correction method, the corrected best fits are:")
            self.log(f"a = {2*nInit-nMean:.2f};\tb = {2*mInit-mMean:.2f}\tTD50 = {2*TD50Init-TD50Mean:.2f} Gy")
        elif self.options.bootstrapCorrectionMethod.get() == 'median':
            self.log("Using the MEDIAN bias correction method, the corrected best fits are:")
            self.log(f"a = {2*nInit-nMedian:.2f};\tb = {2*mInit-mMedian:.2f}\tTD50 = {2*TD50Init-TD50Median:.2f} Gy")
        if self.options.bootstrapCorrectionMethod.get() == 'none':
            self.log("The values of the bootstrapped distributions are shown below as median (CI):")
            self.log(f"a\t= [{nCI[0]:.2f},  {nCI[1]:.2f}]")
            self.log(f"b\t= [{mCI[0]:.2f},  {mCI[1]:.2f}]")
            self.log(f"TD50\t= {TD50CI[0]:.2f} - {TD50CI[1]:.2f}.")
        else:
            self.log("If the \"pivot\" method is applied, the following CIs are obtained:")
            self.log(f"a\t= {nCIPivot[0]:.2f} - {nCIPivot[1]:.2f}")
            self.log(f"b\t= {mCIPivot[0]:.2f} - {mCIPivot[1]:.2f}")
            self.log(f"TD50\t: {TD50CIPivot[0]:.2f} - {TD50CIPivot[1]:.2f}.")
          
        fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(1, 6, sharex=False, sharey=False, figsize=(20,8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        plt.ylabel("Frequency")
        ax1.hist(x=nHist, bins=50)
        ax1.set_xlabel("a value distribution")
        ax2.hist(x=mHist, bins=50)
        ax2.set_xlabel("b value distribution")

        n95 = [2*nInit - np.percentile(nHist, 97.5), 2*nInit - np.percentile(nHist,  2.5)]
        n68 = [2*nInit - np.percentile(nHist, 84.0), 2*nInit - np.percentile(nHist, 16.0)]
        m95 = [2*mInit - np.percentile(mHist, 97.5), 2*mInit - np.percentile(mHist,  2.5)]
        m68 = [2*mInit - np.percentile(mHist, 84.0), 2*mInit - np.percentile(mHist, 16.0)]

        if not self.options.fixA.get():
            nHist95 = nHist[(n95[0] < nHist) & (nHist < n95[1]) & (m95[0] < mHist) & (mHist < m95[1])]
            nHist68 = nHist[(n68[0] < nHist) & (nHist < n68[1]) & (m68[0] < mHist) & (mHist < m68[1])]
        else:
            nHist95 = nHist
            nHist68 = nHist
        
        if not self.options.fixB.get():
            mHist95 = mHist[(n95[0] < nHist) & (nHist < n95[1]) & (m95[0] < mHist) & (mHist < m95[1])]
            mHist68 = mHist[(n68[0] < nHist) & (nHist < n68[1]) & (m68[0] < mHist) & (mHist < m68[1])]
        else:
            mHist95 = mHist
            mHist68 = mHist

        ax3.plot(nHist, mHist, "o", c="red", zorder=0, label="Outside 95% CI")
        ax3.plot(nHist95, mHist95, "o", c="yellow", zorder=5, label="Within 95% CI")
        ax3.plot(nHist68, mHist68, "o", c="green", zorder=10, label="Within 68% CI")
        b, m = np.polynomial.polynomial.polyfit(nHist, mHist, 1)
        ax3.plot(nHist, b + m * nHist, "-")
        ax3.set_xlabel("a values")
        ax3.set_ylabel("b values")
        ax3.legend()
        
        ax4.hist(x=LLHhist, bins=50)
        ax4.set_xlabel("Log Likelihood values")

        ax5.hist(x=TD50Hist, bins=50)
        ax5.set_xlabel("TD50 calculation distribution")

        oldBestParametersNone = self.bestParametersNone
        oldBestParametersMean = self.bestParametersMean
        oldBestParametersMedian = self.bestParametersMedian
        
        self.confidenceIntervalNone = [nCI, mCI, TD50CI]
        self.confidenceIntervalMean = [nCIPivot, mCIPivot, TD50CIPivot]
        self.confidenceIntervalMedian = [nCIPivot, mCIPivot, TD50CIPivot]
        self.bestParametersNone = [nInit, mInit, TD50Init]
        self.bestParametersMean = [2*nInit-nMean, 2*mInit-mMean, 2*TD50Init - TD50Mean]
        self.bestParametersMedian = [2*nInit-nMedian, 2*mInit-mMedian, 2*TD50Init - TD50Median]

        if self.options.bootstrapCorrectionMethod.get() == "none":
            self.bestParameters = self.bestParametersNone
            self.confidenceInterval = self.confidenceIntervalNone

        elif self.options.bootstrapCorrectionMethod.get() == "median":
            self.bestParameters = self.bestParametersMedian
            self.confidenceInterval = self.confidenceIntervalMedian

        elif self.options.bootstrapCorrectionMethod.get() == "mean":
            self.bestParameters = self.bestParametersMean
            self.confidenceInterval = self.confidenceIntervalMean

        self.aHist = nHist
        self.bHist = mHist
        self.LLHhist = LLHhist
        self.correlationValueLogit = m
        self.log("Found correlation %.4f between 'a' and 'b' parameters." % (m))
        
        y11, y12 = ax1.get_ybound()
        y21, y22 = ax2.get_ybound()
        y51, y52 = ax5.get_ybound()
        
        if self.options.bootstrapCorrectionMethod.get() in ["mean", "median"]:
            ax1.plot([nCIPivot[0], nCIPivot[0]], [y11, y12], "k-")            
            ax1.plot([nCIPivot[1], nCIPivot[1]], [y11, y12], "k-")            
            ax2.plot([mCIPivot[0], mCIPivot[0]], [y21, y22], "k-")            
            ax2.plot([mCIPivot[1], mCIPivot[1]], [y21, y22], "k-")
            ax5.plot([TD50CIPivot[0], TD50CIPivot[0]], [y51, y52], "k-")
            ax5.plot([TD50CIPivot[1], TD50CIPivot[1]], [y51, y52], "k-")
        else:
            ax1.plot([nCI[0], nCI[0]], [y11, y12], "k-") 
            ax1.plot([nCI[1], nCI[1]], [y11, y12], "k-")            
            ax2.plot([mCI[0], mCI[0]], [y21, y22], "k-")            
            ax2.plot([mCI[1], mCI[1]], [y21, y22], "k-")
            ax5.plot([TD50CI[0], TD50CI[0]], [y51, y52], "k-")
            ax5.plot([TD50CI[1], TD50CI[1]], [y51, y52], "k-")

        # NOW check the distribution A - B
        if len(self.TD50Hist) and self.options.makeDifferentialCI.get():
            self.log("Found prior TD50 distribution, making differential histogram (prior - this):")
            TD50DiffHist = [TD50A - TD50B for TD50A, TD50B in zip(self.TD50Hist, TD50Hist)]
            ax6.hist(TD50DiffHist, bins=50)
            ax6.set_xlabel("Differential TD50 distribution")
            diffMean = np.mean(TD50DiffHist)
            diffMedian = np.mean(TD50DiffHist)                

            oldInit = oldBestParametersNone[2]
            thisInit = self.bestParametersNone[2]
            diffInit = oldInit - thisInit
            
            print(f"Prior cohort TD50 = {oldInit:.2f} Gy, this TD50 = {thisInit:.2f} Gy")

            if self.options.bootstrapCorrectionMethod.get() == "mean":
                bestParameterDiff = 2 * diffInit - diffMean
                diffCI = [2 * diffInit - np.percentile(TD50DiffHist, upperPercent),
                          2 * diffInit - np.percentile(TD50DiffHist, lowerPercent)]

            elif self.options.bootstrapCorrectionMethod.get() == "median":
                bestParameterDiff = 2 * diffInit - diffMedian
                diffCI = [2 * diffInit - np.percentile(TD50DiffHist, upperPercent),
                          2 * diffInit - np.percentile(TD50DiffHist, lowerPercent)]

            else:
                bestParameterDiff = diffInit
                diffCI = [np.percentile(TD50DiffHist, lowerPercent), np.percentile(TD50DiffHist, upperPercent)]

            y61, y62 = ax6.get_ybound()
            ax6.plot([diffCI[0], diffCI[0]], [y61, y62], "k-")
            ax6.plot([diffCI[1], diffCI[1]], [y61, y62], "k-")

            self.log(f"The differential TD50 distribution (A - B) is {bestParameterDiff:.2f} ({diffCI[0]:.2f} - {diffCI[1]:.2f})")
            
        self.TD50Hist = TD50Hist
        plt.show()

    # Write to output files
    bootstrapOutput = open("Output/bs.csv")
    for idx in range(len(self.mHist)):
        bootstrapOutput.write(f"{mHist[idx]:.4f},{nHist[idx]:.4f},{TD50Hist[idx]:.4f}\n")
    bootstrapOutput.close()