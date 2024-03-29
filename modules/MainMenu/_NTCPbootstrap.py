import random
import time
import numpy as np
import matplotlib.pyplot as plt
import math
from tkinter import *
from ..Patients import *

def calculatePseudoR2(patients, res, options):
    mean_y = model_n = LL0 = LLH = 0
    for name, patient in patients.items():
        tox = patient.getTox() >= options.toxLimit.get()
        mean_y += tox
        model_n += 1
    mean_y /= model_n

    eps = 1e-10

    for name, patient in patients.items():
        if options.NTCPcalculation.get() == "LKB":
            n, m, TD50 = res.x
            NTCP = HPM((patient.getGEUD(n) - TD50) / (m * TD50))
        else:
            a, b = res.x
            NTCP = 1 - 1 / (1 + math.exp(a + b * patient.getDpercent()))

        tox = patient.getTox() >= options.toxLimit.get()
        LL0 += tox * math.log(mean_y+eps) + (1-tox) * math.log(1 - mean_y + eps)
        LLH += tox * math.log(NTCP + eps) + (1-tox) * math.log(1 - NTCP + eps)

    try:
        nagelkerke = (1 - math.exp(-2*(LLH-LL0)/model_n)) / (1 - math.exp(2*LL0/model_n))
    except OverflowError as e:
        print("Overflow Error: ", e)
        print("LLH", LLH, "LL0", LL0, "model_n", model_n)
        nagelkerke = 0

    mcfadden = 1 - LLH/LL0

    return {'Nagelkerke': nagelkerke, 'McFadden': mcfadden, 'LLH': LLH}

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

    # REDO THIS

    self.window.destroy()  # Destroy bootstrap dialog window

    cohortList = list(self.patients.values())
    cohortNames = ", ".join(list(self.patients.keys()))
    nIterations = self.options.confidenceIntervalIterations.get()

    for patients in cohortList:
        patients.resetIdx()
        patients.options = self.options
        patients.pSpace = patients.ParameterSpace(nIterations, self.log, patients.cohort, self.options, patients.idx)

        if self.options.NTCPcalculation.get() == "Logit":
            patients.calculateDpercent(self.options.NTCPcalculationDpercent.get())
            if np.sum(patients.bestParameters):
                a = self.options.fixA.get() and self.options.aFrom.get() or patients.bestParameters[patients.idx['a']]
                b = self.options.fixB.get() and self.options.bFrom.get() or patients.bestParameters[patients.idx['b']]
                patients.pSpace.setParameters({'a': a, 'b': b})

        else:
            if np.sum(patients.bestParameters):
                n = self.options.fixN.get() and self.options.nFrom.get() or patients.bestParameters[patients.idx['n']]
                m = self.options.fixM.get() and self.options.mFrom.get() or patients.bestParameters[patients.idx['m']]
                TD50 = self.options.fixTD50.get() and self.options.TD50From.get() or patients.bestParameters[patients.idx['TD50']]

                patients.pSpace.setParameters({'n': n, 'm': m, 'TD50': TD50})

        if self.options.optimizationScheme.get() == "GradientDescent":
            res = patients.doGradientOptimization(self.progress)
            assert not isinstance(res, int)

        elif self.options.optimizationScheme.get() == "MatrixMinimization":
            res = patients.doMatrixMinimization(self.progress)

        patients.calculateNTCP()
        patients.bestParameters = res.x

        originalR2 = calculatePseudoR2(patients.patients, res, self.options)
        patients.pSpace.addPointOriginalLLH(originalR2['LLH'])
        patients.pSpace.addPointOriginalNagelkerke(originalR2['Nagelkerke'])
        patients.pSpace.addPointOriginalMcFadden(originalR2['McFadden'])

        print("Original McFadden is", originalR2['McFadden'])

    origIt = self.options.basinHoppingIterations.get()
    self.options.basinHoppingIterations.set(2)
    self.progress['maximum'] = nIterations * len(cohortList)

    grade = self.options.toxLimit.get()

    # Loop over cohorts early
    for cohort, patients in self.patients.items():
        print(f"Looping over cohort {cohort} of {cohortNames}")
        time1 = time.time()

        if self.options.confidenceIntervalMethod.get() == "ParametricBootstrapping":
            patients.saveTox()

            for k in range(nIterations):
                self.progress.step(1)
                self.progress.update_idletasks()

                for patient in patients.patients.values():
                    rn = random.random()
                    ntcp = patient.getNTCP()
                    patient.setTox(rn < ntcp and grade or 0)

                if self.options.optimizationScheme.get() == "GradientDescent":
                    res = patients.doGradientOptimization(None)
                elif self.options.optimizationScheme.get() == "MatrixMinimization":
                    res = patients.doMatrixMinimization(None)

                if res.fun < -self.options.confidenceIntervalLikelihoodLimit.get():
                    continue

                patients.pSpace.addPoint(res.x)
                
                trainingR2 = calculatePseudoR2(patients.patients, res, self.options)
                patients.pSpace.addPointTrainingLLH(trainingR2['LLH'])
                patients.pSpace.addPointTrainingNagelkerke(trainingR2['Nagelkerke'])
                patients.pSpace.addPointTrainingMcFadden(trainingR2['McFadden'])

                patients.restoreTox()

                testR2 = calculatePseudoR2(patients.patients, res, self.options)
                patients.pSpace.addPointTestLLH(testR2['LLH'])
                patients.pSpace.addPointTestNagelkerke(testR2['Nagelkerke'])
                patients.pSpace.addPointTestMcFadden(testR2['McFadden'])

        elif self.options.confidenceIntervalMethod.get() == "NonParametricBootstrapping":
            patientZip = list()
            for patientName in patients.patients.keys():
                patientZip.append((cohort, patientName))
            nPatients = len(patientZip)

            ntcp_file = open("Output/bootstrapTestNTCP.csv", "w")
            name_file = open("Output/bootstrapTestNames.csv", "w")

            for k in range(nIterations):
                print(".", end="")
                self.progress.step(1)
                self.progress.update_idletasks()

                newPatientCohort = Patients(self.options)
                newPatientCohort.bestParameters = list(patients.bestParameters)
                newPatientCohort.idx = patients.idx

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
                    res = newPatientCohort.doGradientOptimization(None)
                elif self.options.optimizationScheme.get() == "MatrixMinimization":
                    res = newPatientCohort.doMatrixMinimization(None)

                trainingR2 = calculatePseudoR2(newPatientCohort.patients, res, self.options)
                patients.pSpace.addPointTrainingLLH(trainingR2['LLH'])
                patients.pSpace.addPointTrainingNagelkerke(trainingR2['Nagelkerke'])
                patients.pSpace.addPointTrainingMcFadden(trainingR2['McFadden'])

                # Write the NTCP values of all patients to new file
                for name, patient in newPatientCohort.patients.items():
                    NTCP = patient.getNTCP()
                    ntcp_file.write(f"{NTCP:.3f},")
                    name_file.write(f"{name},")
                ntcp_file.write("\n")
                name_file.write("\n")

                del newPatientCohort

                patients.pSpace.addPoint(res.x)

                testR2 = calculatePseudoR2(patients.patients, res, self.options)
                patients.pSpace.addPointTestLLH(testR2['LLH'])
                patients.pSpace.addPointTestNagelkerke(testR2['Nagelkerke'])
                patients.pSpace.addPointTestMcFadden(testR2['McFadden'])

                if res.fun < -self.options.confidenceIntervalLikelihoodLimit.get():
                    continue
            
            ntcp_file.close()
            name_file.close()

        elif self.options.confidenceIntervalMethod.get() == "ProfileLikelihood":
            res = patients.profileLikelihood()
            if self.options.NTCPcalculation.get() == "Logit":
                patients.pSpace.CI["a"] = res[patients.idx['a']]
                patients.pSpace.CI["b"] = res[patients.idx['b']]
            else:
                patients.pSpace.CI["n"] = res[patients.idx['n']]
                patients.pSpace.CI["m"] = res[patients.idx['m']]
                patients.pSpace.CI["TD50"] = res[patients.idx['TD50']]

            patients.pSpace.printCI()

            patients.confidenceInterval = res
            return

        print("Done\n\n")

        time2 = time.time()
        self.options.basinHoppingIterations.set(origIt)

        ###################################################################
        # Finished performing the bootstrap, analyse the results ... :)   #
        ###################################################################

        patients.pSpace.trim()

        patients.pSpace.setPercentile(self.options.confidenceIntervalPercent.get())
        patients.pSpace.setBootstrapCorrectionMethod(self.options.bootstrapCorrectionMethod.get())
        patients.pSpace.calculateCI()
        patients.pSpace.applyPivot()
        patients.pSpace.writeToFile(patients.patients)

        self.log(f"\nFinished Confidence Interval tests for cohort {cohort} ({(time2-time1)/60:.1f} minutes).")
        self.log(f"{self.options.confidenceIntervalPercent.get()}% CI calculated as "
                 f"the percentiles of a {self.options.confidenceIntervalMethod.get()} procedure.")

        patients.pSpace.printResults(self.log)
        patients.bestParameters = patients.pSpace.getParameters()
        patients.pSpace.plotResults()
        self.options.didBootstrap = True

    plt.show()
