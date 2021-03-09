import os
import re
import pandas as pd
import numpy as np

from tkinter import *

from ..Tools import *
from .Patient import *


class Patients:
    def __init__(self, options):
        self.patients = {}
        self.bestParameters = []
        self.dataFolder = None
        self.structure = None
        self.cohort = None
        self.NTCPTimeDict = None
        self.doseUnit = None
        self.options = options

        """
        self.bestParametersNone = []
        self.bestParametersMean = []
        self.bestParametersMedian = []
        self.confidenceInterval = [[0, 0], [0, 0], [0, 0]]
        self.confidenceIntervalNone = [[0, 0], [0, 0], [0, 0]]
        self.confidenceIntervalMean = [[0, 0], [0, 0], [0, 0]]
        self.confidenceIntervalMedian = [[0, 0], [0, 0], [0, 0]]
        self.correlationLogit = -0.016
        self.calculateMeanDose = IntVar(value=1)
        self.aHist = []
        self.bHist = []
        self.TD50Hist = []
        self.LLHhist = []
        """

        self.pSpace = None

        self.plot = None
        self.ax1 = None

    def __repr__(self):
        return {'Cohort name': self.cohort, 'NPatients': len(patients)}

    def __len__(self):
        return len(patients)

    from ._Optimization import doGradientOptimization, doMatrixMinimization, profileLikelihood
    from ._Plotting import drawSigmoid, drawAUROC
    from ._ParameterSpace import ParameterSpace

    def setDataFolder(self, folder):
        self.dataFolder = folder

    def getFilePath(self, fileName):
        return self.dataFolder + "/" + fileName

    def findStructures(self):
        structureNames = []
        for root, d, f, in os.walk(self.dataFolder):
            for filename in f:
                # if not ".txt" in filename[-4:] and not ".dvh" in filename[-4]: continue
                if not filename[-4:] in [".txt", ".dvh"]:
                    continue
                if "gEUD" in filename:
                    continue
                with open(self.getFilePath(filename), "r") as textin:
                    for line in textin:
                        if self.options.DVHFileType.get() == "ECLIPSE":
                            if "Structure:" in line:
                                structureNames.append(re.sub("\s+", "", line.split(": ")[-1]))
                        elif self.options.DVHFileType.get() == "RayStation":
                            if "#RoiName:" in line:
                                structureNames.append(re.sub("\s+", "", line.split(":")[-1]))

        return list(set(structureNames))

    def findPlans(self):
        planNames = []
        for root, d, f, in os.walk(self.dataFolder):
            for filename in f:
                if not filename[-4:] in [".txt", ".dvh"]:
                    continue
                if "gEUD" in filename:
                    continue
                try:
                    with open(self.getFilePath(filename), "r") as textin:
                        for line in textin:
                            if self.options.DVHFileType.get() == "ECLIPSE":
                                if "Plan:" in line:
                                    planNames.append(re.sub("\s+", "", line.split(": ")[-1]))
                            elif self.options.DVHFileType.get() == "RayStation":
                                if "#Dosename:" in line:
                                    planNames.append(re.sub("\s+", "", line.split(":")[-1]))
                except Exception as e:
                    print(f"Skipping filename {filename} because {e}")
                    continue
        return list(set(planNames))

    def findNPatients(self):
        nPatients = 0
        for root, d, f in os.walk(self.dataFolder):
            for filename in f:
                if not ".txt" in filename[-4:] and not ".dvh" in filename[-4]:
                    continue
                if "gEUD" in filename:
                    continue
                nPatients += 1
        return nPatients

    def loadPatientsECLIPSE(self, progress):
        def match(a, b):
            a += "$"  # Don't match end-of-lines
            a = a.replace("*", ".*")  # Use regex type wildcard
            return re.match(a, b)

        log = []

        if not self.options.loadToxFromFilename.get():
            try:
                with open("%s_tox.csv" % (self.dataFolder), "r") as toxFile:
                    toxDict = {}
                    for line in toxFile.readlines():
                        linesplit = line.split(",")
                        toxDict[linesplit[0]] = int(linesplit[1])
            except:
                log.append(f"Cannot open toxicity file {self.dataFolder}_tox.csv, please include in format (filename,tox).")
                return log

        n = 0

        progress['maximum'] = self.findNPatients()

        printedLogOutput = False

        for root, d, f, in os.walk(self.dataFolder):
            for filename in f:
                if not ".txt" in filename[-4:]:
                    continue
                if "gEUD" in filename:
                    continue
                idx = 0
                structureNames = []
                structureStarts = []
                structureLength = []
                structureMeanDose = []
                structureMinDose = []
                structureMaxDose = []
                structureVolume = []
                planNames = []
                planStarts = []
                planLength = []
                lastLineEmpty = False

                progress.step(1)
                progress.update_idletasks()

                emptyStructure = False
                dontAppendLength = False
                doseHeaderLine = None

                with open(self.getFilePath(filename), "r") as textin:
                    for line in textin:
                        if "Plan: " in line and lastLineEmpty:
                            # Remove summary lines in start of file
                            # They always start with \nPlan, instead of \nStructure
                            idx += 1
                            continue

                        if "Structure: " in line:
                            structureNames.append(re.sub("\s+", "", line.split(": ")[-1]))
                            if structureStarts:
                                if dontAppendLength:
                                    dontAppendLength = False
                                else:
                                    structureLength.append(idx - structureStarts[-1] - 1)

                        if "Plan: " in line:
                            planNames.append(re.sub("\s+", "", line.split(": ")[-1]))
                            if planStarts:
                                if dontAppendLength:
                                    dontAppendLength = False
                                else:
                                    planLength.append(idx - planStarts[-1] - 3)

                        if line == "Min Dose [Gy]: \n":
                            log.append(f"Removing structure/plan {structureNames[-1]}/{planNames[-1]} due to empty DVH data in {filename}.")
                            structureNames.pop()
                            planNames.pop()
                            emptyStructure = True
                            dontAppendLength = True

                        if not emptyStructure:
                            if "Min Dose [" in line:
                                structureMinDose.append(float(line.split(": ")[-1]))

                            if "Max Dose [" in line:
                                structureMaxDose.append(float(line.split(": ")[-1]))

                            if "Mean Dose [" in line:
                                structureMeanDose.append(float(line.split(": ")[-1]))

                            if "Volume [" in line and not "Total Structure" in line and not "Structure Volume" in line:
                                structureVolume.append(float(line.split(": ")[-1]))

                        if "Dose [" in line and "Volume [" in line:
                            if not doseHeaderLine:
                                doseHeaderLine = line

                            if not emptyStructure:
                                structureStarts.append(idx + 1)
                                planStarts.append(idx + 1)
                            else:
                                emptyStructure = False

                        idx += 1

                        if line == "\n":
                            lastLineEmpty = True
                        else:
                            lastLineEmpty = False

                # To get to end of file
                structureLength.append(1000000)
                planLength.append(1000000)

                structures = zip(structureNames, structureStarts, structureLength, structureMeanDose, structureMinDose, structureMaxDose, structureVolume)
                plans = zip(planNames, planStarts, planLength)

                for structure, plan in zip(structures, plans):
                    if match(self.options.structureToUse.get(), structure[0]) and match(self.options.planToUse.get(), plan[0]):
                        if self.options.autodetectDVHHeader.get():
                            header_dict = {"Dose [Gy]": "Dose", "Dose [cGy]": "Dose", "Ratio of Total Structure Volume [%]": "Volume", 
                            "Structure Volume [cmÂ³]": "Volume", "Relative dose [%]": "Relative dose"}
                            
                            possible_headers = list(header_dict.keys())

                            header_index = list()
                            for header in possible_headers:
                                try:
                                    header_index.append(doseHeaderLine.index(header))
                                except:
                                    header_index.append(99)
                                    pass

                            header_index = [header_index.index(k) for k in sorted(header_index) if k != 99]
                            identified_headers = np.array(possible_headers)[header_index]
                            headers = [header_dict[k] for k in np.array(possible_headers)[header_index] if k in header_dict]
                            if not printedLogOutput:
                                #printedLogOutput = True
                                log.append(f"Identified the following ECLIPSE DVH header structure for file {filename}: {headers}")

                        else:
                            try:
                                headers = self.options.customDVHHeader.get().split(",")
                            except:
                                print(f"Could not identify structure: {doseHeaderLine}, please input custom DVH header.")
                                headers = ["Dose", "Relative dose", "Volume"]

                        dvh = pd.read_csv(self.getFilePath(filename), header=None, names=headers,
                                          usecols=["Dose", "Volume"], decimal=".", sep="\s+",
                                          skiprows=structure[1], nrows=structure[2], engine="python")

                        if np.sum(dvh.isnull().values):
                            headers = ["Dose", "Volume"]
                            dvh = pd.read_csv(self.getFilePath(filename), header=None, names=headers,
                                              usecols=["Dose", "Volume"], decimal=".", sep="\s+",
                                              skiprows=structure[1], nrows=structure[2], engine="python")

                        # Create Patient object with DVH data
                        if self.options.doseUnit.get() == 'autodetect' and not self.doseUnit:
                            maxDose = dvh["Dose"].max()
                            if maxDose > 1000:
                                doseUnit = self.options.cGy
                            else:
                                doseUnit = self.options.Gy
                            self.doseUnit = doseUnit
                        elif self.doseUnit:
                            doseUnit = self.doseUnit
                        else:
                            doseUnit = self.options.doseUnit.get() == "cGy" and self.options.cGy or self.options.Gy

                        n += 1
                        dvh["Dose"] = dvh["Dose"] * doseUnit
                        maxVolume = dvh["Volume"].at[0]
                        dvh["Volume"] = dvh["Volume"] * 100 / maxVolume

                        dvh = dvh[dvh["Volume"] > 0]

                        patient = Patient(dvh)
                        if self.options.loadToxFromFilename.get():
                            patient.setTox("tox" in filename.lower())
                        patient.setStructure(structure[0])
                        patient.setPlan(plan[0])
                        patientName = filename[:-4]
                        patient.setCohort(self.dataFolder.split("/")[-1])
                        patient.setDataFolder(self.dataFolder)
                        patient.setID(f"{patientName}_{patient.getPlan()}_{patient.getStructure()}")
                        patient.setMeanDoseFromEclipse(structure[3] * doseUnit)
                        patient.setMinDoseFromEclipse(structure[4] * doseUnit)
                        patient.setMaxDoseFromEclipse(structure[5] * doseUnit)
                        patient.setVolumeFromEclipse(structure[6])

                        # Add object to dictionary
                        self.cohort = self.dataFolder.split("/")[-1]
                        self.structure = structure[0]
                        self.patients[patient.getID()] = patient

        progress['value'] = 0
        return log

    def loadPatientsRayStation(self, progress):
        def match(a, b):
            a += "$"  # Don't match end-of-lines
            a = a.replace("*", ".*")  # Use regex type wildcard
            return re.match(a, b)

        log = []

        if not self.options.loadToxFromFilename.get():
            try:
                with open("%s_tox.csv" % (self.dataFolder), "r") as toxFile:
                    toxDict = {}
                    for line in toxFile.readlines():
                        linesplit = line.split(",")
                        toxDict[linesplit[0]] = int(linesplit[1])
            except:
                log.append(f"Cannot open toxicity file {self.dataFolder}_tox.csv, please include in format (filename,tox).")
                return log

        n = 0

        progress['maximum'] = self.findNPatients()

        for root, d, f, in os.walk(self.dataFolder):
            for filename in f:
                if not filename[-4:] in [".txt", ".dvh"]:
                    continue
                if "gEUD" in filename:
                    continue
                idx = 0
                structureNames = []
                structureStarts = []
                structureLength = []
                structureMeanDose = []
                structureMinDose = []
                structureMaxDose = []
                planNames = []
                planStarts = []
                planLength = []
                lastLineEmpty = False
                firstDoseLine = None

                progress.step(1)
                progress.update_idletasks()

                with open(self.getFilePath(filename), "r") as textin:
                    for line in textin:
                        if "#RoiName:" in line:
                            structureNames.append(re.sub("\s+", "", line.split(":")[-1]))
                            if structureStarts:
                                structureLength.append(idx - structureStarts[-1] - 1)

                        if "#Dosename:" in line:
                            planNames.append(re.sub("\s+", "", line.split(":")[-1]))
                            if planStarts:
                                planLength.append(idx - planStarts[-1] - 3)

                        if "#Dose unit: " in line:
                            structureStarts.append(idx + 1)
                            planStarts.append(idx + 1)

                        if not "#" in line and not firstDoseLine:
                            firstDoseLine = [float(k) for k in line.split("\t")]

                        idx += 1

                # To get to end of file
                structureLength.append(1000000)
                planLength.append(1000000)

                structures = zip(structureNames, structureStarts, structureLength)
                plans = zip(planNames, planStarts, planLength)

                for structure in structures:
                    if match(self.options.structureToUse.get(), structure[0]):

                        if self.options.autodetectDVHHeader.get():
                            if len(firstDoseLine) == 2:
                                if firstDoseLine[0] == 0:
                                    headers = ["Dose", "Volume"]
                                else:
                                    headers = ["Volume", "Dose"]
                            else:
                                print(f"Could not identify structure: {firstDoseLine}, please input custom DVH header.")
                                headers = self.options.customDVHHeader.get().split(",")

                        else:
                            headers = self.options.customDVHHeader.get().split(",")

                        dvh = pd.read_csv(self.getFilePath(filename), header=None, names=headers,
                                          usecols=["Dose", "Volume"], decimal=".", sep="\s+",
                                          skiprows=structure[1], nrows=structure[2], engine="python")

                        # Create Patient object with DVH data
                        if self.options.doseUnit.get() == 'autodetect' and not self.doseUnit:
                            maxDose = dvh["Dose"].max()
                            if maxDose > 1000:
                                doseUnit = self.options.cGy
                            else:
                                doseUnit = self.options.Gy
                            self.doseUnit = doseUnit
                        elif self.doseUnit:
                            doseUnit = self.doseUnit
                        else:
                            doseUnit = self.options.doseUnit.get() == "cGy" and self.options.cGy or self.options.Gy

                        n += 1
                        dvh["Dose"] = dvh["Dose"] * doseUnit
                        maxVolume = dvh["Volume"].at[0]
                        dvh["Volume"] = dvh["Volume"] * 100 / maxVolume

                        dvh = dvh[dvh["Volume"] > 0]

                        patient = Patient(dvh)
                        if self.options.loadToxFromFilename.get():
                            patient.setTox("tox" in filename.lower())
                        patient.setStructure(structure[0])
                        # patient.setPlan(plan[0])
                        patient.setPlan(planNames[0])
                        patientName = filename[:-4]
                        patient.setCohort(self.dataFolder.split("/")[-1])
                        patient.setDataFolder(self.dataFolder)
                        patient.setID(f"{patientName}_{patient.getPlan()}_{patient.getStructure()}")

                        # Add object to dictionary
                        self.cohort = self.dataFolder.split("/")[-1]
                        self.structure = structure[0]
                        self.patients[patient.getID()] = patient

        progress['value'] = 0
        return log

    def loadPatientsSIMPLE(self, progress):
        log = []
        if not self.options.loadToxFromFilename.get():
            try:
                with open("%s_tox.csv" % (self.dataFolder), "r") as toxFile:
                    toxDict = {}
                    for line in toxFile.readlines():
                        linesplit = line.split(",")
                        toxDict[linesplit[0]] = int(linesplit[1])
            except:
                log.append("Cannot open toxicity file %s_tox.csv, please include in format (filename,tox)." % (self.dataFolder))
                return log

        nFiles = 0
        for root, d, f, in os.walk(self.dataFolder):
            for filename in f:
                if not "gEUD" in filename:
                    nFiles += 1
        progress['maximum'] = nFiles

        printedLogOutput = False

        for root, d, f, in os.walk(self.dataFolder):  # GENERALIZE THIS IF MORE 'simple' FILES ARE TO BE USED ...
            for filename in f:
                progress.step(1)
                progress.update_idletasks()
                if "gEUD" in filename:
                    continue

                if self.options.CSVStyle.get() == "autodetect":
                    with open(self.getFilePath(filename)) as autodetectFile:
                        for _ in range(self.options.skipRows.get()):
                            autodetectFile.readline()
                        line = autodetectFile.readline()
                        if ";" in line:
                            dec = ","
                            sep = ";"
                        else:
                            dec = "."
                            sep = ","

                elif self.options.CSVStyle.get() == "periodComma":
                    dec = "."
                    sep = ","
                elif self.options.CSVStyle.get() == "commaSemicolon":
                    dec = ","
                    sep = ";"

                if self.options.autodetectDVHHeader.get():
                    with open(self.getFilePath(filename)) as autodetectFile:
                        for _ in range(self.options.skipRows.get()):
                            autodetectFile.readline()
                        line = autodetectFile.readline()

                        firstDoseLine = [float(k) for k in line.split(sep)]
                        if len(firstDoseLine) == 2:
                            if firstDoseLine[0] == 0:
                                headers = ["Dose", "Volume"]
                            else:
                                headers = ["Volume", "Dose"]

                            if not printedLogOutput:
                                printedLogOutput = True
                                log.append(f"Identified the following DVH header structure: {headers}")
                        else:
                            print(f"Could not identify structure: {firstDoseLine}, please input custom DVH header.")
                            headers = self.options.customDVHHeader.get().split(",")

                else:
                    headers = self.options.customDVHHeader.get().split(",")

                headers = self.options.customDVHHeader.get().split(",")
                dvh = pd.read_csv(self.getFilePath(filename), decimal=dec, sep=sep, names=headers,
                                  skiprows=self.options.skipRows.get(), usecols=["Dose", "Volume"], dtype=np.float64)  # skiprows = 2 if name

                if self.options.doseUnit.get() == 'autodetect' and not self.doseUnit:
                    maxDose = dvh["Dose"].max()
                    if maxDose > 1000:
                        doseUnit = self.options.cGy
                    else:
                        doseUnit = self.options.Gy
                    self.doseUnit = doseUnit
                elif self.doseUnit:
                    doseUnit = self.doseUnit
                else:
                    doseUnit = self.options.doseUnit.get() == "cGy" and self.options.cGy or self.options.Gy

                dvh["Dose"] = dvh["Dose"] * doseUnit
                dvh = dvh.dropna(axis=0, how='any')

                if dvh["Volume"].at[0] == 0:
                    dvh = dvh[::-1].reset_index(drop=True)

                maxVolume = dvh["Volume"].at[0]
                dvh["Volume"] = dvh["Volume"] * 100 / maxVolume

                patient = Patient(dvh)
                patientName = filename[:-4]

                patient.setCohort(self.dataFolder.split("/")[-1])
                patient.setDataFolder(self.dataFolder)
                patient.setStructure("CSV")
                patient.setPlan("TPS")
                patient.setID("%s" % (patientName))

                self.cohort = self.dataFolder.split("/")[-1]
                self.structure = "CSV"
                self.patients[patient.getID()] = patient

                if self.options.loadToxFromFilename.get():
                    patient.setTox("tox" in filename.lower())
                else:
                    try:
                        patient.setTox(toxDict[patient.ID])
                    except:
                        log.append("Could not find patient %s in %s_tox.csv, skipping patient." % (patient.ID, self.dataFolder))
                        continue

        progress['value'] = 0
        return log

    def getPatient(self, name):
        return self.patients[name]

    def getNPatients(self):
        return len(self.patients)

    def getNPatientsWithTox(self):
        n = 0
        for name, patient in list(self.patients.items()):
            n += patient.getTox() >= self.options.toxLimit.get()
        return n

    def calculateNTCP(self):
        if not len(self.bestParameters):
            print("Could not find best parameters...")
            return

        if self.options.NTCPcalculation.get() == "LKB":
            n = self.bestParameters[0]
            m = self.bestParameters[1]
            TD50 = self.bestParameters[2]

            for patient in list(self.patients.values()):
                patient.NTCP = HPM((patient.getGEUD(n) - TD50) / (m * TD50))
        else:
            a = self.bestParameters[0]
            b = self.bestParameters[1]

            for patient in list(self.patients.values()):
                patient.NTCP = 1 - 1 / (1 + exp(a + b * patient.getDpercent()))

    def createDifferentialDVHs(self):
        for name, patient in list(self.patients.items()):
            patient.createDifferentialDVH()

    def createGEUDs(self, progress):
        self.createDifferentialDVHs()
        progress['maximum'] = len(self.patients)
        for patient in list(self.patients.values()):
            progress.step(1)
            progress.update_idletasks()
            patient.createGEUDspline(self.options)
        progress['value'] = 0

    def getTox(self, name):
        return self.getPatient(name).getTox()

    def saveTox(self):
        for patient in list(self.patients.values()):
            patient.storedTox = patient.tox

    def restoreTox(self):
        for patient in list(self.patients.values()):
            patient.tox = patient.storedTox

    def getGEUD(self, name, n):
        return self.getPatient(name).getGEUD(n)

    def checkGEUDsplines(self):
        hasGEUDs = True
        for patient in list(self.patients.values()):
            hasGEUDs *= patient.checkGEUDsplines(self.options)

        return hasGEUDs

    def getDoseAtVolume(self, volume):
        for patient in list(self.patients.values()):
            try:
                dose = patient.getDoseAtVolume(volume)
                print("For patient %s, at volume %.2f cc, the dose is %.2f cGy." % (patient.getID(), volume, dose))
            except:
                print("Could not find volume %.2f cm3 in DVH interval." % (volume))

    def getVolumeAtDose(self, dose):
        for patient in list(self.patients.values()):
            try:
                volume = patient.getVolumeAtDose(dose)
                print("For patient %s, at dose %.2f cGy, the volume is %.2f cc." % (patient.getID(), dose, volume))
            except:
                print("Could not find dose %.2f cGy in DVH interval." % (dose))

    def getVolumeAtRelativeDose(self, relativeDose):
        """ Requires ECLIPSE metadata """

        for patient in list(self.patients.values()):
            try:
                maxdose = patient.getMaxDoseFromEclipse()
                volume = patient.getVolumeAtDose(relativeDose / 100 * maxdose)
                print("For patient %s, at dose %.2f cGy, the volume is %.2f cc." % (patient.getID(), dose, volume))
            except:
                print("Could not find dose %.2f cGy in DVH interval." % (dose))

    def calculateDpercent(self, percent):
        for patient in list(self.patients.values()):
            if self.options.useNTCPcc.get() and isinstance(self.volumeFromEclipse, float):
                patient.calculateDcc(percent)  # percent = cc now
            else:
                patient.calculateDpercent(percent)

    def calculateTDxFromLogit(self, percent, a=None, b=None):
        X_test = np.linspace(0, 100, 300)
        if not a:
            a = self.bestParameters[0]
        if not b:
            b = self.bestParameters[1]
        NTCP = np.array([1 - 1 / (1 + exp(a + b * k)) for k in X_test])
        TDx = X_test[np.argmin(abs(NTCP - percent / 100))]
        return TDx
