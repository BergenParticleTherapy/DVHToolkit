import os
from math import *
import pandas as pd
import numpy as np
from bisect import bisect_left
from tkinter import *


class Patient:
    def __init__(self, dvh):
        self.dvh = dvh
        self.gEUD = {}
        self.bestGEUD = None
        self.structure = None
        self.plan = None
        self.cohort = None
        self.NTCP = None
        self.dataFolder = None
        self.ID = None
        self.tox = None
        self.storedTox = None
        self.meanDoseFromEclipse = None
        self.minDoseFromEclipse = None
        self.maxDoseFromEclipse = None
        self.volumeFromEclipse = None
        self.lastN = None
        self.lastGEUD = None
        self.GEUDlist = None
        self.nList = None

    def __repr__(self):
        return {'cohort': self.cohort, 'plan': self.plan, 'structure': self.structure}

    def setTox(self, tox):
        self.tox = tox

    def getTox(self):
        return self.tox

    def getMeanDoseFromEclipse(self):
        return self.meanDoseFromEclipse

    def getMinDoseFromEclipse(self):
        return self.minDoseFromEclipse

    def getMaxDoseFromEclipse(self):
        return self.maxDoseFromEclipse

    def getVolumeFromEclipse(self):
        return self.volumeFromEclipse

    def setMeanDoseFromEclipse(self, d):
        self.meanDoseFromEclipse = d

    def setMinDoseFromEclipse(self, d):
        self.minDoseFromEclipse = d

    def setMaxDoseFromEclipse(self, d):
        self.maxDoseFromEclipse = d

    def setVolumeFromEclipse(self, d):
        self.volumeFromEclipse = d

    def getNTCP(self):
        return self.NTCP

    def setStructure(self, structure):
        self.structure = structure

    def setPlan(self, plan):
        self.plan = plan

    def getStructure(self):
        return self.structure

    def getPlan(self):
        return self.plan

    def setCohort(self, cohort):
        self.cohort = cohort

    def setID(self, ID):
        self.ID = ID

    def setDataFolder(self, folder):
        self.dataFolder = folder

    def getDataFolder(self):
        return self.dataFolder

    def getID(self):
        return self.ID

    def createDifferentialDVH(self):
        if "Avg. Dose" in self.dvh.columns and "Diff. Volume" in self.dvh.columns:
            return

        nrows = self.dvh.shape[0]
        avgDoseList = np.zeros(nrows)
        diffVolumeList = np.zeros(nrows)
        doses = self.dvh["Dose"]
        volumes = self.dvh["Volume"]

        maxVolume = volumes.at[0]
        volumes = volumes / maxVolume

        for idx in range(nrows):
            try:
                avgDose = np.mean([doses[idx], doses[idx + 1]])
            except:
                avgDose = doses[idx]
            try:
                diffVolume = volumes[idx] - volumes[idx + 1]
            except:
                diffVolume = 0

            avgDoseList[idx] = avgDose
            diffVolumeList[idx] = diffVolume

        # Add the two new columns (avg dose, differential volume) to the dataframe
        self.dvh = self.dvh.assign(avgDose=avgDoseList, diffVolume=diffVolumeList)
        print(self.dvh.columns)
        self.dvh.columns = ["Dose", "Volume", "Avg. Dose", "Diff. Volume"]
        #self.dvh.columns.append("Avg. Dose")
        #self.dvh.columns.append("Diff. Volume")

        # Remove rows where the differential volume is zero, not needed for calculation, reduce data set by 10%-20%
        self.dvh = self.dvh[self.dvh["Diff. Volume"] > 0]
        self.dvh = self.dvh.reset_index()

    def checkGEUDsplines(self, options):
        try:
            with open(f"{self.dataFolder}/gEUD/gEUD_{self.cohort}_{self.structure}_{self.ID}.csv", "r") as fh:
                GEUDlist = []
                nList = []
                for line in fh.readlines():
                    linesplit = line.split(",")
                    nList.append(float(linesplit[0]))
                    GEUDlist.append(float(linesplit[1]))
                self.nList = np.array(nList)
                self.GEUDlist = np.array(GEUDlist)
                options.nFrom.set(nList[0])
                options.nTo.set(nList[-1] - 0.02)

            return True
        except:
            return False

    def createGEUDspline(self, options):
        def integrate_eud(row):
            return row["Avg. Dose"] ** ninv * row["Diff. Volume"]

        if not self.checkGEUDsplines(options):
            nList = np.arange(min(0.05, options.nFrom.get() * 0.8), max(options.nTo.get() * 1.2, 2), 0.02)
            # nList = np.arange(0.05, 2.0, 0.05)
            GEUDlist = []
            for n in nList:
                ninv = 1 / n
                EUD = self.dvh.apply(integrate_eud, axis=1)
                GEUD = EUD.sum() ** n
                GEUDlist.append(GEUD)
            try:
                self.nList = np.array(nList)
                self.GEUDlist = np.array(GEUDlist)
            except:
                print(f"Could not create gEUD list for patient {self.ID}.")
                self.nList = None
                self.GEUDlist = None

            # Save result to CSV files
            if not os.path.exists(f"{self.dataFolder}/gEUD/"):
                os.makedirs(f"{self.dataFolder}/gEUD/")
            with open(f"{self.dataFolder}/gEUD/gEUD_{self.cohort}_{self.structure}_{self.ID}.csv", "wb") as fh:
                for n, GEUD in zip(nList, GEUDlist):
                    fh.write(b"%.3f,%.3f\n" % (n, GEUD))

    def getGEUD(self, n):
        if self.lastN == n:
            return self.lastGEUD

        else:
            # Fast interpolation (much quicker than np.interp1d)
            # stackoverflow.com/questions/35508571/multiple-1d-interpolations-in-python
            try:
                index = bisect_left(self.nList, n)
                _xrange = self.nList[index] - self.nList[index - 1]
                xdiff = n - self.nList[index - 1]
                modolo = xdiff / _xrange
                ydiff = self.GEUDlist[index] - self.GEUDlist[index - 1]
                self.lastN = n
                self.lastGEUD = self.GEUDlist[index - 1] + modolo * ydiff
            except IndexError as e:
                print(f"Error: The wanted value n={n} is not pre-calculated: {e}")
                print(f"The pre-calculated range is [{self.nList[0]} - {self.nList[-1]}].")
                return 0

        return self.lastGEUD

    def getDoseAtVolume(self, volume):
        dvh_sort = self.dvh.sort_values(by=["Volume"])
        return np.interp(volume, dvh_sort["Volume"], dvh_sort["Dose"])

    def getVolumeAtDose(self, dose):
        if self.dvh.loc[self.dvh["Dose"] > dose]["Volume"].sum() == 0:
            return 0
        return np.interp(dose, self.dvh["Dose"], self.dvh["Volume"])

    def calculateDpercent(self, percent):
        Dpercent = self.getDoseAtVolume(percent)
        self.Dpercent = Dpercent

    def calculateDcc(self, cc):
        percent = cc / self.volumeFromEclipse * 100
        Dpercent = self.getDoseAtVolume(percent)
        self.Dpercent = Dpercent

    def getDpercent(self):
        return self.Dpercent
