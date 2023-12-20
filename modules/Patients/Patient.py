import os
from math import *
import pandas as pd
import numpy as np
import numba as nb
from bisect import bisect_left
from tkinter import *
from matplotlib import pyplot as plt

@nb.njit
def njitGetGEUD(n, nList, GEUDlist):
	# Fast interpolation (much quicker than np.interp1d)
	# stackoverflow.com/questions/35508571/multiple-1d-interpolations-in-python
	index = np.searchsorted(nList, n)
	_xrange = nList[index] - nList[index - 1]
	xdiff = n - nList[index - 1]
	modolo = xdiff / _xrange
	ydiff = GEUDlist[index] - GEUDlist[index - 1]
	lastGEUD = GEUDlist[index - 1] + modolo * ydiff
	return n, lastGEUD

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
		self.nListFirst = 0
		self.nListLast = 2

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
		self.dvh.columns = ["Dose", "Volume", "Avg. Dose", "Diff. Volume"]

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
					self.nListFirst = self.nList[0]
					self.nListLast = self.nList[-1]
					options.nFrom.set(nList[0])
					options.nTo.set(nList[-1] - 0.02)

				return True
		except:
				return False

	def createGEUDspline(self, options):
		def integrate_eud(row):
				return row["Avg. Dose"] ** ninv * row["Diff. Volume"]

		# if not self.checkGEUDsplines(options):

		if options.nIsLinear.get():
				nList = np.linspace(options.nFrom.get() * 0.8, options.nTo.get() * 1.2, options.nGrid.get())
		else:
				nList = np.logspace(np.log10(options.nFrom.get() * 0.8), np.log10(options.nTo.get() * 1.2), options.nGrid.get())

		GEUDlist = []
		nListGood = []
		for n in nList:
				ninv = 1 / n
				EUD = self.dvh.apply(integrate_eud, axis=1)
				GEUD = EUD.sum() ** n
				if not isinf(GEUD):
					GEUDlist.append(GEUD)
					nListGood.append(n)
		try:
				self.nList = np.array(nListGood)
				self.GEUDlist = np.array(GEUDlist)
				self.nListFirst = self.nList[0]
				self.nListLast = self.nList[-1]
		except:
				print(f"Could not create gEUD list for patient {self.ID}.")
				self.nList = None
				self.GEUDlist = None

		# Save result to CSV files
		if not os.path.exists(f"{self.dataFolder}/gEUD/"):
				os.makedirs(f"{self.dataFolder}/gEUD/")
		with open(f"{self.dataFolder}/gEUD/gEUD_{self.cohort}_{self.structure}_{self.ID}.csv", "wb") as fh:
				for n, GEUD in zip(self.nList, self.GEUDlist):
					fh.write(b"%.3f,%.3f\n" % (n, GEUD))

	def fastGetGEUD(self, n):
		if self.lastN == n:
				return self.lastGEUD

		elif not self.nListFirst <= n <= self.nListLast:
				return 0
						
		self.lastN, self.lastGEUD = njitGetGEUD(n, self.nList, self.GEUDlist)
		return self.lastGEUD

	def getGEUD(self, n):
		if self.lastN == n:
				return self.lastGEUD

		elif not self.nListFirst <= n <= self.nListLast:
				return 0

		else:
				# Fast interpolation (much quicker than np.interp1d)
				# stackoverflow.com/questions/35508571/multiple-1d-interpolations-in-python
				index = bisect_left(self.nList, n)
				_xrange = self.nList[index] - self.nList[index - 1]
				xdiff = n - self.nList[index - 1]
				modolo = xdiff / _xrange
				ydiff = self.GEUDlist[index] - self.GEUDlist[index - 1]
				self.lastN = n
				self.lastGEUD = self.GEUDlist[index - 1] + modolo * ydiff

		return self.lastGEUD

	def getDoseAtVolume(self, volume):
		dvh_sort = self.dvh.sort_values(by=["Volume"])
		return np.interp(volume, dvh_sort["Volume"], dvh_sort["Dose"])

	def getVolumeAtDose(self, dose):
		if self.dvh.loc[self.dvh["Dose"] > dose]["Volume"].sum() == 0:
				return 0
		return np.interp(dose, self.dvh["Dose"], self.dvh["Volume"])

	def getVolumeAtRelativeDose(self, relativeDose):
		""" Requires ECLIPSE metadata """

		maxdose = self.getMaxDoseFromEclipse()
		volume = self.getVolumeAtDose(relativeDose / 100 * maxdose)
		return volume

	def calculateDpercent(self, percent):
		Dpercent = self.getDoseAtVolume(percent)
		self.Dpercent = Dpercent

	def calculateDcc(self, cc):
		percent = cc / self.volumeFromEclipse * 100
		Dpercent = self.getDoseAtVolume(percent)
		self.Dpercent = Dpercent

	def getDpercent(self):
		return self.Dpercent

	def getSIndex(self):
		self.createDifferentialDVH()

		D = self.dvh["Dose"].values
		V = self.dvh["Diff. Volume"].values

		D_linspace = np.linspace(np.min(D), np.max(D), len(D))

		Dsd = 0
		Dmean = 0
		V_interp_list = np.zeros(len(D))

		D_width = D_linspace[1] - D_linspace[0]

		for i in range(len(D_linspace)):
			Di = D_linspace[i]
			Dlo = Di - D_width/2
			Dhi = Di + D_width/2
			mask = (D > Dlo) & (D < Dhi)
			V_interp_list[i] = np.sum(V[mask])
			Dmean += V_interp_list[i] * D_linspace[i]
			
		Dmean /= np.sum(V_interp_list)
	
		for i in range(len(D_linspace)):
			Dsd += np.sqrt((D_linspace[i] - Dmean)**2) * V_interp_list[i]

		return Dsd
