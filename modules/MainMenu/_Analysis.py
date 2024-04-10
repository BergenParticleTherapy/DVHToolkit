import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from operator import itemgetter
from math import *
from tkinter import *
from ..Patients import *

def perc(num):
	if num % 10 == 1:
		return f"{num:.0f}st"
	elif num % 10 == 2:
		return f"{num:.0f}nd"
	elif num % 10 == 3:
		return f"{num:.0f}rd"
	else:
		return f"{num:.0f}th"

def calculateDVHvalues(self):
	# Add here for specific organ-specific gEUD calculations to DVH output file
	nValues = {'Bladder': 1 / 8, 'Rectum': 1 / 12, 'Intestine': 1 / 4}

	patients = set()
	for patientsInCohort in self.patients.values():
		for name in patientsInCohort.patients.keys():
				patients.add(name)

	columns = ["Cohort", "Name", "Plan", "Structure"]
	csv = pd.DataFrame(index=sorted(list(patients)), columns=columns)
	for cohortName, patientsInCohort in self.patients.items():
		for name, patient in list(patientsInCohort.patients.items()):
				csv.loc[name, "Plan"] = patient.getPlan()
				csv.loc[name, "Structure"] = patient.getStructure()
				csv.loc[name, "Cohort"] = cohortName
				csv.loc[name, "Name"] = name.split("_")[0]

				# ECLIPSE Dose Metrics
				if self.calculateMeanDose.get():
					csv.loc[name, "ECLIPSEMeanDose [Gy]"] = patient.getMeanDoseFromEclipse()
					csv.loc[name, "ECLIPSEMinDose [Gy]"] = patient.getMinDoseFromEclipse()
					csv.loc[name, "ECLIPSEMaxDose [Gy]"] = patient.getMaxDoseFromEclipse()
					csv.loc[name, "ECLIPSEVolume [cc]"] = patient.getVolumeFromEclipse()

				# gEUD
				if self.dvhCheckVarIncludeGEUD.get():
					if not isinstance(patient.GEUDlist, list):
						patient.createDifferentialDVH()
						patient.createGEUDspline(self.options)

					csv.loc[name, "Fixed n-value"] = self.options.nSet.get()
					csv.loc[name, "gEUD [Gy]"] = patient.getGEUD(self.options.nSet.get())

					"""
					elif patient.getStructure().capitalize() in nValues:
						csv.loc[name, "Structure n-value"] = nValues[patient.getStructure()]
						csv.loc[name, "gEUD [Gy]"] = patient.getGEUD(nValues[patient.getStructure()])
					else:
						self.log(f"No seriality parameter found or set for patient/structure {name}/{patient.getStructure()}.")
						self.log("(Other structures are still saved).")
					"""

				# DVH Dose Metrics
				if self.dvhCheckVarDoseAtVolume.get():
					for volume in self.dvhEntryVar1.get().split(","):
						csv.loc[name, f"D{float(volume):g}%"] = patient.getDoseAtVolume(float(volume))
				if self.dvhCheckVarVolumeAtDose.get():
					for dose in self.dvhEntryVar2.get().split(","):
						csv.loc[name, f"V{float(dose):g}Gy"] = patient.getVolumeAtDose(float(dose))
				if self.dvhCheckVarVolumeAtRelDose.get():
					for relativeDose in self.dvhEntryVar3.get().split(","):
						csv.loc[name, f"V{float(dose):g}%"] = patient.getVolumeAtRelativeDose(float(relativeDose))
				
				# NTCP
				if self.dvhCheckVarIncludeNTCP.get():
					if len(self.bestParameters) == 0:
						if self.options.NTCPcalculation.get() == "LKB":
								if self.options.fixN.get() + self.options.fixM.get() + self.options.fixTD50.get() < 3:
									self.log("Must choose fixed N, M, TD50 for LKB/NTCP calculation when no fit is performed")
									self.dvhCheckVarIncludeNTCP.set(0)
								else:
									n = self.options.nFrom.get()
									m = self.options.mFrom.get()
									TD50 = self.options.TD50From.get()
									patient.NTCP = HPM((patient.getGEUD(n) - TD50) / (m * TD50))

						elif self.options.NTCPcalculation.get() == "Logit":
								if self.options.fixA.get() + self.options.fixB.get() < 2:
									self.log("Must choose fixed a, b for logit NTCP calculation when no fit is performed")
									self.dvhCheckVarIncludeNTCP.set(0)

								else:
									a = self.options.aFrom.get()
									b = self.options.bFrom.get()
									patient.NTCP = 1 - 1 / (1 + exp(a + b * patient.getDpercent()))

					csv.loc[name, f"NTCP ({self.options.NTCPcalculation.get()})"] = patient.NTCP

				# HI / S-index
				if self.dvhCheckVarCalculateSIndex.get():
					csv.loc[name, "S-index"] = patient.getSIndex()

	directory = "/".join(self.outputFileNameVar.get().split("/")[:-1])
	if not os.path.exists(directory):
		self.log(f"Making directory {directory}.")
		os.makedirs(directory)

	if "csv" in self.outputFileNameVar.get():
		csv.to_csv(self.outputFileNameVar.get(), index=False)
	elif "xlsx" in self.outputFileNameVar.get():
		csv.to_excel(self.outputFileNameVar.get())
	else:
		self.log("Please specify valid output file type (.csv or .xlsx).")

	self.window.destroy()


def calculateAggregatedDVH(self):
	self.window.destroy()
	cohortDVH = {}
	tox = {}
	notox = {}

	# CALCULATION
	#############

	doses = {}
	for cohort_, patientsInCohort in self.patients.items():
		for name, patient in patientsInCohort.patients.items():
				plan = patient.getPlan()
				structure = patient.getStructure()
				cohort = f"{structure}/{plan}"
				if cohort not in doses:
					doses[cohort] = pd.DataFrame({"Dose": [0.0]})
					doses[cohort].set_index("Dose", inplace=True)

				doses[cohort] = doses[cohort].merge(patient.dvh["Dose"], how="outer", on="Dose")
				doses[cohort] = doses[cohort].sort_values("Dose").drop_duplicates().reset_index(drop=True)
				
				"""
				structure = cohort.split("/")[0]
				doses[cohort].to_excel(f"Output/doses_{structure}.xlsx")
				"""

	if self.dvhStyleVar2.get() == "showAll": # compare tox vs no tox
		first = True
		for cohort_, patientsInCohort in list(self.patients.items()):
				for name, patient in list(patientsInCohort.patients.items()):
					plan = list(patientsInCohort.patients.values())[0].getPlan()
					structure = list(patientsInCohort.patients.values())[0].getStructure()
					cohort = f"{structure}/{plan}"
					if cohort not in tox:
						tox[cohort] = []
						notox[cohort] = []
					if patient.getTox() >= self.options.toxLimit.get():
						tox[cohort].append(f"Volume_{name}")
					else:
						notox[cohort].append(f"Volume_{name}")

					namePx = name.split("_")[0]
					namePx = namePx[:-1]

					if first:
						cohortDVH[cohort] = pd.DataFrame({"Dose": patient.dvh["Dose"], f"Volume_{name}": patient.dvh["Volume"]})
						cohortDVH[cohort].set_index("Dose", inplace=True)
						first = False
					else:
						newDVH = pd.DataFrame({"Dose": patient.dvh["Dose"], f"Volume_{name}": patient.dvh["Volume"]})
						newDVH.set_index("Dose", inplace=True)
						cohortDVH[cohort] = cohortDVH[cohort].merge(newDVH, how="outer", right_index=True, left_index=True)

				cohortDVH[cohort] = cohortDVH[cohort].interpolate(method='index', limit_direction='backward', limit=100).fillna(0)

				cohortTox = cohortDVH[cohort][tox[cohort]]
				cohortNoTox = cohortDVH[cohort][notox[cohort]]

				if self.dvhStyleVar1.get() == "mean":
					cohortDVH[cohort]["Volume agg"] = cohortDVH[cohort].mean(axis=1)
					cohortDVH[cohort]["Volume agg tox"] = cohortTox.mean(axis=1)
					cohortDVH[cohort]["Volume agg notox"] = cohortNoTox.mean(axis=1)
				else:
					cohortDVH[cohort]["Volume agg"] = cohortDVH[cohort].median(axis=1)
					cohortDVH[cohort]["Volume agg tox"] = cohortTox.median(axis=1)
					cohortDVH[cohort]["Volume agg notox"] = cohortNoTox.median(axis=1)

				ql = (1 - self.dvhConfidenceInterval.get() / 100) / 2
				qu = 1 - ql
				cohortDVH[cohort]["Volume agg tox CI lower"] = cohortTox.quantile(ql, axis=1)
				cohortDVH[cohort]["Volume agg tox CI upper"] = cohortTox.quantile(qu, axis=1)
				cohortDVH[cohort]["Volume agg notox CI lower"] = cohortNoTox.quantile(ql, axis=1)
				cohortDVH[cohort]["Volume agg notox CI upper"] = cohortNoTox.quantile(qu, axis=1)

				qps = self.dvhPercentileInterval.get().split(",")
				for qp_raw in qps:
					qp = int(qp_raw) / 100
					qp_str = perc(qp*100)
					cohortDVH[cohort][f"Volume agg tox {qp_str} percentile"] = cohortTox.quantile(qp, axis=1)
					cohortDVH[cohort][f"Volume agg notox {qp_str} percentile"] = cohortNoTox.quantile(qp, axis=1)

		if self.dvhSaveAggDVH.get() == 1:
				if not os.path.exists("Output/aggDVH"):
					os.makedirs("Output/aggDVH")
				
				for cohort in cohortDVH.keys():
					st = self.dvhStyleVar1.get()
					cohortStructure, cohortPlan = cohort.split("/")
					q = self.dvhConfidenceInterval.get()
					filename = f"Output/aggDVH/{cohortPlan}_{cohortStructure}_{st}.csv"
					splice = cohortDVH[cohort][["Volume agg", "Volume agg tox", "Volume agg notox", 
														 "Volume agg tox CI lower", "Volume agg tox CI upper",
														 "Volume agg tox percentile", "VOlume agg notox percentile"]]

					splice.to_csv(filename, sep=";", decimal=".", index=True, header=True, na_rep="")

	elif self.dvhStyleVar2.get() == "comparePlans":  # Compare within Patient
		for cohort_, patientsInCohort in list(self.patients.items()):
				cohorts = set()
				for name, patient in list(patientsInCohort.patients.items()):
					structure = patient.getStructure()
					namePx = patient.getID().split("_")[0]
					cohort = f"{structure}/{namePx}"
					cohorts.add(cohort)

					if cohort not in cohortDVH:
						cohortDVH[cohort] = pd.DataFrame({"Dose": patient.dvh["Dose"], f"Volume_{name}": patient.dvh["Volume"]})
						cohortDVH[cohort].set_index("Dose", inplace=True)
					else:
						# newDVH gets very very large
						newDVH = pd.DataFrame({"Dose": patient.dvh["Dose"], f"Volume_{name}": patient.dvh["Volume"]})
						newDVH.set_index("Dose", inplace=True)
						cohortDVH[cohort] = cohortDVH[cohort].merge(newDVH, how="outer", right_index=True, left_index=True)

				for cohort in cohorts:
					cohortDVH[cohort] = cohortDVH[cohort].interpolate(method='index', limit_direction='backward', limit=100).fillna(0)

					if self.dvhStyleVar1.get() == "mean":
						cohortDVH[cohort]["Volume agg"] = cohortDVH[cohort].mean(axis=1)
					else:
						cohortDVH[cohort]["Volume agg"] = cohortDVH[cohort].median(axis=1)

					ql = (1 - self.dvhConfidenceInterval.get() / 100) / 2
					qu = 1 - ql
					cohortDVH[cohort]["Volume agg CI lower"] = cohortDVH[cohort].quantile(ql, axis=1)
					cohortDVH[cohort]["Volume agg CI upper"] = cohortDVH[cohort].quantile(qu, axis=1)

					qps = self.dvhPercentileInterval.get().split(",")
					for qp_raw in qps:
						qp = int(qp_raw) / 100
						qp_str = perc(qp*100)
						cohortDVH[cohort][f"Volume agg {qp_str} percentile"] = cohortDVH[cohort].quantile(qp, axis=1)

		if self.dvhSaveAggDVH.get() == 1:
				if not os.path.exists("Output/aggDVH"):
					os.makedirs("Output/aggDVH")
				
				for cohort in cohortDVH.keys():
					st = self.dvhStyleVar1.get()
					q = self.dvhConfidenceInterval.get()
					cohortStructure, cohortPlan = cohort.split("/")
					q = self.dvhConfidenceInterval.get()
					filename = f"Output/aggDVH/{cohortPlan}_{cohortStructure}_{st}.csv"
					splice = cohortDVH[cohort][["Volume agg", "Volume agg CI lower", 
														 "Volume agg CI upper", "Volume agg percentile"]]
					splice.to_csv(filename, sep=";", decimal=".", 
										index=True, header=True, na_rep="")

	else:  # Compare across patients
		for cohort_, patientsInCohort in self.patients.items():
				for name, patient in patientsInCohort.patients.items():
					plan = patient.getPlan()
					structure = patient.getStructure()
					cohort = f"{structure}/{plan}"

					if cohort not in cohortDVH:
						cohortDVH[cohort] = pd.DataFrame({"Dose": doses[cohort]["Dose"]})
						cohortDVH[cohort].set_index("Dose", inplace=True)

					newDVH = pd.DataFrame({"Dose": patient.dvh["Dose"], f"Volume_{name}": patient.dvh["Volume"]})
					newDVH.set_index("Dose", inplace=True)
					cohortDVH[cohort] = cohortDVH[cohort].merge(newDVH, how="left", on="Dose")

		for cohort in cohortDVH.keys():
				cohortDVH[cohort] = cohortDVH[cohort].interpolate(method='index', limit_direction='backward', limit = 3000).fillna(0)
				
				structure = cohort.split("/")[0]
				cohortDVH[cohort].to_excel(f"Output/{structure}.xlsx")

				if self.dvhStyleVar1.get() == "mean":
					cohortDVH[cohort]["Volume agg"] = cohortDVH[cohort].mean(axis=1)  # tox and notox
				else:
					cohortDVH[cohort]["Volume agg"] = cohortDVH[cohort].median(axis=1)

				ql = (1 - self.dvhConfidenceInterval.get() / 100) / 2
				qu = 1 - ql
				cohortDVH[cohort]["Volume agg CI lower"] = cohortDVH[cohort].quantile(ql, axis=1)
				cohortDVH[cohort]["Volume agg CI upper"] = cohortDVH[cohort].quantile(qu, axis=1)

				qps = self.dvhPercentileInterval.get().split(",")
				for qp_raw in qps:
					qp = int(qp_raw) / 100
					qp_str = perc(qp*100)
					cohortDVH[cohort][f"Volume agg {qp_str} percentile"] = cohortDVH[cohort].quantile(qp, axis=1)

		if self.dvhSaveAggDVH.get() == 1:
				if not os.path.exists("Output/aggDVH"):
					os.makedirs("Output/aggDVH")
				
				for cohort in cohortDVH.keys():
					st = self.dvhStyleVar1.get()
					q = self.dvhConfidenceInterval.get()
					cohortStructure, cohortPlan = cohort.split("/")

					filename = f"Output/aggDVH/{cohortPlan}_{cohortStructure}_{st}.csv"
					splice = cohortDVH[cohort][["Volume agg", "Volume agg CI lower", 
														 "Volume agg CI upper", "Volume agg percentile"]]
					splice.to_csv(filename, sep=";", decimal=".", index=True, header=True, na_rep="")

	###############
	# PLOTTING    #
	###############

	if self.dvhStyleVar2.get() == "showAll":  # Show all aggregated cohorts
		for k, v in cohortDVH.items():
				v["Volume agg tox"].plot(use_index=True, linestyle="-", color="r", label=f"{k} with toxicity")
				v["Volume agg notox"].plot(use_index=True, linestyle="-", color="k", label=f"{k} no toxicity")
		plt.xlabel("Dose [Gy]")
		plt.ylabel("Volume [%]")

		if self.dvhStyleVar3.get():
				plt.fill_between(v.index, v["Volume agg notox CI lower"], v["Volume agg notox CI upper"], color="k", alpha=0.2)
				plt.plot(v["Volume agg notox CI lower"].index, v["Volume agg notox CI lower"], linestyle="-", color="k", linewidth=1, alpha=0.7)
				plt.plot(v["Volume agg notox CI upper"].index, v["Volume agg notox CI upper"], linestyle="-", color="k", linewidth=1, alpha=0.7)

				plt.fill_between(v.index, v["Volume agg tox CI lower"], v["Volume agg tox CI upper"], color="r", alpha=0.2)
				plt.plot(v["Volume agg tox CI lower"].index, v["Volume agg tox CI lower"], linestyle="-", color="r", linewidth=1, alpha=0.7)
				plt.plot(v["Volume agg tox CI upper"].index, v["Volume agg tox CI upper"], linestyle="-", color="r", linewidth=1, alpha=0.7)

		if self.dvhStyleVarPercentile.get():
			qps = self.dvhPercentileInterval.get().split(",")
			for qp_raw in qps:
				qp = int(qp_raw) / 100
				qp_str = perc(qp*100)
				plt.plot(v[f"Volume agg tox {qp_str} percentile"].index, 
							v[f"Volume agg tox {qp_str} percentile"], 
							linestyle="--", 
							color="r", 
							linewidth=2, 
							label=f"{qp_str} percentile")

				plt.plot(v[f"Volume agg notox {qp_str} percentile"].index, 
							v[f"Volume agg notox {qp_str} percentile"], 
							linestyle="--", 
							color="k", 
							linewidth=2, 
							label="{qp_str} percentile")

			if len(qps) > 1:
				qp_lower = int(qps[0]) / 100
				qp_lower_str = perc(qp_lower*100)
				qp_upper = int(qps[-1]) / 100
				qp_upper_str = perc(qp_upper*100)
				plt.fill_between(v.index, 
									  v["Volume agg notox {qp_lower_str} percentile"], 
									  v["Volume agg notox {qp_upper_str} percentile"], 
									  color="k", 
									  alpha=0.2)

				plt.fill_between(v.index, 
									  v["Volume agg tox {qp_lower_str} percentile"], 
									  v["Volume agg tox {qp_upper_str} percentile"], 
									  color="r", 
									  alpha=0.2)

		# plt.xlim([0, 75])
		plt.legend()
		plt.show()

	elif self.dvhStyleVar2.get() == "comparePlans":  # Within patient
		for k, v in cohortDVH.items():
				c = self.colorVarList[k.split("/")[0]].get()
				if ":" in c:
					alpha1 = float(c.split(":")[-1])
					alpha2 = min(1, alpha1 + 0.2)
					c = c.split(":")[0]
				else:
					alpha1 = 0.3
					alpha2 = 0.7

				if self.dvhStyleBlackPlot.get():
					c = "k"

				structure = k.split("/")[0]
				if self.dvhStyleSinglePlot.get():
					fignum = f"{self.dvhStyleVar1.get().capitalize()} DVH for all structures"
				else:
					fignum = f"{self.dvhStyleVar1.get().capitalize()} DVH for {structure}"
				plt.figure(figsize=(10, 7.5), num=fignum)

				v["Volume agg"].plot(use_index=True, linestyle="-", color=c, label=k.split("/")[1])

				if self.dvhStyleVar3.get():
					plt.fill_between(v.index, v["Volume agg CI lower"], v["Volume agg CI upper"], color=c, alpha=alpha1)
					plt.plot(v["Volume agg CI lower"].index, v["Volume agg CI lower"], linestyle="-", color=c, linewidth=1, alpha=alpha2)
					plt.plot(v["Volume agg CI upper"].index, v["Volume agg CI upper"], linestyle="-", color=c, linewidth=1, alpha=alpha2)

					# Calculate integral between lower and upper
					# BUGFIX 2023-10-10: Changed from v["Dose"] to v.index, since this is how the information was stored
					diffDose = v.index[1] - v.index[0]
					v["CI difference"] = (v["Volume agg CI upper"] - v["Volume agg CI lower"]) * diffDose
					area = v["CI difference"].sum()
					self.log(f"{self.dvhConfidenceInterval.get()} % CI area for {k}: {area:.3f}")

				if self.dvhStyleVarPercentile.get():
					qps = self.dvhPercentileInterval.get().split(",")
					for qp_raw in qps:
						qp = int(qp_raw) / 100
						qp_str = perc(qp*100)

						plt.plot(v[f"Volume agg {qp_str} percentile"].index, 
									v[f"Volume agg {qp_str} percentile"], 
									linestyle="--", 
									color=c, 
									linewidth=2, 
									label=f"{qp_str} percentile")
					
					if len(qps) > 1:
						qp_lower = int(qps[0]) / 100
						qp_lower_str = perc(qp_lower*100)
						qp_upper = int(qps[-1]) / 100
						qp_upper_str = perc(qp_upper*100)

						plt.fill_between(v.index, 
											  v[f"Volume agg {qp_lower_str} percentile"],
											  v[f"Volume agg {qp_upper_str} percentile"], 
											  color=c, 
											  alpha=alpha1)

		plt.xlabel("Dose [Gy]")
		plt.ylabel("Volume [%]")
		plt.legend()
		plt.show()

	elif self.dvhStyleVar2.get() == "compare":  # Across patients
		"""Use this to create a separate plot window for each structure in the cohort, with lines
				for mean / median values per plan for all patients."""


		plans_sorted = [k.split("/")[1] for k in cohortDVH.keys()]
		structures_sorted = [k.split("/")[0] for k in cohortDVH.keys()]
		
		plans = sorted(set(plans_sorted), key=plans_sorted.index)
		structures = sorted(set(structures_sorted), key=structures_sorted.index)

		styleIdx = {k:idx for idx,k in enumerate(plans)}
		stylesToUse = {
			'solid': "-", 
			'dashed': (0,(5,5)), 
			'dashdot': (0,(3,5,1,5,1,5)),
			'loosely dashed': (0, (5, 10)),
			'dotted': 'dotted', 
			'loosely dotted': (0, (1, 10)), 
			'loosely dashdotted': (0, (3, 10, 1, 10))
		}

		style = list(stylesToUse.values())
		
		plotsStructure = list()
		plotsPlan = list()
		figs = dict()

		if self.useCustomAggregateDVHPlot:
			fig, axs = plt.subplots(self.aggregateNrows.get(), self.aggregateNcols.get(), figsize=(13,10), squeeze=False)
			# axs[y][x]
		
		colors = ["red", "orange", "darkgreen", "blue", "lightsalmon", "darkviolet",
				"gold", "darkred", "indianred", "seagreen", "magenta", "goldenrod"]

		for k,v in cohortDVH.items():
				plan = k.split("/")[1]
				structure = k.split("/")[0]
				
				if self.dvhStyleSinglePlot.get():
					fignum = f"{self.dvhStyleVar1.get().capitalize()} DVH for all structures"
					if len(list(set(structures))) == 1:
						fignum = f"{self.dvhStyleVar1.get().capitalize()} DVH for {structure}"
					elif len(list(set(plans))) == 1:
						fignum = f"{self.dvhStyleVar1.get().capitalize()} DVH for {plan}"
					if len(list(set(structures))) == 1 and len(list(set(plans))) == 1:
						fignum = f"{self.dvhStyleVar1.get().capitalize()} DVH for {plan} / {structure}"

					planLabel = None
				else:
					fignum = structure
					planLabel = plan

				if not self.useCustomAggregateDVHPlot:
					plt.figure(figsize=(10,7.5), num = fignum)

				try:
					ls = style[styleIdx[plan]]
				except:
					ls = '-'

				c = self.colorVarList[structure].get()

				# if len(plans) > 1 and len(structures) == 1:
				# 	c = colors[styleIdx[plan]]

				if ":" in c:
					alpha1 = float(c.split(":")[-1])
					alpha2 = min(1, alpha1+0.6)
					c = c.split(":")[0]
				else:
					alpha1 = 0.3
					alpha2 = 0.7

				if self.dvhStyleBlackPlot.get():
					c = "k"
					
				if self.dvhStyleVar3.get():
					if not self.useCustomAggregateDVHPlot:
						plt.fill_between(v.index, v["Volume agg CI lower"], v["Volume agg CI upper"], color=c, alpha=alpha1)
						plt.plot(v["Volume agg CI lower"].index, v["Volume agg CI lower"], linestyle=ls, color=c, linewidth=2, alpha=alpha2)
						plt.plot(v["Volume agg CI upper"].index, v["Volume agg CI upper"], linestyle=ls, color=c, linewidth=2, alpha=alpha2)
					else:
						for x in range(self.aggregateNcols.get()):
								for y in range(self.aggregateNrows.get()):
									xy = f"{x}{y}"
									planMatches = self.aggregateGridOptions[xy]["planMatches"]
									structureMatches = self.aggregateGridOptions[xy]["structureMatches"]
									
									if plan in planMatches and structure in structureMatches:
										axs[y][x].fill_between(v.index,
																		v["Volume agg CI lower"], v["Volume agg CI upper"], color=c, alpha=alpha1)
										axs[y][x].plot(v["Volume agg CI lower"].index,
															v["Volume agg CI lower"], linestyle=ls, color=c, linewidth=1.5, alpha=alpha2)
										axs[y][x].plot(v["Volume agg CI upper"].index,
															v["Volume agg CI upper"], linestyle=ls, color=c, linewidth=1.5, alpha=alpha2)

					# Calculate integral between lower and upper
					diffDose = v.index[1] - v.index[0]
					v["CI difference"] = (v["Volume agg CI upper"] - v["Volume agg CI lower"]) * diffDose
					area = v["CI difference"].sum()
					self.log(f"{self.dvhConfidenceInterval.get()} % CI area for {k}: {area:.3f}")
					
					if not os.path.exists("Output/CI"):
						os.makedirs("Output/CI")
						
					v.to_excel(f"Output/CI/{plan}_{structure}.xlsx")

				if self.dvhStyleVarPercentile.get():
					qps = self.dvhPercentileInterval.get().split(",")
					for qp_raw in qps:
						qp = int(qp_raw) / 100
						qp_str = perc(qp*100)

						plt.plot(v[f"Volume agg {qp_str} percentile"].index, 
									v[f"Volume agg {qp_str} percentile"], 
									linestyle=ls,
									color=c, 
									linewidth=0.5, 
									label=f"{qp_str} percentile")
					
					if len(qps) > 1:
						qp_lower = int(qps[0]) / 100
						qp_lower_str = perc(qp_lower*100)
						qp_upper = int(qps[-1]) / 100
						qp_upper_str = perc(qp_upper*100)

						plt.fill_between(v.index, 
											  v[f"Volume agg {qp_lower_str} percentile"],
											  v[f"Volume agg {qp_upper_str} percentile"], 
											  color=c, 
											  alpha=alpha1)

					if not os.path.exists("Output/CI"):
						os.makedirs("Output/CI")

					v.to_excel(f"Output/CI/percentiles_{plan}_{structure}.xlsx")

				if not self.useCustomAggregateDVHPlot:                
					plt.plot(v["Volume agg"].index, 
								v["Volume agg"], 
								linestyle=ls, 
								color=c,
								linewidth=2, 
								label=planLabel, 
								alpha=1)

					plt.xlabel("Dose [Gy]", fontsize=12)
					plt.ylabel("Volume [%]", fontsize=12)
					plt_type = self.dvhStyleVar1.get().capitalize()

					if self.dvhStyleSinglePlot.get() and len(list(set(structures))) > 1:
						plt.title(f"{plt_type} DVH for all structures")
						if len(list(set(plans))) == 1:
								plt.title(f"{plt_type} DVH for {plan}")
					else:
						if len(list(set(plans))) == 1:
								plt.title(f"{plt_type} DVH for {plan} / {structure}")
						else:
								plt.title(f"{plt_type} DVH for {structure}")
				else:
					for x in range(self.aggregateNcols.get()):
						for y in range(self.aggregateNrows.get()):
							xy = f"{x}{y}"
							planMatches = self.aggregateGridOptions[xy]["planMatches"]
							structureMatches = self.aggregateGridOptions[xy]["structureMatches"]

							if plan in planMatches and structure in structureMatches:
								axs[y][x].plot(v["Volume agg"].index, v["Volume agg"],
													linestyle=ls, color=c, linewidth=2, label=planLabel, alpha=1)
					
								axs[y][x].set_xlabel("Dose [Gy]", fontsize=12)
								axs[y][x].set_ylabel("Volume [%]", fontsize=12)
								structureMatchesStr = ", ".join(structureMatches)
								plt_type = self.dvhStyleVar1.get().capitalize()
								axs[y][x].set_title(f"{plt_type} DVH for {structureMatchesStr}")
								
		cleaned_colorVarDict = { k:v.get().split(":")[0] for k,v in self.colorVarList.items() }

		if self.dvhStyleSinglePlot.get():
			if not self.useCustomAggregateDVHPlot:
				if self.dvhStyleBlackPlot.get():
					custom_lines = {k:Line2D([0], [0], color="k", ls=style[styleIdx[k]], lw=2) for k,v in zip(plans, colors)}
					custom_lines2 = {k:Line2D([0], [0], color="k", ls="-", lw=2) for k,v in cleaned_colorVarDict.items()}

				else:
					if len(list(set(structures))) == 1: # Plan lines @ structure colors if one structure, black if more

						custom_lines = {k:Line2D([0], [0], 
											color=self.colorVarList[structure].get(), 
											ls=style[styleIdx[k]], lw=2) \
											for k,v in zip(plans, colors)}

					else:
						custom_lines = {k:Line2D([0], [0], color="k", ls=style[styleIdx[k]], lw=2) \
									for k,v in zip(plans, colors)}

					custom_lines2 = {k:Line2D([0], [0], color=v, ls="-", lw=2) \
									for k,v in cleaned_colorVarDict.items()}

					if self.dvhStyleVarPercentile.get():
						qps = self.dvhPercentileInterval.get().split(",")

						if len(list(set(plans))) > 1:
							custom_lines[""] = Line2D([],[],linestyle='')
							for qp_raw in sorted(qps, reverse=True, key=int):
								qp_str = perc(int(qp_raw))
								custom_lines[f"{qp_str} percentile ({k})"] \
										= Line2D([0], [0], color="k", ls='-', lw=0.5)
						else:
							custom_lines2[""] = Line2D([],[],linestyle='')
							qp_raw = self.dvhPercentileInterval.get().split(",")
							for qp_raw in sorted(qps, reverse=True, key=int):
								qp_str = perc(int(qp_raw))
								custom_lines2[f"{qp_str} percentile ({k})"] \
											= Line2D([0], [0], color="k", ls='-', lw=0.5)

				all_legends = list(custom_lines.values()) + [Line2D([],[],linestyle='')] \
									+ list(custom_lines2.values())
				all_labels = list(custom_lines.keys()) + [''] + list(custom_lines2.keys())

				if len(list(set(structures))) == 1:
					all_legends = list(custom_lines.values())
					all_labels = list(custom_lines.keys())
				elif len(list(set(plans))) == 1:
					all_legends = list(custom_lines2.values())
					all_labels = list(custom_lines2.keys())

				# if len(list(set(structures))) > 1 or len(list(set(plans))) > 1 or self.dvhStyleVarPercentile.get():
				plt.legend(all_legends, all_labels, handlelength=3)

			else: # customAggregatePlots == True
				if self.dvhStyleBlackPlot.get():
					custom_lines = {k:Line2D([0], [0], color="k", ls=style[styleIdx[k]], lw=2) \
								for k,v in zip(plans, colors)}
					custom_lines2 = {k:Line2D([0], [0], color="k", ls="-", lw=2) \
								for k,v in cleaned_colorVarDict.items()}    
				else: # Use colors
					custom_lines = {k:Line2D([0], [0], color=v, ls=style[styleIdx[k]], lw=2) \
								for k,v in zip(plans, colors)}
					custom_lines2 = {k:Line2D([0], [0], color=v, ls="-", lw=2) \
								for k,v in cleaned_colorVarDict.items()}

				all_legends = list(custom_lines.values()) + [Line2D([],[],linestyle='')] \
								+ list(custom_lines2.values())
				all_labels = list(custom_lines.keys()) + [''] + list(custom_lines2.keys())
				if len(list(set(structures))) == 1:
					all_legends = list(custom_lines.values())
					all_labels = list(custom_lines.keys())
				elif len(list(set(plans))) == 1:
					all_legends = list(custom_lines2.values())
					all_labels = list(custom_lines2.keys())


				for x in range(self.aggregateNcols.get()):
					for y in range(self.aggregateNrows.get()):
							xy = f"{x}{y}"
							planMatches = self.aggregateGridOptions[xy]["planMatches"]
							structureMatches = self.aggregateGridOptions[xy]["structureMatches"]
							idxs = [ all_labels.index(k) for k in planMatches ] \
										+ [all_labels.index('')] + [ all_labels.index(k) for k in structureMatches ]
							axs[y][x].legend(itemgetter(*idxs)(all_legends), 
													itemgetter(*idxs)(all_labels), 
													handlelength=3)

		else: # self.dvhStyleSinglePlot.get() == False
			if self.dvhStyleBlackPlot.get():
				custom_lines = {k:Line2D([0], [0], color="k", ls=style[styleIdx[k]], lw=2) for k,v in zip(plans, colors)}

			else:
				# Loop through plot numbers
				for structure in structures_sorted:
					c = self.colorVarList[structure].get()
					custom_lines = {k:Line2D([0], [0], color=c, ls=style[styleIdx[k]], lw=2) \
								for k,v in zip(plans, colors)}

					if self.dvhStyleVarPercentile.get():
						custom_lines[""] = Line2D([],[],linestyle='')
							
						qps = self.dvhPercentileInterval.get().split(",")
						for qp_raw in sorted(qps, reverse=True, key=int):
							qp_str = perc(int(qp_raw))
							custom_lines[f"{qp_str} percentile"] \
									= Line2D([0], [0], color=c, ls='-', lw=0.5)
					
					all_legends = list(custom_lines.values())
					all_labels = list(custom_lines.keys())

					# Change to correct plot number
					plt.figure(structure)
					plt.legend(all_legends, all_labels, handlelength=3)
		plt.show()
				
	else: # Subtract two cohorts
		if ":" in c:
				alpha = c.split(":")[-1]
				c = c.split(":")[0]
				
		style = ['-', '--']
		fig = plt.figure(figsize=(6*1.5,8*1.5))
		colorSet = {'PTV72.5': 'darkred', 'PTV67.5': 'indianred', 'PTV50': 'red', 'PTV60': 'salmon',
		'Rectum': 'seagreen', 'Intestine': 'magenta', 'Bowel': 'magenta', 'Bladder': 'goldenrod' }
		
		custom_lines = {k:Line2D([0], [0], color="k", ls=k, lw=2) for k in style}
		if self.dvhStyleBlackPlot.get():
				custom_lines2 = {k:Line2D([0], [0], color="k", ls='-', lw=2) for k in colorSet.values()}
		else:
				custom_lines2 = {k:Line2D([0], [0], color=k, ls='-', lw=2) for k in colorSet.values()}

		plans_sorted = [k.split("/")[1] for k in cohortDVH.keys()]
		structures_sorted = [k.split("/")[0] for k in cohortDVH.keys()]
		
		plans = sorted(set(plans_sorted), key=plans_sorted.index)
		structures = sorted(set(structures_sorted), key=structures_sorted.index)

		cohortSum = {}
		for plan in plans:
				cohortSum[plan] = {}
				for cohort in cohortDVH.keys():
					if plan not in cohort: continue
					cohortSum[plan][cohort] = pd.DataFrame(cohortDVH[cohort].index \
															* cohortDVH[cohort]["Volume agg"])
					cohortSum[plan][cohort] = cohortSum[plan][cohort].sum().sum()

				cohortSum[plan] = np.sum(list(cohortSum[plan].values()))

		sorted_by_value = sorted(cohortSum.items(), key=lambda kv: kv[1], reverse=True)
		kLarge = sorted_by_value[1][0]
		kSmall = sorted_by_value[0][0]

		pd.set_option('display.max_columns', None)  

		cohortDiff = {}
		cohortDiffPerPatient = {}
		for structure in structures:
				kLargeCohort = f"{structure}/{kLarge}"
				kSmallCohort = f"{structure}/{kSmall}"
				rename_dict = {k1:k2 for k1,k2 in zip(*[cohortDVH[kLargeCohort].columns, 
																  cohortDVH[kSmallCohort].columns])}

				cohortDiff[structure] = cohortDVH[kLargeCohort]["Volume agg"] \
												- cohortDVH[kSmallCohort]["Volume agg"]
				cohortDiff[structure] = cohortDiff[structure].interpolate(method='linear', 
																				  limit_direction='backward', 
																				  limit=1).fillna(0)
				
				cohortDiffPerPatient[structure] = \
							cohortDVH[kLargeCohort].rename(columns=rename_dict) \
							- cohortDVH[kSmallCohort]

				_df = cohortDiffPerPatient[structure]

				_df.drop(["Volume agg", 
							"Volume agg CI lower", 
							"Volume agg CI upper"], 
							axis=1,									  
							inplace=True)

				if self.dvhStyleVarPercentile.get():
					qps = self.dvhPercentileInterval.get().split(",")
					for qp_raw in qps:
						qp_str = perc(int(qp_raw))
						_df.drop(["Volume agg {qp_str} percentile"],
																		 axis=1,
																		 inplace=True)

				if self.dvhStyleVar1.get() == "mean":
					_df["Volume agg"] = _df.mean(axis=1)
				else:
					_df["Volume agg"] = _df.median(axis=1)

				if self.dvhStyleVar3.get():
					ql = (1 - self.dvhConfidenceInterval.get()/100)/2
					qu = 1-ql
					_df["Volume agg CI lower"] = _df.quantile(ql, axis=1)
					_df["Volume agg CI upper"] = _df.quantile(qu, axis=1)
				
				if self.dvhStyleVarPercentile.get():
					qps = self.dvhPercentileInterval.get().split(",")
					for qp_raw in qps:
						qp =  int(qp_raw) / 100
						qp_str = perc(int(qp_raw))
						_df["Volume agg {qp_str} percentile"] = _df.quantile(qp, axis=1)

		plotsStructure = list()
		for structure in structures:
				c = colorSet[structure]
				if self.dvhStyleBlackPlot.get():
					c = "k"
				ls = '-'
				
				if self.dvhStyleVar2.get() == "subtract":
					cohortDiff[structure].plot(use_index=True, color=c, label=structure)
					plt_type = self.dvhStyleVar1.get().capitalize()
					plt.title(f"{plt_type} DVH {kLarge} - {plt_type} DVH {kSmall}")
				elif self.dvhStyleVar2.get() == "subtractPerPatient":
					_df["Volume agg"].plot(use_index=True, color=c, label=structure)

					if self.dvhStyleVar3.get():
						plt.fill_between(_df.index,
												_df["Volume agg CI lower"],
												_df["Volume agg CI upper"],
												color=structure=="Bladder" and "gold" or c,
												alpha=structure=="Bladder" and 0.5 or 0.3)
						
						plt.plot(_df["Volume agg CI lower"].index,
									_df["Volume agg CI lower"],
									linestyle=ls, color=c, linewidth=1, alpha=structure=="Bladder" and 0.9 or 0.7)
						
						plt.plot(_df["Volume agg CI upper"].index,
									_df["Volume agg CI upper"], linestyle=ls,
									color=c, linewidth=1, alpha=structure=="Bladder" and 0.9 or 0.7) 
						
					if self.dvhStyleVarPercentile.get():
						qps = self.dvhPercentileInterval.get().split(",")
						for qp_raw in qps:
							qp = int(qp_raw) / 100
							qp_str = perc(qp*100)

							plt.plot(_df[f"Volume agg {qp_str} percentile"].index,
										_df[f"Volume agg {qp_str} percentile"],
										linestyle="--", 
										color=c, 
										linewidth=2, 
										label="{qp_str} percentile")
						
						if len(qps) > 1:
							qp_lower = int(qps[0]) / 100
							qp_lower_str = perc(qp_lower*100)
							qp_upper = int(qps[-1]) / 100
							qp_upper_str = perc(qp_upper*100)
							
							plt.fill_between(_df[f"Volume agg {qp_str} percentile"].index, 
												  _df[f"Volume agg {qp_lower_str} percentile"],
												  _df[f"Volume agg {qp_upper_str} percentile"], 
												  color=c,
												  alpha=alpha1)
		
							plt.rcParams['font.size'] = '12'
		plt.tick_params(labelsize=12)
		plt.xlabel("Dose [Gy]", fontsize=12)
		plt.ylabel("Volume [%]", fontsize=12)
		plt.legend()
		plt.show()
