import os
import random
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from operator import itemgetter
from math import *
from tkinter import *
from ..Patients import *


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
            if not cohort in doses:
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
                if not cohort in tox:
                    tox[cohort] = []
                    notox[cohort] = []
                if patient.getTox() >= self.options.toxLimit.get():
                    print(f"Adding tox for {name}")
                    tox[cohort].append(f"Volume_{name}")
                else:
                    print(f"Adding notox for {name}")
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

        if self.dvhSaveAggDVH.get() == 1:
            if not os.path.exists("Output/aggDVH"):
                os.makedirs("Output/aggDVH")
            
            for cohort in cohortDVH.keys():
                st = self.dvhStyleVar1.get()
                cohortStructure, cohortPlan = cohort.split("/")
                q = self.dvhConfidenceInterval.get()
                filename = f"Output/aggDVH/{cohortPlan}_{cohortStructure}_{st}.csv"
                splice = cohortDVH[cohort][["Volume agg", "Volume agg tox", "Volume agg notox", "Volume agg tox CI lower", "Volume agg tox CI upper"]]
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
                    print(newDVH)
                    newDVH.set_index("Dose", inplace=True)
                    cohortDVH[cohort] = cohortDVH[cohort].merge(newDVH, how="outer", right_index=True, left_index=True)

                print(name, len(cohortDVH[cohort]))

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

        if self.dvhSaveAggDVH.get() == 1:
            if not os.path.exists("Output/aggDVH"):
                os.makedirs("Output/aggDVH")
            
            for cohort in cohortDVH.keys():
                st = self.dvhStyleVar1.get()
                q = self.dvhConfidenceInterval.get()
                cohortStructure, cohortPlan = cohort.split("/")
                q = self.dvhConfidenceInterval.get()
                filename = f"Output/aggDVH/{cohortPlan}_{cohortStructure}_{st}.csv"
                splice = cohortDVH[cohort][["Volume agg", "Volume agg CI lower", "Volume agg CI upper"]]
                splice.to_csv(filename, sep=";", decimal=".", index=True, header=True, na_rep="")

    else:  # Compare across patients
        for cohort_, patientsInCohort in self.patients.items():
            for name, patient in patientsInCohort.patients.items():
                plan = patient.getPlan()
                structure = patient.getStructure()
                cohort = f"{structure}/{plan}"

                if cohort not in cohortDVH:
                    cohortDVH[cohort] = pd.DataFrame({"Dose": doses[cohort]["Dose"]})
                    cohortDVH[cohort].set_index("Dose", inplace=True)

                print("cohortDVH[cohort]:", cohortDVH[cohort])

                newDVH = pd.DataFrame({"Dose": patient.dvh["Dose"], f"Volume_{name}": patient.dvh["Volume"]})
                print("newDVH", newDVH)
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

        if self.dvhSaveAggDVH.get() == 1:
            if not os.path.exists("Output/aggDVH"):
                os.makedirs("Output/aggDVH")
            
            for cohort in cohortDVH.keys():
                st = self.dvhStyleVar1.get()
                q = self.dvhConfidenceInterval.get()
                cohortStructure, cohortPlan = cohort.split("/")
                q = self.dvhConfidenceInterval.get()
                filename = f"Output/aggDVH/{cohortPlan}_{cohortStructure}_{st}.csv"
                splice = cohortDVH[cohort][["Volume agg", "Volume agg CI lower", "Volume agg CI upper"]]
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
                diffDose = v["Dose"][1] - v["Dose"][0]
                v["CI difference"] = (v["Volume agg CI upper"] - v["Volume agg CI lower"]) * diffDose
                area = v["CI difference"].sum()
                self.log(f"{self.dvhConfidenceInterval.get()} % CI area for {k}: {area:.3f}")

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
        print(styleIdx)
        stylesToUse = {'solid': "-", 'dotted': 'dotted', 'dashed': (0,(5,5)), 'dashdotted': (0,(3,5,1,5,1,5)),
                       'loosely dotted': (0, (1, 10)), 'loosely dashdotted': (0, (3, 10, 1, 10)), 'loosely dashed': (0, (5, 10))}
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

            ls = style[styleIdx[plan]]
            c = self.colorVarList[structure].get()

            if ":" in c:
                alpha1 = float(c.split(":")[-1])
                alpha2 = min(1, alpha1+0.6)
                c = c.split(":")[0]
            else:
                alpha1 = 0.3
                alpha2 = 0.7

            # c = colors.pop(0)
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

            if not self.useCustomAggregateDVHPlot:                
                plt.plot(v["Volume agg"].index, v["Volume agg"], linestyle=ls, color=c, linewidth=2, label=planLabel, alpha=1)
                plt.xlabel("Dose [Gy]", fontsize=12)
                plt.ylabel("Volume [%]", fontsize=12)
                if self.dvhStyleSinglePlot.get() and len(list(set(structures))) > 1:
                    plt.title(f"{self.dvhStyleVar1.get().capitalize()} DVH for all structures")
                    if len(list(set(plans))) == 1:
                        plt.title(f"{self.dvhStyleVar1.get().capitalize()} DVH for {plan}")
                else:
                    if len(list(set(plans))) == 1:
                        plt.title(f"{self.dvhStyleVar1.get().capitalize()} DVH for {plan} / {structure}")
                    else:
                        plt.title(f"{self.dvhStyleVar1.get().capitalize()} DVH for {structure}")
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
                            axs[y][x].set_title(f"{self.dvhStyleVar1.get().capitalize()} DVH for {structureMatchesStr}")
                        
            if not self.dvhStyleSinglePlot.get():
                plt.legend()

        cleaned_colorVarDict = { k:v.get().split(":")[0] for k,v in self.colorVarList.items() }

        if self.dvhStyleSinglePlot.get():
            if not self.useCustomAggregateDVHPlot:
                if self.dvhStyleBlackPlot.get():
                    custom_lines = {k:Line2D([0], [0], color="k", ls=style[styleIdx[k]], lw=2) for k,v in zip(plans, colors)}
                    custom_lines2 = {k:Line2D([0], [0], color="k", ls="-", lw=2) for k,v in cleaned_colorVarDict.items()}
                else:
                    if len(list(set(structures))) == 1: # Plan lines @ structure colors if one structure, black if more
                        custom_lines = {k:Line2D([0], [0], color=c, ls=style[styleIdx[k]], lw=2) for k,v in zip(plans, colors)}
                    else:
                        custom_lines = {k:Line2D([0], [0], color="k", ls=style[styleIdx[k]], lw=2) for k,v in zip(plans, colors)}
                    custom_lines2 = {k:Line2D([0], [0], color=v, ls="-", lw=2) for k,v in cleaned_colorVarDict.items()}
                all_legends = list(custom_lines.values()) + [Line2D([],[],linestyle='')] + list(custom_lines2.values())
                all_labels = list(custom_lines.keys()) + [''] + list(custom_lines2.keys())
                if len(list(set(structures))) == 1:
                    all_legends = list(custom_lines.values())
                    all_labels = list(custom_lines.keys())
                elif len(list(set(plans))) == 1:
                    all_legends = list(custom_lines2.values())
                    all_labels = list(custom_lines2.keys())

                if len(list(set(structures))) > 1 or len(list(set(plans))) > 1:
                    plt.legend(all_legends, all_labels, handlelength=3)

            else:
                if self.dvhStyleBlackPlot.get():
                    custom_lines = {k:Line2D([0], [0], color="k", ls=style[styleIdx[k]], lw=2) for k,v in zip(plans, colors)}
                    custom_lines2 = {k:Line2D([0], [0], color="k", ls="-", lw=2) for k,v in cleaned_colorVarDict.items()}    
                else:
                    custom_lines = {k:Line2D([0], [0], color=v, ls=style[styleIdx[k]], lw=2) for k,v in zip(plans, colors)}
                    custom_lines2 = {k:Line2D([0], [0], color=v, ls="-", lw=2) for k,v in cleaned_colorVarDict.items()}
                all_legends = list(custom_lines.values()) + [Line2D([],[],linestyle='')] + list(custom_lines2.values())
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
                        idxs = [ all_labels.index(k) for k in planMatches ] + [all_labels.index('')] + [ all_labels.index(k) for k in structureMatches ]
                        axs[y][x].legend(itemgetter(*idxs)(all_legends), itemgetter(*idxs)(all_labels), handlelength=3)

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
                if not plan in cohort: continue
                cohortSum[plan][cohort] = pd.DataFrame(cohortDVH[cohort].index * cohortDVH[cohort]["Volume agg"])
                cohortSum[plan][cohort] = cohortSum[plan][cohort].sum().sum()
                print(f"Cohort {cohort} has a sum of {cohortSum[plan][cohort]}")

            cohortSum[plan] = np.sum(list(cohortSum[plan].values()))
            print(f"Plan {plan} has a total sum of {cohortSum[plan]}")

        sorted_by_value = sorted(cohortSum.items(), key=lambda kv: kv[1], reverse=True)
        kLarge = sorted_by_value[1][0]
        kSmall = sorted_by_value[0][0]
        print(f"Highest-dose plan (across all structures): {kLarge}")

        pd.set_option('display.max_columns', None)  

        cohortDiff = {}
        cohortDiffPerPatient = {}
        for structure in structures:
            kLargeCohort = f"{structure}/{kLarge}"
            kSmallCohort = f"{structure}/{kSmall}"
            rename_dict = {k1:k2 for k1,k2 in zip(*[cohortDVH[kLargeCohort].columns, cohortDVH[kSmallCohort].columns])}
            cohortDiff[structure] = cohortDVH[kLargeCohort]["Volume agg"] - cohortDVH[kSmallCohort]["Volume agg"]
            cohortDiff[structure] = cohortDiff[structure].interpolate(method='linear', limit_direction='backward', limit=1).fillna(0)
            
            cohortDiffPerPatient[structure] = cohortDVH[kLargeCohort].rename(columns=rename_dict) - cohortDVH[kSmallCohort]
            cohortDiffPerPatient[structure].drop(["Volume agg", "Volume agg CI lower", "Volume agg CI upper"], axis=1, inplace=True)
            if self.dvhStyleVar1.get() == "mean":
                cohortDiffPerPatient[structure]["Volume agg"] = cohortDiffPerPatient[structure].mean(axis=1)
            else:
                cohortDiffPerPatient[structure]["Volume agg"] = cohortDiffPerPatient[structure].median(axis=1)

            if self.dvhStyleVar3.get():
                ql = (1 - self.dvhConfidenceInterval.get()/100)/2
                qu = 1-ql
                cohortDiffPerPatient[structure]["Volume agg CI lower"] = cohortDiffPerPatient[structure].quantile(ql, axis=1)
                cohortDiffPerPatient[structure]["Volume agg CI upper"] = cohortDiffPerPatient[structure].quantile(qu, axis=1)

        plotsStructure = list()
        for structure in structures:
            c = colorSet[structure]
            if self.dvhStyleBlackPlot.get():
                c = "k"
            ls = '-'
            
            if self.dvhStyleVar2.get() == "subtract":
                cohortDiff[structure].plot(use_index=True, color=c, label=structure)
                plt.title(f"{self.dvhStyleVar1.get().capitalize()} DVH {kLarge} - self.dvhStyleVar1.get().capitalize() DVH {kSmall}")
            elif self.dvhStyleVar2.get() == "subtractPerPatient":
                cohortDiffPerPatient[structure]["Volume agg"].plot(use_index=True, color=c, label=structure)
                if self.dvhStyleVar3.get():
                    plt.fill_between(cohortDiffPerPatient[structure].index,
                                     cohortDiffPerPatient[structure]["Volume agg CI lower"],
                                     cohortDiffPerPatient[structure]["Volume agg CI upper"],
                                     color=structure=="Bladder" and "gold" or c,
                                     alpha=structure=="Bladder" and 0.5 or 0.3)
                    
                    plt.plot(cohortDiffPerPatient[structure]["Volume agg CI lower"].index,
                             cohortDiffPerPatient[structure]["Volume agg CI lower"],
                             linestyle=ls, color=c, linewidth=1, alpha=structure=="Bladder" and 0.9 or 0.7)
                    
                    plt.plot(cohortDiffPerPatient[structure]["Volume agg CI upper"].index,
                             cohortDiffPerPatient[structure]["Volume agg CI upper"], linestyle=ls,
                             color=c, linewidth=1, alpha=structure=="Bladder" and 0.9 or 0.7)                    
            
        txt = "Figure 2: Population mean and 90% CI of the absolute difference in " \
              "relative volume as a function of dose between RP and CP for each patient."
        plt.figtext(0.5,0.01,txt,wrap=True,horizontalalignment='center',fontsize=12)


        plt.rcParams['font.size'] = '12'
        plt.tick_params(labelsize=12)
        plt.xlabel("Dose [Gy]", fontsize=12)
        plt.ylabel("RapidPlan Volume [%] - Clinical Plan Volume [%]", fontsize=12)
        plt.legend()
        plt.show()