import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
from functools import partial

from tkinter import *
from tkinter import ttk
from tkinter import filedialog

from ..Patients import *
from ._Tooltip import Tooltip

def myQuit(self):
    self.options.saveOptions() # to config file
    self.parent.destroy()
    plt.close("all")
    self.quit()


def log(self, text):
    self.logtext.insert(END, text + "\n")
    self.logtext.see(END)


def loadPatientsCommand(self):
    newdir = filedialog.askdirectory(title="Get output file directory", initialdir=self.options.dataFolder.get())
    if not newdir:
        self.log("No folder selected, aborting.")
        return
        
    self.options.dataFolder.set(newdir)
        
    cohortName = newdir.split("/")[-1]
    if cohortName in list(self.patients.keys()):
        maxCohort = -1
        for patientCohort in list(self.patients.keys()):
            if cohortName == patientCohort[:-3] and patientCohort[-3] == "_" and patientCohort[-2:] in [f"{k:02d}" for k in range(10)]:
                maxCohort = max(maxCohort, int(patientCohort[-2:]))
        
        if maxCohort < 0:
            cohortName += "_01"
        else:
            cohortName += f"_{maxCohort + 1:02d}"

        self.log(f"New cohort already in list of patient cohorts, renaming to {cohortName}.")
    
    self.log("Loading cohort %s." % cohortName)
    self.patients[cohortName] = Patients(self.options)
    self.patients[cohortName].setDataFolder(newdir)
    
    if self.options.DVHFileType.get() in ["ECLIPSE", "RayStation"]:
        plans = self.patients[cohortName].findPlans()
        structures = self.patients[cohortName].findStructures()

        self.options.structureToUse.set("")
        self.options.planToUse.set("")
        if len(structures)>1:
            self.window = Toplevel(self)
            sList = ", ".join(sorted(structures))
            Label(self.window, text=f"The following structures are available, choose one (* for wildcard, | for multiple):\n\n{sList}\n", wraplength=750).pack()
            e1 = Entry(self.window, textvariable=self.options.structureToUse, width=60)
            e1.pack()
            e1.focus()
        else:
            self.options.structureToUse.set(structures[0])
        
        if len(plans)>1:
            if not self.window:
                self.window = Toplevel(self)
            pList = ", ".join(sorted(plans))
            Label(self.window, text=f"\n\nThe following plans are available, choose one (* for wildcard, | for multiple):\n\n{pList}\n", wraplength=750).pack()
            e2 = Entry(self.window, textvariable=self.options.planToUse, width=60)
            e2.pack()
            if not len(structures): e2.focus()
        else:
            self.options.planToUse.set(plans[0])
            
        if len(structures)>1 or len(plans)>1:
            b = Button(self.window, text="OK", command=partial(self.chooseStructureCommand, cohortName))
            b.pack()
            self.window.bind('<Return>', lambda event=None: b.invoke())
        
        else:
            self.log("Using structure: %s\n" % (structures[0]))
            self.options.structureToUse.set(structures[0])
            self.log("Using plan: %s\n" % (plans[0]))
            self.options.plansToUse.set(plans[0])

        self.progress['maximum'] = len(plans) * len(structures) * self.patients[cohortName].findNPatients()
        
    else:
        self.options.planToUse.set("")
        self.options.structureToUse.set("")
        log = self.patients[cohortName].loadPatientsSIMPLE(self.progress)
        self.log("\n".join(log))
        self.progress['maximum'] = np.sum([len(k.patients) for k in self.patients.values()])
        self.addPatientCohort(cohortName, self.options.structureToUse.get(), self.options.planToUse.get(), 
                              self.patients[cohortName].getNPatients(), self.patients[cohortName].getNPatientsWithTox())


def addPatientCohort(self, cohortName, structureName, planName, nPatients, nPatientsWithTox):
    cohortStr = f"{cohortName}\n"
    cohortStr += planName and f"{planName}\n" or ""
    cohortStr += structureName and f"{structureName}\n" or ""
    cohortStr += f"{nPatients} patients\n({nPatientsWithTox} with tox)"

    self.cohortList[cohortName] = [Frame(self.middleLeftLowerContainer)]
    self.cohortList[cohortName][0].pack(anchor=N)
    self.cohortList[cohortName].append(Button(self.cohortList[cohortName][0], 
                                        text="DEL", command=partial(self.removePatientCohort, cohortName), width=7, height=3, bg='white'))
    self.cohortList[cohortName].append(Label(self.cohortList[cohortName][0], 
                                        text=cohortStr, bg='white', width=15))
    for i in range(1,3):
        self.cohortList[cohortName][i].pack(side=LEFT, pady=i and 0 or 5)
    
    self.buttonShowDVH['state'] = 'normal'
    self.buttonCalculateGEUD['state'] = 'normal'
    self.changeNamingButton['state'] = 'normal'
    
    hasGEUDs = True
    for patientsInCohort in list(self.patients.values()):
        hasGEUDs *= patientsInCohort.checkGEUDsplines()
    if len(self.patients) == 0:
        hasGEUDs = False
    
    if hasGEUDs:
        self.buttonCalculateNTCP['state'] = 'normal'
        self.buttonCalculateAUROC['state'] = 'normal'
        self.buttonShowGEUDvsN['state'] = 'normal'
        self.buttonCalculateDVH['state'] = 'normal'
        self.buttonAggregateDVH['state'] = 'normal'
        self.log("Found n-vs-gEUD files for all patients\n\t(and adjusted n parameter accordingly)")
        self.buttonCalculateGEUD['state'] = 'normal'
        self.buttonLKBuncert['state'] = 'normal'
    
    elif not hasGEUDs and self.options.NTCPcalculation.get() == "Logit":
        self.buttonCalculateNTCP['state'] = 'normal'
        self.buttonCalculateAUROC['state'] = 'normal'
        self.buttonCalculateDVH['state'] = 'normal'
        self.buttonAggregateDVH['state'] = 'normal'
        self.buttonCalculateGEUD['state'] = 'normal'
        self.buttonLKBuncert['state'] = 'normal'
    else:
        self.buttonCalculateNTCP['state'] = 'disabled'
    #    self.buttonLKBuncert['state'] = 'disabled'
        self.buttonCalculateAUROC['state'] = 'disabled'
        self.buttonCalculateDVH['state'] = 'normal'
        self.buttonAggregateDVH['state'] = 'normal'
        self.buttonCalculateGEUD['state'] = 'normal'


def removePatientCohort(self, cohortName):
    del self.patients[cohortName]
    for i in [2,1,0]:
        self.cohortList[cohortName][i].destroy()
    del self.cohortList[cohortName]
    self.log(f"Deleted cohort {cohortName}.")
    self.bestParameters = []
    if not self.patients:
        self.buttonShowDVH['state'] = 'disabled'
        self.buttonCalculateGEUD['state'] = 'disabled'
        self.changeNamingButton['state'] = 'disabled'


def chooseStructureCommand(self,cohortName):
    self.log(f"Choosing structure: {self.options.structureToUse.get()}\n")
    self.window.destroy()
    if self.options.DVHFileType.get() == "ECLIPSE":
        res = self.patients[cohortName].loadPatientsECLIPSE(self.progress)
    elif self.options.DVHFileType.get() == "RayStation":
        res = self.patients[cohortName].loadPatientsRayStation(self.progress)
    self.addPatientCohort(cohortName, self.options.structureToUse.get(), self.options.planToUse.get(),
                          self.patients[cohortName].getNPatients(), self.patients[cohortName].getNPatientsWithTox())
    self.log("\n".join(res))


def calculateGEUDWindow(self):
    self.window = Toplevel(self)
    self.window.title("Calculate gEUD")
    self.window.focus()
    self.styleContainer = Frame(self.window)
    self.styleContainer.pack(anchor=W)
    self.styleContainer2 = Frame(self.window)
    self.styleContainer2.pack(anchor=W)
    self.styleContainer3 = Frame(self.window)
    self.styleContainer3.pack(anchor=W)
    self.styleContainer4 = Frame(self.window)
    self.styleContainer4.pack(anchor=W)


    self.buttonContainer = Frame(self.window)
    self.buttonContainer.pack(anchor=W, fill=X, expand=1)

    Label(self.styleContainer, text="Lowest n value: ").pack(side=LEFT)
    Entry(self.styleContainer, textvariable=self.options.nFrom, width=5).pack(side=LEFT, anchor=W)

    Label(self.styleContainer2, text="Highest n value: ").pack(side=LEFT)
    Entry(self.styleContainer2, textvariable=self.options.nTo, width=5).pack(side=LEFT, anchor=W)

    Label(self.styleContainer3, text="Calculation grid size: ").pack(side=LEFT)
    Entry(self.styleContainer3, textvariable=self.options.nGrid, width=5).pack(side=LEFT, anchor=W)

    for text, mode in [["Linear n", True], ["Logarithmic n", False]]:
        Radiobutton(self.styleContainer4, text=text, value=mode, variable=self.options.nIsLinear).pack(side=LEFT, anchor=W)
    
    b1 = Button(self.buttonContainer, text="(Re)calculate gEUD look-up-table", command=self.calculateGEUDCommand, width=self.button_width)
    b1.pack(side=LEFT, anchor=W)

    b2 = Button(self.buttonContainer, text="Cancel", command=self.cancelCalculateGEUDCommand, width=self.button_width)
    b2.pack(side=LEFT, anchor=W)

    self.window.bind('<Return>', lambda event=None: b1.invoke())
    self.window.bind('<Escape>', lambda event=None: b2.invoke())


def calculateGEUDCommand(self):
    self.window.destroy()

    for k,v in list(self.patients.items()):
        self.log(f"Calculating gEUD values for cohort {k}...")
        v.createGEUDs(self.progress)
    self.log("Done.")

    self.buttonCalculateNTCP['state'] = 'normal'
    self.buttonCalculateAUROC['state'] = 'normal'
    self.buttonShowGEUDvsN['state'] = 'normal'
    self.buttonCalculateDVH['state'] = 'normal'

def cancelCalculateGEUDCommand(self):
    self.window.destroy()

def toxLimitChange(self):
    if self.cohortList:
        for cohortName, patientCohort in list(self.cohortList.items()):
            cohortStr = f"{cohortName}\n"
            cohortStr += self.options.planToUse.get() and f"{self.options.planToUse.get()}\n" or ""
            cohortStr += self.options.structureToUse.get() and f"{self.options.structureToUse.get()}\n" or ""
            cohortStr += f"{self.patients[cohortName].getNPatients()} patients\n"
            cohortStr += f"({self.patients[cohortName].getNPatientsWithTox()} with tox)"
            patientCohort[-1].config(text=cohortStr)


def toxFromFilenameCommand(self):
    if self.options.loadToxFromFilename.get() == 1:
        self.options.toxLimit.set(1)
        self.toxLimitChange()


def customDVHHeaderCommand(self):
    headerList = self.options.customDVHHeader.split(",")
    volumeInHeaderList = "Volume" in headerList
    doseInHeaderList = "Dose [cGy]" in headerList or "Dose [mGy]" in headerList or "Dose [Gy]" in headerList
    if volumeInHeaderList and doseInHeaderList:
        self.customDVHHeaderEntry.configure({'bg':'white'})
        return True
    else:
        return False
        self.customDVHHeaderEntry.configure({'bg':'coral1'})


def selectDVHFileTypeCommand(self):
    pass
    """
    if self.options.DVHFileType.get() in ["ECLIPSE", "RayStation"]:
        self.customDVHHeaderEntry['state'] = 'disabled'
    else:
        self.customDVHHeaderEntry['state'] = 'normal'
    """


def showDVHCommand(self):
    self.window = Toplevel(self)
    self.window.title("Plot DVH values")
    self.styleContainer = Frame(self.window)
    self.styleContainer.pack(anchor=W)
    self.styleContainer2 = Frame(self.window)
    self.styleContainer2.pack(anchor=W)
    self.styleContainer3 = Frame(self.window)
    self.styleContainer3.pack(anchor=W)
    self.styleContainer4 = Frame(self.window)
    self.styleContainer4.pack(anchor=W)
    self.styleContainer5 = Frame(self.window)
    self.styleContainer5.pack(anchor=W)
    self.styleContainer6 = Frame(self.window)
    self.styleContainer6.pack(anchor=W)
    self.buttonContainer = Frame(self.window)
    self.buttonContainer.pack(anchor=W, fill=X, expand=1)
    
    Label(self.styleContainer, text="Line color style: ").pack(anchor=W)
    for text, mode in [["Red (tox) / Black (noTox)", True], ["Grouped by property", False]]:
        Radiobutton(self.styleContainer, text=text, value=mode, variable=self.options.dvhPlotUseToxAsColor).pack(side=LEFT, anchor=W)
        
    Label(self.styleContainer2, text="Legend items: ").pack(anchor=W)
    for text, mode in [["None", "none"], ["Structure name", "structureName"],
                        ["Plan name", "planName"], ["Folder name", "folderName"], ["File name", "fileName"]]:
        Radiobutton(self.styleContainer2, text=text, value=mode, variable=self.options.dvhPlotLegendMarker).pack(side=LEFT, anchor=W)
        
    Label(self.styleContainer3, text="Legend font size: ").pack(side=LEFT, anchor=W)
    Entry(self.styleContainer3, textvariable=self.options.dvhPlotLegendSize, width=10).pack(side=LEFT)

    Label(self.styleContainer4, text="Line style grouping: ").pack(anchor=W)
    for text,mode in [["Plan", "plan"], ["Structure", "structure"], ["Folder", "folder"], ["None", "none"]]:
        Radiobutton(self.styleContainer4, text=text, value=mode, variable=self.options.dvhPlotLineStyleGrouping).pack(side=LEFT, anchor=W)
    
    Label(self.styleContainer5, text="Separate plot windows? ").pack(anchor=W)
    for text,mode in [["For plans", "plan"], ["For structures", "structure"], ["No", "none"]]:
        Radiobutton(self.styleContainer5, text=text, value=mode, variable=self.options.dvhPlotSeparatePlots).pack(side=LEFT, anchor=W)

    Checkbutton(self.styleContainer6, text="Save plots? (Outputfiles/DVHfigures/): ", variable=self.options.dvhPlotsSaveFigs).pack(anchor=W)
    
    b = Button(self.buttonContainer, text="Show DVHs", command=self.showDVHPlotCommand, width=self.button_width)
    b.pack(side=LEFT, anchor=W)
    b2 = Button(self.buttonContainer, text="Exit", command=self.cancelCalculateDVHvalues, width=self.button_width)
    b2.pack(side=LEFT, anchor=W)
    
    self.window.bind('<Return>', lambda event=None: b.invoke())
    self.window.bind('<Escape>', lambda event=None: b2.invoke())


def showDVHPlotCommand(self):
    styles = ["-", "--", "-.", ":", (0,(1,10)), (0,(1,1)), (0,(1,1)), (0,(5,1)), (0,(5,5))]

    plans = set()
    structures = set()
    folders = set()

    for patientsInCohort in list(self.patients.values()):
        for patient in list(patientsInCohort.patients.values()):
            plans.add(patient.getPlan())
            structures.add(patient.getStructure())
            folders.add(patient.getDataFolder().split("/")[-1])

    plotDict = dict()
    idx = 0
    if self.options.dvhPlotSeparatePlots.get() == "plan":
        for plan in plans:
            plotDict[plan] = idx
            plt.figure(idx, figsize=(8,5))
            plt.title(f"DVHs for {plan}")
            plt.xlabel("Dose [Gy]")
            plt.ylabel("Volume [%]")
            plt.ylim([0, 105])
            idx += 1
    elif self.options.dvhPlotSeparatePlots.get() == "structure":
        for structure in structures:
            plotDict[structure] = idx
            plt.figure(idx, figsize=(8,5))
            plt.title(f"DVHs for {structure}")
            plt.xlabel("Dose [Gy]")
            plt.ylabel("Volume [%]")
            plt.ylim([0, 105])
            idx += 1
    else:
        plotDict[0] = 1
        plt.figure(1, figsize=(8,5))
        plt.title("All DVHs")
        plt.xlabel("Dose [Gy]")
        plt.ylabel("Volume [%]")
        plt.ylim([0, 105])

    styleDict = dict()
    custom_lines = dict()
    idx = 0
    if self.options.dvhPlotLineStyleGrouping.get() == "plan":
        for plan in plans:
            styleDict[plan] = styles[idx]
            custom_lines[plan] = Line2D([0], [0], color="k", ls=styles[idx], lw=2)
            idx += 1
    elif self.options.dvhPlotLineStyleGrouping.get() == "structure":
        for structure in structures:
            styleDict[structure] = styles[idx]
            custom_lines[structure] = Line2D([0], [0], color="k", ls=styles[idx], lw=2)
            idx += 1
    elif self.options.dvhPlotLineStyleGrouping.get() == "folder":
        for folder in folders:
            styleDict[folder] = styles[idx]
            custom_lines[folder] = Line2D([0], [0], color="k", ls=styles[idx], lw=2)
            idx += 1
    
    for patientsInCohort in list(self.patients.values()):
        for patient in list(patientsInCohort.patients.values()):
            if self.options.dvhPlotLineStyleGrouping.get() == "plan":
                style = styleDict[patient.getPlan()]
            elif self.options.dvhPlotLineStyleGrouping.get() == "structure":
                style = styleDict[patient.getStructure()]
            elif self.options.dvhPlotLineStyleGrouping.get() == "folder":
                style = styleDict[patient.getDataFolder().split("/")[-1]]
            else:
                style = styles[0]

            if self.options.dvhPlotSeparatePlots.get() == "plan":
                plt.figure(plotDict[patient.getPlan()])
            elif self.options.dvhPlotSeparatePlots.get() == "structure":
                plt.figure(plotDict[patient.getStructure()])
            else:
                plt.figure(1)

            name = None
            if self.options.dvhPlotLineStyleGrouping.get() == "none":
                if self.options.dvhPlotLegendMarker.get() == "planName":
                    name = patient.getPlan()
                elif self.options.dvhPlotLegendMarker.get() == "structureName":
                    name = patient.getStructure()
                elif self.options.dvhPlotLegendMarker.get() == "folderName":
                    name = patient.getDataFolder().split("/")[-1]
                elif self.options.dvhPlotLegendMarker.get() == "fileName":
                    name = patient.getID().split("_")[0]
    
            if self.options.dvhPlotUseToxAsColor.get():
                plt.plot(patient.dvh["Dose"], patient.dvh["Volume"], linestyle=style, lw=1,
                         color=patient.getTox() >= self.options.toxLimit.get() and "red" or "black")
                custom_lines_tox = dict()
                custom_lines_tox["Toxicity"] = Line2D([0], [0], color="r", ls="-", lw=2)
                custom_lines_tox["No toxicity"] = Line2D([0], [0], color="k", ls="-", lw=2)
                plt.legend(custom_lines_tox.values(), custom_lines_tox.keys(),
                           prop={'size':float(self.options.dvhPlotLegendSize.get())})
            else:
                plt.plot(patient.dvh["Dose"], patient.dvh["Volume"], linestyle=style, lw=1, label=name)
    
    if self.options.dvhPlotLegendMarker.get() in ['planName', 'structureName', 'folderName','fileName'] \
                and not self.options.dvhPlotUseToxAsColor.get():
        if self.options.dvhPlotLineStyleGrouping.get() in ["plan", "structure", "folder"]:
            for plotIdx in plotDict.values():
                plt.figure(plotIdx)
                plt.legend(custom_lines.values(), custom_lines.keys(),
                           prop={'size':float(self.options.dvhPlotLegendSize.get())})
        else:                
            plt.legend(prop={'size':float(self.options.dvhPlotLegendSize.get())})

    if self.options.dvhPlotsSaveFigs.get():
        for plotName, plotIdx in plotDict.items():
            plt.figure(plotIdx)
            if not os.path.exists("Output/DVHfigures/"): 
                os.makedirs("Output/DVHfigures/")
            
            plt.savefig(f"Output/DVHfigures/{plotName}.pdf")

    plt.show()


def bootstrapCorrectionMethodCommand(self):
    if self.bestParameters:
        if self.options.bootstrapCorrectionMethod.get() == "none":
            self.log("Changing bootstrap correction method to None.")
            self.log("Best parameters: {} -> {}".format([float(f"{k:.2f}") for k in self.bestParameters], 
                                                         [float(f"{k:.2f}") for k in self.bestParametersNone]))
            self.log("Confidence Interval: {} -> {}".format([f"{k1:.2f} - {k2:.2f}" for k1,k2 in self.confidenceInterval],
                                                            [f"{k1:.2f} - {k2:.2f}" for k1,k2 in self.confidenceIntervalNone]))
            self.bestParameters = self.bestParametersNone
            self.confidenceInterval = self.confidenceIntervalNone
            
        elif self.options.bootstrapCorrectionMethod.get() == "mean":
            self.log("Changing bootstrap correction method to Mean.")
            self.log("Best parameters: {} -> {}".format([float(f"{k:.2f}") for k in self.bestParameters], 
                                                         [float(f"{k:.2f}") for k in self.bestParametersMean]))
            self.log("Confidence Interval: {} -> {}".format([f"{k1:.2f} - {k2:.2f}" for k1,k2 in self.confidenceInterval],
                                                            [f"{k1:.2f} - {k2:.2f}" for k1,k2 in self.confidenceIntervalMean]))
            self.bestParameters = self.bestParametersMean
            self.confidenceInterval = self.confidenceIntervalMean
            
        elif self.options.bootstrapCorrectionMethod.get() == "median":
            self.log("Changing bootstrap correction method to Median.")
            self.log("Best parameters: {} -> {}".format([float(f"{k:.2f}") for k in self.bestParameters], 
                                                         [float(f"{k:.2f}") for k in self.bestParametersMedian]))
            self.log("Confidence Interval: {} -> {}".format([f"{k1:.2f} - {k2:.2f}" for k1,k2 in self.confidenceInterval],
                                                            [f"{k1:.2f} - {k2:.2f}" for k1,k2 in self.confidenceIntervalMedian]))
            self.bestParameters = self.bestParametersMedian
            self.confidenceInterval = self.confidenceIntervalMedian


def calculateNTCPWindow(self, draw=True):
    self.window = Toplevel(self)
    self.window.title("Calculate NTCP model")
    self.window.focus()
    
    self.NTCPContainer = Frame(self.window)
    self.NTCPButtonContainer = Frame(self.window)

    self.NTCPContainer.pack(anchor=W)
    self.NTCPButtonContainer.pack(anchor=N)

    self.paramLambdaEntry = list()
    self.paramGammaEntry = list()
    self.paramNEntry = list()
    self.paramMEntry = list()
    self.paramTD50Entry = list()
    self.paramBHsizeLabel = list()
    self.paramBHsizeEntry = list()
    
    self.NTCPcalculationContainer = Frame(self.NTCPContainer)
    self.NTCPpercentCcContainer = Frame(self.NTCPContainer)
    self.NTCPpercentContainer = Frame(self.NTCPContainer)
    self.optimizationSchemeContainer = Frame(self.NTCPContainer)
    self.optimizationMetricContainer = Frame(self.NTCPContainer)
    self.matrixSizeContainer = Frame(self.NTCPContainer)
    self.basinHoppingsIterationsContainer = Frame(self.NTCPContainer)
    self.basinHoppingsTemperatureContainer = Frame(self.NTCPContainer)
    self.basinHoppingsNsizeContainer = Frame(self.NTCPContainer)
    self.basinHoppingsMsizeContainer = Frame(self.NTCPContainer)
    self.basinHoppingsTD50sizeContainer = Frame(self.NTCPContainer)
    self.basinHoppingsJumpLenghtsContainer = Frame(self.NTCPContainer)
    self.NTCPBoundWeightContainer = Frame(self.NTCPContainer)
    self.NTCPBoundLowerContainer = Frame(self.NTCPContainer)
    self.NTCPBoundUpperContainer = Frame(self.NTCPContainer)
    self.timeDependencyContainer = Frame(self.NTCPContainer)
    self.timeDependencyLambdaContainer = Frame(self.NTCPContainer)
    self.timeDependencyGammaContainer = Frame(self.NTCPContainer)

    Label(self.NTCPContainer, text="NTCP OPTIONS", font=("Helvetica", 12) ).pack(pady=(15,0))
    Frame(self.NTCPContainer, bg="grey", relief=SUNKEN).pack(fill=X,expand=1,anchor=W)
    
    self.NTCPcalculationContainer.pack(anchor=W)
    Label(self.NTCPcalculationContainer, text="NTCP calculation: ").pack(side=LEFT, anchor=W)
    for text, mode in [("D%/cc + logit", "Logit"), ("LKB", "LKB")]:
        Radiobutton(self.NTCPcalculationContainer, text=text, variable=self.options.NTCPcalculation, 
                    value=mode, command=self.NTCPcalculationCommand).pack(side=LEFT, anchor=W)
    Tooltip(self.NTCPcalculationContainer, text="The parameter limits for the LKB model are set below. "
            "For the D% + logit model, only the percentage value has to be set.", wraplength=self.wraplength)

    self.NTCPpercentCcContainer.pack(anchor=W)
    Checkbutton(self.basinHoppingsNsizeContainer, text="Use cc?", variable=self.options.useNTCPcc, command=self.switchToNTCPcc).pack(anchor=W)
    
    self.NTCPpercentContainer.pack(anchor=W)
    Label(self.NTCPpercentContainer, text="NTCP calculation Dose value: ").pack(side=LEFT, anchor=W)
    self.NTCPpercentLabel = Entry(self.NTCPpercentContainer, textvariable=self.options.NTCPcalculationDpercent, width=5)
    self.NTCPpercentLabel.pack(side=LEFT, anchor=W)
    self.NTCPpercentLabelPercentLabel = StringVar(value="%")
    self.NTCPpercentLabelPercent = Label(self.NTCPpercentContainer, textvariable=self.NTCPpercentLabelPercentLabel)
    self.NTCPpercentLabelPercent.pack(side=LEFT, anchor=W)
    Tooltip(self.NTCPcalculationContainer, text="The D% value from the DVHs to be used in the logit model.", wraplength=self.wraplength)

    self.basinHoppingsNsizeContainer.pack(anchor=W)
    self.paramNLabel = Label(self.basinHoppingsNsizeContainer, text="gEUD n parameter: ")
    self.paramNLabel.pack(side=LEFT, anchor=W)
    self.paramNEntry.append(Checkbutton(self.basinHoppingsNsizeContainer, text="Fixed", variable=self.options.fixN, command=self.switchNto))
    self.paramNEntry[0].pack(side=LEFT, anchor=W)
    for text,mode in [("", self.options.nFrom), ("→", self.options.nTo)]:
        Label(self.basinHoppingsNsizeContainer, text=text).pack(side=LEFT, anchor=W)
        self.paramNEntry.append(Entry(self.basinHoppingsNsizeContainer, textvariable=mode, width=5))
        self.paramNEntry[-1].pack(side=LEFT, anchor=W)
    Tooltip(self.basinHoppingsNsizeContainer, text="The boundaries of the n/a parameter. "
            "It will affect new gEUD spline calculations. If it is fixed, the lower value is applied", wraplength=self.wraplength)

    self.basinHoppingsMsizeContainer.pack(anchor=W)
    self.paramMLabel = Label(self.basinHoppingsMsizeContainer, text="LKB m parameter: ")
    self.paramMLabel.pack(side=LEFT, anchor=W)
    self.paramMEntry.append(Checkbutton(self.basinHoppingsMsizeContainer, text="Fixed", variable=self.options.fixM, command=self.switchMto))
    self.paramMEntry[0].pack(side=LEFT, anchor=W)
    for text,mode in [("", self.options.mFrom), ("→", self.options.mTo)]:
        Label(self.basinHoppingsMsizeContainer, text=text).pack(side=LEFT, anchor=W)
        self.paramMEntry.append(Entry(self.basinHoppingsMsizeContainer, textvariable=mode, width=5))
        self.paramMEntry[-1].pack(side=LEFT, anchor=W)
    Tooltip(self.basinHoppingsMsizeContainer, text="The boundaries of the m/b parameter. "
            "If it is fixed, the lower value is applied", wraplength=self.wraplength)
    
    self.basinHoppingsTD50sizeContainer.pack(anchor=W)
    self.paramTD50Label = Label(self.basinHoppingsTD50sizeContainer, text="LKB TD50 parameter [Gy]: ")
    self.paramTD50Label.pack(side=LEFT, anchor=W)
    self.paramTD50Entry.append(Checkbutton(self.basinHoppingsTD50sizeContainer, text="Fixed",
                                           variable=self.options.fixTD50, command=self.switchTD50to))
    self.paramTD50Entry[0].pack(side=LEFT, anchor=W)
    for text,mode in [("", self.options.TD50From), ("→", self.options.TD50To)]:
        Label(self.basinHoppingsTD50sizeContainer, text=text).pack(side=LEFT, anchor=W)
        self.paramTD50Entry.append(Entry(self.basinHoppingsTD50sizeContainer, textvariable=mode, width=5))
        self.paramTD50Entry[-1].pack(side=LEFT, anchor=W)
    Tooltip(self.basinHoppingsTD50sizeContainer, text="The boundaries of the TD50 parameter [Gy]. "
            "If it is fixed, the lower value is applied", wraplength=self.wraplength)
                        
    self.optimizationSchemeContainer.pack(anchor=W)
    Label(self.optimizationSchemeContainer, text="NTCP optimization scheme: ").pack(side=LEFT, anchor=W)
    for text, mode in [("Gradient Descent Optimization","GradientDescent"), ("Matrix Minimization", "MatrixMinimization")]:
        Radiobutton(self.optimizationSchemeContainer, text=text, variable=self.options.optimizationScheme, value=mode).pack(side=LEFT, anchor=W)
    Tooltip(self.optimizationSchemeContainer, text="The optimization method. Gradient Descent is a \"smart\" algorithm searching for "
            "optimal n,m,TD50 values towards the negative gradient of the error (between tox and LKB). Matrix minimization creates a "
            "3D grid of values to be tested, the cell with the lowest error containing the optimal values. The former method is the "
            "best in terms of computational time, however the latter is more intuitive.", wraplength=self.wraplength)

    self.optimizationMetricContainer.pack(anchor=W)
    Label(self.optimizationMetricContainer, text="NTCP optimization metric: ").pack(side=LEFT, anchor=W)
    for text, mode in [("Least Squares", "LS"), ("Log Likelihood", "LLH")]:
        Radiobutton(self.optimizationMetricContainer, text=text, variable=self.options.optimizationMetric, value=mode).pack(side=LEFT, anchor=W)
    Tooltip(self.optimizationMetricContainer, text="Which metric to be minimized: a least-squares method sum_i (tox_i - LKB_i)^2, or "
            "a log likelihood sum_tox ln NTCP + sum_notox ln (1-NTCP).", wraplength=self.wraplength)
    
    self.matrixSizeContainer.pack(anchor=W)
    Label(self.matrixSizeContainer, text="Matrix Size: ").pack(side=LEFT, anchor=W)
    Entry(self.matrixSizeContainer, textvariable=self.options.matrixSize, width=5).pack(side=LEFT, anchor=W)
    Tooltip(self.matrixSizeContainer, text="The linear size of the grid in the Matrix Minimization optimization scheme. "
            "The CPU requirements increase as the cube of this number.", wraplength=self.wraplength)
    
    self.basinHoppingsIterationsContainer.pack(anchor=W)
    Label(self.basinHoppingsIterationsContainer, text="G.D. # of \"basin hoppings\": ").pack(side=LEFT, anchor=W)
    Entry(self.basinHoppingsIterationsContainer, textvariable=self.options.basinHoppingIterations, width=5).pack(side=LEFT, anchor=W)
    Tooltip(self.basinHoppingsIterationsContainer, text="The number of times to perturb, through \"basin hoppings\", the optimal parameter set in order to check if "
            "the found values are local or global minima.", wraplength=self.wraplength)
    
    self.basinHoppingsTemperatureContainer.pack(anchor=W)
    Label(self.basinHoppingsTemperatureContainer, text="G.D. \"basin hopping\" temperature: ").pack(side=LEFT, anchor=W)
    Entry(self.basinHoppingsTemperatureContainer, textvariable=self.options.basinHoppingTemperature, width=5).pack(side=LEFT, anchor=W)
    Tooltip(self.basinHoppingsTemperatureContainer, text="The expected difference in error between the different local minima.",
            wraplength=self.wraplength)
    
    self.basinHoppingsJumpLenghtsContainer.pack(anchor=W)
    Label(self.basinHoppingsJumpLenghtsContainer, text="G.D. \"basin hopping\" jump lengths: ").pack(side=LEFT, anchor=W)
    for text, mode in [("n: ", self.options.basinHoppingNsize), ("m: ",self.options.basinHoppingMsize), ("TD50: ", self.options.basinHoppingTD50size)]:
        self.paramBHsizeLabel.append(Label(self.basinHoppingsJumpLenghtsContainer, text=text))
        self.paramBHsizeLabel[-1].pack(side=LEFT, anchor=W)
        self.paramBHsizeEntry.append(Entry(self.basinHoppingsJumpLenghtsContainer, textvariable=mode, width=5))
        self.paramBHsizeEntry[-1].pack(side=LEFT, anchor=W)

    Tooltip(self.basinHoppingsJumpLenghtsContainer, 
            text="For each \"basin hopping\", each parameter is perturbed as a random number within ± these values.", wraplength=self.wraplength)

    self.NTCPBoundWeightContainer.pack(anchor=W)
    Label(self.NTCPBoundWeightContainer, text="Boundary conditions: Weight ").pack(side=LEFT, anchor=W)
    Entry(self.NTCPBoundWeightContainer, textvariable=self.options.NTCPBoundWeight, width=5).pack(side=LEFT, anchor=W)
    Tooltip(self.NTCPBoundWeightContainer,
            text="Add artificial data points at Lower (no tox) and Upper (tox) positions to avoid too broad distributions."
            " The weight is the number of points to put at the Lower / Upper positions: NPoints = Dataset size * weight", wraplength=self.wraplength)

    self.NTCPBoundLowerContainer.pack(anchor=W)
    Label(self.NTCPBoundLowerContainer, text="Boundary conditions: Lower edge [Gy] ").pack(side=LEFT, anchor=W)
    Entry(self.NTCPBoundLowerContainer, textvariable=self.options.NTCPBoundLower, width=5).pack(side=LEFT, anchor=W)

    self.NTCPBoundUpperContainer.pack(anchor=W)
    Label(self.NTCPBoundUpperContainer, text="Boundary conditions: Upper edge [Gy] ").pack(side=LEFT, anchor=W)
    Entry(self.NTCPBoundUpperContainer, textvariable=self.options.NTCPBoundUpper, width=5).pack(side=LEFT, anchor=W)

    self.timeDependencyContainer.pack(anchor=W)
    Checkbutton(self.timeDependencyContainer, text="Use time dependent NTCP", variable=self.options.NTCPTimeDependent).pack(anchor=W, side=LEFT)

    Tooltip(self.timeDependencyContainer,
        text="See Niyazi et al (2000) (10.1016/j.radonc.2019.09.008). Correct for time of last follow-up using "
        "the Weibull function F(t) = exp ((-\lambda t)^{\gamma}). Then NTCP(t) = NTCP(1-F(t)). Need a CSV file "
        "in the parent folder (named as [data folder].csv) with patient name + NTCP status time in months. Also "
        "give the fixed paramaters lambda, gamma or a range if they should be fitted to the data.",
        wraplength=self.wraplength)
    
    self.timeDependencyLambdaContainer = Frame(self.NTCPContainer)
    self.timeDependencyLambdaContainer.pack(anchor=W)
    self.paramLambdaLabel = Label(self.timeDependencyLambdaContainer, text="NTCP(t) Lambda parameter: ")
    self.paramLambdaLabel.pack(side=LEFT, anchor=W)
    self.paramLambdaEntry.append(Checkbutton(self.timeDependencyLambdaContainer, text="Fixed", variable=self.options.fixLambda, command=self.switchNto))
    self.paramLambdaEntry[0].pack(side=LEFT, anchor=W)
    for text,mode in [("", self.options.lambdaFrom), ("→", self.options.lambdaTo)]:
        Label(self.timeDependencyLambdaContainer, text=text).pack(side=LEFT, anchor=W)
        self.paramLambdaEntry.append(Entry(self.timeDependencyLambdaContainer, textvariable=mode, width=5))
        self.paramLambdaEntry[-1].pack(side=LEFT, anchor=W)
    Tooltip(self.timeDependencyLambdaContainer, text="The boundaries of the gamma parameter. "
            "If it is fixed, the lower value is applied", wraplength=self.wraplength)

    self.timeDependencyGammaContainer = Frame(self.NTCPContainer)
    self.timeDependencyGammaContainer.pack(anchor=W)
    self.paramGammaLabel = Label(self.timeDependencyGammaContainer, text="NTCP(t) Gamma parameter: ")
    self.paramGammaLabel.pack(side=LEFT, anchor=W)
    self.paramGammaEntry.append(Checkbutton(self.timeDependencyGammaContainer, text="Fixed", variable=self.options.fixGamma, command=self.switchNto))
    self.paramGammaEntry[0].pack(side=LEFT, anchor=W)
    for text, mode in [("", self.options.gammaFrom), ("→", self.options.gammaTo)]:
        Label(self.timeDependencyGammaContainer, text=text).pack(side=LEFT, anchor=W)
        self.paramGammaEntry.append(Entry(self.timeDependencyGammaContainer, textvariable=mode, width=5))
        self.paramGammaEntry[-1].pack(side=LEFT, anchor=W)
    Tooltip(self.timeDependencyGammaContainer, text="The boundaries of the gamma parameter. "
            "If it is fixed, the lower value is applied", wraplength=self.wraplength)

    Frame(self.NTCPContainer, bg="grey", relief=SUNKEN).pack(fill=X, expand=1, anchor=W, pady=5)

    if self.options.NTCPcalculation.get() == "Logit":
        self.NTCPpercentLabel['state'] = 'normal'
        self.paramNLabel['text'] = "Logit a parameter: "
        self.paramMLabel['text'] = "Logit b parameter: "
        self.paramNEntry[0]['variable'] = self.options.fixA
        self.paramMEntry[0]['variable'] = self.options.fixB
        self.paramNEntry[1]['textvariable'] = self.options.aFrom
        self.paramMEntry[1]['textvariable'] = self.options.bFrom
        self.paramNEntry[2]['textvariable'] = self.options.aTo
        self.paramMEntry[2]['textvariable'] = self.options.bTo
        
        self.paramBHsizeLabel[0]['text'] = "a: "
        self.paramBHsizeLabel[1]['text'] = "b: "
        self.paramBHsizeEntry[0]['textvariable'] = self.options.basinHoppingAsize
        self.paramBHsizeEntry[1]['textvariable'] = self.options.basinHoppingBsize
        self.paramBHsizeEntry[2]['state'] = 'disabled'
    
        for k in self.paramTD50Entry:
            k['state'] = 'disabled'
    else:
        self.NTCPpercentLabel['state'] = 'disabled'

    b1 = Button(self.NTCPButtonContainer, text="NTCP model fit", command=self.calculateNTCPCommand, width=self.button_width)
    b2 = Button(self.NTCPButtonContainer, text="Cancel", command=self.calculateNTCPWindowCancel, width=self.button_width)
    b1.pack(side=LEFT, anchor=N, pady=10)
    b2.pack(side=LEFT, anchor=N, pady=10)

    self.window.bind('<Return>', lambda event=None: b1.invoke())
    self.window.bind('<Escape>', lambda event=None: b2.invoke())


def calculateNTCPWindowCancel(self):
    plt.close("all")
    self.window.destroy()


def switchToNTCPcc(self):
    if self.options.useNTCPcc.get():
        self.NTCPpercentLabelPercentLabel.set("cc")
    else:
        self.NTCPpercentLabelPercentLabel.set("%")


def calculateNTCPCommand(self, draw=True):
    self.window.destroy()
    cohortList = list(self.patients.values())
    
    primaryCohort = cohortList[0]
    secondaryCohorts = len(cohortList) > 1 and cohortList[1:] or {}
    
    for cohort in cohortList:
        cohort.options = self.options

    if not primaryCohort.didBootstrap: # Don't recalculate after bootstrap (= scrap pivot correction)
        if self.options.NTCPTimeDependent.get():
            for patients in cohortList:
                df = pd.read_csv(f"{patients.dataFolder}.csv", sep=None, engine="python")
                patients.NTCPTimeDict = dict(zip(df.Name, df.Time))
        
        if self.options.NTCPcalculation.get() == "Logit":
            for patients in cohortList:
                patients.calculateDpercent(self.options.NTCPcalculationDpercent.get())

        if self.options.optimizationScheme.get() == "GradientDescent":
            for patients in cohortList:
                self.log(f"Performing {self.options.optimizationScheme.get()} ({self.options.optimizationMetric.get()}) optimization on cohort: {patients.cohort}.")

                res = patients.doGradientOptimization(self.progress)
                
                patients.bestParameters = res.x
                self.log("\n".join([f"{k}={v}" for k,v in list(res.items())]))
            
        elif self.options.optimizationScheme.get() == "MatrixMinimization":
            for patients in cohortList:
                self.log("Calculating toxicity array with matrix size {self.options.matrixSize.get()}")
                res = patients.doMatrixMinimization(self.progress)
                
                patients.bestParameters = res.x
                self.log("\n".join([f"{k}={v}" for k,v in list(res.items())]))
        
    if draw:
        plt.figure(figsize=(10,6))
        self.style1 = ["darkred", "darkblue", "k"] * 100
        self.style2 = ["r", "b", "k"] * 100
        for patients in cohortList:
            patients.drawSigmoid(self.log, self.style1, self.style2)
        plt.show()

def NTCPcalculationCommand(self):
    if self.options.NTCPcalculation.get() == "Logit":
        self.NTCPpercentLabel['state'] = 'normal'
        self.buttonCalculateNTCP['state'] = 'normal'
        self.buttonCalculateAUROC['state'] = 'normal'
        self.buttonCalculateDVH['state'] = 'normal'

        self.paramNLabel['text'] = "Logit a parameter: "
        self.paramMLabel['text'] = "Logit b parameter: "
        self.paramNEntry[0]['variable'] = self.options.fixA
        self.paramMEntry[0]['variable'] = self.options.fixB
        self.paramNEntry[1]['textvariable'] = self.options.aFrom
        self.paramMEntry[1]['textvariable'] = self.options.bFrom
        self.paramNEntry[2]['textvariable'] = self.options.aTo
        self.paramMEntry[2]['textvariable'] = self.options.bTo
        for k in self.paramTD50Entry:
            k['state'] = 'disabled'
            
        self.paramBHsizeLabel[0]['text'] = "a: "
        self.paramBHsizeLabel[1]['text'] = "b: "
        self.paramBHsizeEntry[0]['textvariable'] = self.options.basinHoppingAsize
        self.paramBHsizeEntry[1]['textvariable'] = self.options.basinHoppingBsize
        self.paramBHsizeEntry[2]['state'] = 'disabled'
    else:
        self.NTCPpercentLabel['state'] = 'disabled'
        hasGEUDs = True
        for patientsInCohort in list(self.patients.values()):
            hasGEUDs *= patientsInCohort.checkGEUDsplines()
        
        for k in self.paramTD50Entry:
            k['state'] = 'normal'            
        
        self.paramNLabel['text'] = "LKB n parameter: "
        self.paramMLabel['text'] = "LKB m parameter: "
        self.paramNEntry[0]['variable'] = self.options.fixN
        self.paramMEntry[0]['variable'] = self.options.fixM
        self.paramTD50Entry[0]['variable'] = self.options.fixTD50
        self.paramNEntry[1]['textvariable'] = self.options.nFrom
        self.paramMEntry[1]['textvariable'] = self.options.mFrom
        self.paramTD50Entry[1]['textvariable'] = self.options.TD50From
        self.paramNEntry[2]['textvariable'] = self.options.nTo
        self.paramMEntry[2]['textvariable'] = self.options.mTo
        self.paramTD50Entry[2]['textvariable'] = self.options.TD50To
        
        self.paramBHsizeLabel[0]['text'] = "n: "
        self.paramBHsizeLabel[1]['text'] = "m: "
        self.paramBHsizeEntry[0]['textvariable'] = self.options.basinHoppingNsize
        self.paramBHsizeEntry[1]['textvariable'] = self.options.basinHoppingMsize
        self.paramBHsizeEntry[2]['state'] = 'normal'
        
        if hasGEUDs and len(self.patients):
            self.buttonCalculateNTCP['state'] = 'normal'
            self.buttonCalculateAUROC['state'] = 'normal'
            self.buttonShowGEUDvsN['state'] = 'normal'
            self.buttonCalculateDVH['state'] = 'normal'
            self.log("Found gEUD spline files for all patients.")
            self.buttonCalculateGEUD['state'] = 'disabled'
            self.buttonLKBuncert['state'] = 'normal'

        
        elif not hasGEUDs and self.options.NTCPcalculation.get() == "Logit":
            self.buttonCalculateNTCP['state'] = 'normal'
            self.buttonCalculateAUROC['state'] = 'normal'
            self.buttonCalculateDVH['state'] = 'normal'
            self.buttonCalculateGEUD['state'] = 'normal'
            self.buttonLKBuncert['state'] = 'normal'

        else:
            self.buttonCalculateNTCP['state'] = 'disabled'
            self.buttonCalculateAUROC['state'] = 'disabled'
            self.buttonCalculateDVH['state'] = 'normal'
            self.buttonCalculateGEUD['state'] = 'normal'
            # self.buttonLKBuncert['state'] = 'normal'


def calculateBootstrapWindow(self):
    self.window = Toplevel(self)
    self.window.title("Calculate NTCP model")
    self.window.focus()
    
    self.bootstrapContainer = Frame(self.window)
    self.bootstrapButtonContainer = Frame(self.window)

    self.bootstrapContainer.pack(anchor=W)
    self.bootstrapButtonContainer.pack(anchor=N)

    self.confidenceIntervalPercentContainer = Frame(self.bootstrapContainer)
    self.confidenceIntervalSchemeContainer = Frame(self.bootstrapContainer)
    self.confidenceIntervalIterationsContainer = Frame(self.bootstrapContainer)
    self.bootstrapCorrectionMethodContainer = Frame(self.bootstrapContainer)
    self.makeDifferentialCIContainer = Frame(self.bootstrapContainer)
    self.confidenceIntervalLikelihoodLimitContainer = Frame(self.bootstrapContainer)
    self.confidenceIntervalShowModelsContainer = Frame(self.bootstrapContainer)
    self.nIsLogOrLinContainer = Frame(self.bootstrapContainer)

    Label(self.bootstrapContainer, text="BOOTSTRAP OPTIONS", font=("Helvetica", 12)).pack(pady=(15, 0))
    Frame(self.bootstrapContainer, bg="grey", relief=SUNKEN).pack(fill=X, expand=1, anchor=W)

    self.confidenceIntervalPercentContainer.pack(anchor=W)
    Label(self.confidenceIntervalPercentContainer, text="Confidence Interval percentage: ").pack(side=LEFT, anchor=W)
    Entry(self.confidenceIntervalPercentContainer, textvariable=self.options.confidenceIntervalPercent, width=7).pack(side=LEFT, anchor=W)
    Label(self.confidenceIntervalPercentContainer, text="%").pack(side=LEFT, anchor=W)
    Tooltip(self.confidenceIntervalPercentContainer, text="1 sigma: 68.75, 2 sigma: 95.45, 3 sigma: 99.73. 2-way p<0.05 @ = 83.4", wraplength=self.wraplength)

    self.nIsLogOrLinContainer.pack(anchor=W)
    for text, mode in [["Linear n", True], ["Logarithmic n", False]]:
        Radiobutton(self.nIsLogOrLinContainer, text=text, value=mode, variable=self.options.nIsLinear).pack(side=LEFT, anchor=W)

    self.confidenceIntervalSchemeContainer.pack(anchor=W)
    Label(self.confidenceIntervalSchemeContainer, text="CI technique: ").pack(side=LEFT)
    for text, mode in [("Non-parametric BS", "NonParametricBootstrapping"), ("Parametric BS", "ParametricBootstrapping"), ("Profile likelihood", "ProfileLikelihood")]:
        Radiobutton(self.confidenceIntervalSchemeContainer, text=text, variable=self.options.confidenceIntervalMethod, value=mode).pack(side=LEFT, anchor=W)
    Tooltip(self.confidenceIntervalSchemeContainer, text="Choose the scheme for calculating the confidence interval: \n\n"
            "Non-parametric bootstrapping: Randomly choose patients, with replacement, from sample 1000-2000 times and calculate n,m,TD50 from sample "
            "The confidence interval is chosen from the distribution of n,m,TD50 values "
            "Assumptions: The cohort covers the different parameters, a robust method. "
            "My assumption: We choose the same number of patients. \n\n"
            "Parametric bootstrapping: A synthesised population is derived from the found n,m,TD50 "
            "For each patient, calculate the NTCP using the optimized n,m,TD50 values. "
            "Generate a random number rn [0,1] for each patient; rn < NTCP -> tox, rn > NTCP -> no tox "
            "Find n,m,TD50 from the synthesized population using these values, repeat 1000-2000 times "
            "Assumption: The cohort describes the real population. \n\n"
            "Profile likelihood method: Vary each parameter individually until the LLH is decreased by an amount equal to half the "
            "the critical value of chi2(1) distribution at the desired significance level", wraplength=self.wraplength)

    self.bootstrapCorrectionMethodContainer.pack(anchor=W)
    Label(self.bootstrapCorrectionMethodContainer, text="Bootstrapping subtraction correction method: ").pack(side=LEFT)
    for text, mode in [("None", "none"), ("Mean", "mean"), ("Median", "median")]:
        Radiobutton(self.bootstrapCorrectionMethodContainer, text=text, variable=self.options.bootstrapCorrectionMethod, value=mode,
                    command=self.bootstrapCorrectionMethodCommand).pack(side=LEFT, anchor=W)
    Tooltip(self.bootstrapCorrectionMethodContainer, text="Use pivotal confidence interval estimate (a bit tighter than regular CI). The CI is now calculated "
            "as [2p - p_U, 2p - p_L] insted of [p_L, p_U] (for choices mean / median). In addition, recalculate the best parameter as "
            "new = [2*old - mean/median of bootstrap distribution]. Both methods are valid: Mean is preferred for large datasets, "
            "Median is more stable for smaller datasets. Can be changed after BS calculation and before model display.")

    self.makeDifferentialCIContainer.pack(anchor=W)
    Label(self.makeDifferentialCIContainer, text="Calculate differential TD50 distribution: ").pack(side=LEFT)
    for text, mode in [("Yes", 1), ("No", 0)]:
        Radiobutton(self.makeDifferentialCIContainer, text=text, variable=self.options.makeDifferentialCI, value=mode).pack(side=LEFT, anchor=W)
    Tooltip(self.makeDifferentialCIContainer, text="1) Calculate the bootstrap distribution of a cohort (all patient groups in the left panel). 2) Remove "
            "the cohort, and add a new one. When the TD50 distribution is calculated, the prior distribution is subtracted from the new one a pairwise "
            "basis, and a differential mean / CI is displayed. Remember: For a double cohort-comparison the 83% CI should be separated to achieve p<0.05 "
            " if the variance is approx. equal (Payton et al., J Insect Sci, 2003).")

    self.confidenceIntervalShowModelsContainer.pack(anchor=W)
    Label(self.confidenceIntervalShowModelsContainer, text="Display each BS model within NTCP CI: ").pack(side=LEFT)
    for text, mode in [("Yes", 1), ("No", 0)]:
        Radiobutton(self.confidenceIntervalShowModelsContainer, text=text, variable=self.options.confidenceIntervalShowModels, value=mode).pack(side=LEFT, anchor=W)
    Tooltip(self.confidenceIntervalShowModelsContainer, text="Draw all the bootstrapped models on top of the shaded area containing the NTCP confidence interval. "
            "Only models having a log likelihood within the given % CI are drawn. They should converge to the original line for CI -> 0%!")

    self.confidenceIntervalIterationsContainer.pack(anchor=W)
    Label(self.confidenceIntervalIterationsContainer, text="Number of CI bootstrapping iterations: ").pack(side=LEFT, anchor=W)
    Entry(self.confidenceIntervalIterationsContainer, textvariable=self.options.confidenceIntervalIterations, width=5).pack(side=LEFT)
    Tooltip(self.confidenceIntervalIterationsContainer, text="Each iterations gives parameter set n,m,TD50 from a synthetic cohort. "
            "Usually 1000-2000 for good statistics, this step might be a bit time demanding (esp. for LKB).", wraplength=self.wraplength)

    self.confidenceIntervalLikelihoodLimitContainer.pack(anchor=W)
    Label(self.confidenceIntervalLikelihoodLimitContainer, text="Log Likelihood threshold: ").pack(side=LEFT, anchor=W)
    Entry(self.confidenceIntervalLikelihoodLimitContainer, textvariable=self.options.confidenceIntervalLikelihoodLimit, width=7).pack(side=LEFT, anchor=W)
    Tooltip(self.confidenceIntervalLikelihoodLimitContainer, text="This is to filter out \"bad\" or at-limit parameter sets. Check LLH histogram to find this value, "
            "typically -1 to 0 to remove parameters where the LLH is >2 too high and close to 0.", wraplength=self.wraplength)

    Frame(self.bootstrapContainer, bg="grey", relief=SUNKEN).pack(fill=X, expand=1, anchor=W, pady=5)

    b1 = Button(self.bootstrapButtonContainer, text="Run Bootstrap", command=self.calculateLKBuncert, width=self.button_width)
    b2 = Button(self.bootstrapButtonContainer, text="Cancel", command=self.calculateBootstrapWindowCancel, width=self.button_width)
    b1.pack(side=LEFT, anchor=N, pady=10)
    b2.pack(side=LEFT, anchor=N, pady=10)

    self.window.bind('<Return>', lambda event=None: b1.invoke())
    self.window.bind('<Escape>', lambda event=None: b2.invoke())


def calculateBootstrapWindowCancel(self):
    plt.close("all")
    self.window.destroy()


def calculateAUROCCommand(self):
    cohortList = list(self.patients.values())
    primaryCohort = cohortList[0]
    secondaryCohorts = len(cohortList) > 1 and cohortList[1:] or {}
    primaryCohort.drawAUROC(secondaryCohorts, self.progress)


def calculateDVHCommand(self):
    self.window = Toplevel(self)
    self.window.title("Calculate DVH values")
    self.fileContainer = Frame(self.window)
    self.fileContainer.pack(anchor=W)
    self.inputContainer = Frame(self.window)
    self.inputContainer.pack(anchor=W)
    self.radioContainer1 = Frame(self.inputContainer)
    self.radioContainer1.pack(anchor=W)
    self.radioContainer2 = Frame(self.inputContainer)
    self.radioContainer2.pack(anchor=W)
    self.radioContainer3 = Frame(self.inputContainer)
    self.radioContainer3.pack(anchor=W)
    #self.radioContainer4 = Frame(self.inputContainer)
    #self.radioContainer4.pack(anchor=W)
    self.ntcpContainer = Frame(self.inputContainer)
    self.ntcpContainer.pack(anchor=W)
    self.ntcpContainerType = Frame(self.inputContainer)
    self.ntcpContainerType.pack(anchor=W)
    self.ntcpContainerParamsLogit = Frame(self.inputContainer)
    self.ntcpContainerParamsLogit.pack(anchor=W)
    self.ntcpContainerParamsLKB = Frame(self.inputContainer)
    self.ntcpContainerParamsLKB.pack(anchor=W)
    self.calculateMeanDoseContainer = Frame(self.inputContainer)
    self.calculateMeanDoseContainer.pack(anchor=W)
    self.calculateGEUDContainer = Frame(self.inputContainer)
    self.calculateGEUDContainer.pack(anchor=W)
    
    self.buttonContainer = Frame(self.window)
    self.buttonContainer.pack(anchor=N, fill=X, expand=1)        

    self.dvhCheckVarDoseAtVolume = IntVar(value=0)
    self.dvhCheckVarVolumeAtDose = IntVar(value=0)
    self.dvhCheckVarVolumeAtRelDose = IntVar(value=0)
    self.dvhCheckVarIncludeNTCP = IntVar(value=0)
    self.dvhCheckVarIncludeGEUD = IntVar(value=0)
    self.dvhCheckVarCalculateMeanDose = IntVar(value=0)
     
    self.dvhEntryVar1 = StringVar(value=0)
    self.dvhEntryVar2 = StringVar(value=0)
    self.dvhEntryVar3 = StringVar(value=0)
    self.outputFileNameVar = StringVar(value="Output/dvhValues.xlsx")
    
    Label(self.fileContainer, text="Output file name: ").pack(side=LEFT)
    Entry(self.fileContainer, textvariable=self.outputFileNameVar, width=30).pack(side=LEFT, anchor=W)
    Tooltip(self.fileContainer, text="Filetype must be .csv (CSV file) or .xlsx (Excel file)", wraplength=self.wraplength)
    
    Checkbutton(self.radioContainer1, text="Find dose [Gy] at volumes [%] ", variable=self.dvhCheckVarDoseAtVolume).pack(anchor=W, side=LEFT)
    Entry(self.radioContainer1, textvariable=self.dvhEntryVar1).pack(anchor=W, side=LEFT)
    Tooltip(self.radioContainer1, text="Calculate dose at multiple volumes by giving a comma-separated list ( e.g. 2, 5, 10, 50, 98)",
            wraplength=self.wraplength)
    
    Checkbutton(self.radioContainer2, text="Find volume [%] at doses [Gy] ", variable=self.dvhCheckVarVolumeAtDose).pack(anchor=W, side=LEFT)
    Entry(self.radioContainer2, textvariable=self.dvhEntryVar2).pack(anchor=W, side=LEFT)
    Tooltip(self.radioContainer2, text="Calculate the volume for multiple dose values by giving a comma-separated list ( e.g. 10, 20, 30)",
            wraplength=self.wraplength)

    Checkbutton(self.radioContainer3, text="Find volume [%] at doses [% of max] ", variable=self.dvhCheckVarVolumeAtRelDose).pack(anchor=W, side=LEFT)
    Entry(self.radioContainer3, textvariable=self.dvhEntryVar3).pack(anchor=W, side=LEFT)
    Tooltip(self.radioContainer3, text="Calculate the volume for multiple dose values by giving a comma-separated list ( e.g. 10, 20, 30)",
            wraplength=self.wraplength)

    Checkbutton(self.ntcpContainer, text="Calculate NTCP", variable=self.dvhCheckVarIncludeNTCP).pack(anchor=W)
    Tooltip(self.ntcpContainer, text="Calculate the NTCP values (0->1) for all patients. If a model fit has been performed ,"
           f" the resulting parameters will be applied here ({self.bestParameters}). If not, the parameters must be specified.", 
           wraplength=self.wraplength)

    Label(self.ntcpContainerType, text="NTCP calculation: ").pack(side=LEFT, anchor=W)
    for text, mode in [("D%/cc + logit", "Logit"), ("LKB", "LKB")]:
        Radiobutton(self.ntcpContainerType, text=text, variable=self.options.NTCPcalculation, 
                    value=mode).pack(side=LEFT, anchor=W)
    Tooltip(self.ntcpContainerType, text="The parameter limits for the LKB model are set below. "
            "For the D% + logit model, only the percentage value has to be set.", wraplength=self.wraplength)

    Label(self.ntcpContainerParamsLogit, text="Logit/a = ").pack(side=LEFT, anchor=W)
    Entry(self.ntcpContainerParamsLogit, textvariable=self.options.aFrom).pack(side=LEFT, anchor=W)
    Label(self.ntcpContainerParamsLogit, text="Logit/b = ").pack(side=LEFT, anchor=W)
    Entry(self.ntcpContainerParamsLogit, textvariable=self.options.bFrom).pack(side=LEFT, anchor=W)    
    Label(self.ntcpContainerParamsLKB, text="LKB/n = ").pack(side=LEFT, anchor=W)
    Entry(self.ntcpContainerParamsLKB, textvariable=self.options.nFrom).pack(side=LEFT, anchor=W)
    Label(self.ntcpContainerParamsLKB, text="LKB/m = ").pack(side=LEFT, anchor=W)
    Entry(self.ntcpContainerParamsLKB, textvariable=self.options.mFrom).pack(side=LEFT, anchor=W)        
    Label(self.ntcpContainerParamsLKB, text="LKB/TD50 = ").pack(side=LEFT, anchor=W)
    Entry(self.ntcpContainerParamsLKB, textvariable=self.options.TD50From).pack(side=LEFT, anchor=W)        

    Checkbutton(self.calculateMeanDoseContainer, text="Include dose metrics from ECLIPSE ", variable=self.dvhCheckVarCalculateMeanDose).pack(anchor=W)
    Checkbutton(self.calculateGEUDContainer, text="Calculate gEUD ", variable=self.dvhCheckVarIncludeGEUD).pack(anchor=W, side=LEFT)
    Label(self.calculateGEUDContainer, text="n = ").pack(anchor=W, side=LEFT)
    Entry(self.calculateGEUDContainer, textvariable=self.options.nSet, width=10).pack(anchor=W, side=LEFT)

    Tooltip(self.calculateGEUDContainer, text="Calculate the gEUD based on the pre-set fixed N value. Set n=1 for the mean dose.",
            wraplength=self.wraplength)


    b = Button(self.buttonContainer, text="Calculate", command=self.calculateDVHvalues, width=self.button_width)
    b.pack(side=LEFT, anchor=N)
    b2 = Button(self.buttonContainer, text="Cancel", command=self.cancelCalculateDVHvalues, width=self.button_width)
    b2.pack(side=LEFT, anchor=N)

    self.window.bind('<Return>', lambda event=None: b.invoke())
    self.window.bind('<Escape>', lambda event=None: b2.invoke())


def cancelCalculateDVHvalues(self):
    plt.close("all")
    self.window.destroy()


def aggregateDVHCommand(self):
    # Mean / median
    # Subtract or show all?
    
    self.window = Toplevel(self)
    self.window.title("Aggregate DVH plots")
    self.window.focus()
    self.styleContainer = Frame(self.window)
    self.styleContainer.pack(anchor=W)
    self.styleContainer2 = Frame(self.window)
    self.styleContainer2.pack(anchor=W)
    self.styleContainer3 = Frame(self.window)
    self.styleContainer3.pack(anchor=W)
    self.styleContainer35 = Frame(self.window)
    self.styleContainer35.pack(anchor=W)
    self.styleContainer4 = Frame(self.window)
    self.styleContainer4.pack(anchor=W)
    self.styleContainer5 = Frame(self.window)
    self.styleContainer5.pack(anchor=W, fill=X, expand=1)
    self.styleContainer6 = Frame(self.window)
    self.styleContainer6.pack(anchor=W, fill=X, expand=1)

    self.buttonContainer = Frame(self.window)
    self.buttonContainer.pack(anchor=W, fill=X, expand=1)
    
    self.dvhStyleVar1 = StringVar(value="median")
    self.dvhStyleVar2 = StringVar(value="compare")
    self.dvhStyleVar3 = IntVar(value=0)
    self.dvhSaveAggDVH = IntVar(value=0)
    self.dvhStyleSinglePlot = IntVar(value=1)
    self.dvhStyleBlackPlot = IntVar(value=0)
    self.dvhConfidenceInterval = DoubleVar(value=95)

    colorVarList = list()

    structures = set()
    for patientsInCohort in self.patients.values():
        for patient in patientsInCohort.patients.values():
            structures.add(patient.getStructure())
    structures = sorted(list(structures))
    
    colors = ["red", "orange", "darkgreen", "blue", "lightsalmon", "darkviolet",
              "gold", "darkred", "indianred", "seagreen", "magenta", "goldenrod"]
    self.colorVarList = { s : StringVar(value=k) for s,k in zip(structures,colors) }

    self.colorContainer = [ Frame(self.styleContainer5) for k in range(len(structures)//3+1) ]
    for k in self.colorContainer: k.pack(anchor=W)
    
    if len(colors) < len(structures):
        raise IndexError("Please specify more colors in MainMenu/_Actions.py::aggregateDVHCommand()")
                            
    Label(self.styleContainer, text="Cohort aggregation style: ").pack(side=LEFT)
    for text, mode in [["Median", "median"], ["Mean", "mean"]]:
        Radiobutton(self.styleContainer, text=text, value=mode, variable=self.dvhStyleVar1).pack(side=LEFT, anchor=W)
        
    for text, mode in [["Compare plans per structure (agg patients)", "compare"], ["Aggregate plans within patient", "comparePlans"], ["Compare tox vs no tox", "showAll"], ["Subtract [mean/median of all px]", "subtract"],
                        ["Mean/median [subtract per patient]", "subtractPerPatient"]]:
        Radiobutton(self.styleContainer2, text=text, value=mode, variable=self.dvhStyleVar2).pack(anchor=W)
    
    Label(self.styleContainer3, text="Draw ").pack(side=LEFT, anchor=W)
    Entry(self.styleContainer3, textvariable=self.dvhConfidenceInterval, width=4).pack(side=LEFT, anchor=W)
    Label(self.styleContainer3, text="% confidence interval: ").pack(side=LEFT, anchor=W)
    for text, mode in [["Yes", 1], ["No",0]]:
        Radiobutton(self.styleContainer3, text=text, value=mode, variable=self.dvhStyleVar3).pack(side=LEFT, anchor=W)

    Checkbutton(self.styleContainer35, text="Black/white plot? ", variable=self.dvhStyleBlackPlot).pack(anchor=W)
    Checkbutton(self.styleContainer4, text="Single plot window? ", variable=self.dvhStyleSinglePlot).pack(anchor=W)
    
    Label(self.colorContainer[0], text="Color scheme: ").pack(side=LEFT)
    for nColors, s in enumerate(structures):
       Label(self.colorContainer[(nColors+1)//3], text=f"{s}: ").pack(side=LEFT)
       Entry(self.colorContainer[(nColors+1)//3], textvariable=self.colorVarList[s], width=15).pack(side=LEFT)
    cStr2 = ", ".join(mcolors.CSS4_COLORS)
    Tooltip(self.styleContainer5, text=f"The available colors are (use color:alpha to set opacity of the {self.dvhConfidenceInterval.get()}% confidence envelope): \n\n {cStr2}.",
            wraplength=self.wraplength)

    Checkbutton(self.styleContainer6, text="Save DVH (-> Output/aggDVH)", variable=self.dvhSaveAggDVH).pack(anchor=W)
    
    b2 = Button(self.buttonContainer, text="Show aggregated DVH", command=self.calculateAggregatedDVH, width=self.button_width)
    b2.pack(side=LEFT, anchor=W)
    self.window.bind('<Return>', lambda event=None: b2.invoke())
    Button(self.buttonContainer, text="Custom plot placement", command=self.customAggregatedDVHCommand, width=self.button_width).pack(side=LEFT)
    b3 = Button(self.buttonContainer, text="Cancel", command=self.cancelCalculateDVHvalues, width=self.button_width)
    b3.pack(side=LEFT, anchor=W)
    self.window.bind('<Escape>', lambda event=None: b3.invoke())


def customAggregatedDVHCommand(self):
    self.aggWindow = Toplevel(self)
    self.aggWindow.title("Setup aggregated plot")
    self.aggWindow.focus()

    self.dropdownContainer = Frame(self.aggWindow)
    self.dropdownContainer.pack(anchor=N)

    ## GRID SETTINGS BASED ON INPUT DROPDOWN MENUS
    self.aggregateNrows = IntVar(value=1)
    self.aggregateNcols = IntVar(value=1)

    self.aggregateGridContainers = dict()
    self.aggregateGridOptions = dict()
    self.aggregateGridMatches = dict()
    self.aggregateYlines = list()

    Label(self.dropdownContainer, text="Number of plot rows: ").pack(side=LEFT)
    OptionMenu(self.dropdownContainer, self.aggregateNrows, *range(1,6), command=self.packCustomAggregateDVHCommand).pack(side=LEFT)
    Label(self.dropdownContainer, text="Number of plot cols: ").pack(side=LEFT)
    OptionMenu(self.dropdownContainer, self.aggregateNcols, *range(1,6), command=self.packCustomAggregateDVHCommand).pack(side=LEFT)

    Frame(self.aggWindow, bg="grey", relief=SUNKEN).pack(anchor=W, fill=X, expand=1, padx=5, pady=5)

    self.aggregateGridContainer = Frame(self.aggWindow)
    self.aggregateGridContainer.pack(anchor=W, fill=X, expand=1)

    self.aggregateButtonContainer = Frame(self.aggWindow)
    self.aggregateButtonContainer.pack(anchor=W, fill=X, expand=1)

    b2 = Button(self.aggregateButtonContainer, text="Find matches", command=self.matchCustomAggregateDVHCommand, width=self.button_width)
    b2.pack(side=LEFT, anchor=S)
    self.aggWindow.bind("<Return>", lambda event=None: b2.invoke())

    b3 = Button(self.aggregateButtonContainer, text="Save & Close", command=self.saveCustomAggregateDVHCommand, width=self.button_width)
    b3.pack(side=LEFT, anchor=S)

    b1 = Button(self.aggregateButtonContainer, text="Cancel", command=self.cancelCustomAggregateDVHCommand, width=self.button_width)
    b1.pack(side=LEFT, anchor=S)
    self.aggWindow.bind("<Escape>", lambda event=None: b1.invoke())

    self.packCustomAggregateDVHCommand(1)


def packCustomAggregateDVHCommand(self, value):
    for k,v in self.aggregateGridContainers.items():
        v.pack_forget()

    for v in self.aggregateYlines:
        v.pack_forget()
        
    self.aggregateGridContainers = dict()
    self.aggregateGridOptions = dict()
    self.aggregateGridMatches = dict()
    
    for y in range(self.aggregateNrows.get()):
        if y>0:
            self.aggregateYlines.append(Frame(self.aggregateGridContainer, bg="grey", relief=SUNKEN))
            self.aggregateYlines[-1].pack(anchor=W, fill=X, expand=1, padx=5, pady=5)
        self.aggregateGridContainers[y] = Frame(self.aggregateGridContainer)
        self.aggregateGridContainers[y].pack(anchor=W, fill=X, expand=1)
        for x in range(self.aggregateNcols.get()):
            xy = f"{x}{y}"
            if x>0:
                Frame(self.aggregateGridContainers[y], bg="grey", relief=SUNKEN).pack(anchor=W, side=LEFT, fill=Y, expand=1, padx=5)
            self.aggregateGridContainers[xy] = Frame(self.aggregateGridContainers[y])
            self.aggregateGridContainers[xy].pack(anchor=W, side=LEFT, fill=X, expand=1)
            self.aggregateGridOptions[xy] = { "structureRegex" : StringVar(value=""), "planRegex" : StringVar(value=""),
                                                 "structureMatches" : set(), "planMatches" : set(),
                                                 "structureMatchString" : StringVar(value=""), "planMatchString" : StringVar(value="")}

            Label(self.aggregateGridContainers[xy], text=f"PLOT X={x}, Y={y}").pack(anchor=W)
            
            structureContainer = Frame(self.aggregateGridContainers[xy])
            structureContainer.pack(anchor=W)
            Label(structureContainer, text="Match structures: ").pack(anchor=W, side=LEFT)
            Entry(structureContainer, textvariable=self.aggregateGridOptions[xy]["structureRegex"], width=25).pack(anchor=W, side=LEFT)
            Tooltip(structureContainer, text=f"Match using * (wildcard) and | (multiple)", wraplength=self.wraplength)
            Label(structureContainer, textvariable=self.aggregateGridOptions[xy]["structureMatchString"]).pack(side=LEFT, anchor=W)
            
            planContainer = Frame(self.aggregateGridContainers[xy])
            planContainer.pack(anchor=W)
            Label(planContainer, text="Match plans: ").pack(anchor=W, side=LEFT)
            Entry(planContainer, textvariable=self.aggregateGridOptions[xy]["planRegex"], width=25).pack(anchor=W, side=LEFT)
            Tooltip(planContainer, text=f"Match using * (wildcard) and | (multiple)", wraplength=self.wraplength)
            Label(planContainer, textvariable=self.aggregateGridOptions[xy]["planMatchString"]).pack(side=LEFT, anchor=W)


def cancelCustomAggregateDVHCommand(self):
    self.aggWindow.destroy()
    self.useCustomAggregateDVHPlot = False


def matchCustomAggregateDVHCommand(self):
    def match(a,b):
        a += "$" # Don't match end-of-lines
        a = a.replace("*", ".*") # Use regex type wildcard
        return re.match(a, b)

    for k,v in self.aggregateGridOptions.items():
        v["structureMatches"].clear()
        v["planMatches"].clear()
    
    for cohort_, patientsInCohort in self.patients.items():
        for name, patient in patientsInCohort.patients.items():
            plan = patient.getPlan()
            structure = patient.getStructure()

            for k,v in self.aggregateGridOptions.items():
                if match(v["structureRegex"].get(), structure): v["structureMatches"].add(structure)
                if match(v["planRegex"].get(), plan): v["planMatches"].add(plan)

    for k,v in self.aggregateGridOptions.items():
        v["structureMatchString"].set(", ".join(list(v["structureMatches"])))
        v["planMatchString"].set(", ".join(list(v["planMatches"])))


def saveCustomAggregateDVHCommand(self):
    self.aggWindow.destroy()
    self.useCustomAggregateDVHPlot = True


def changeNamingCommand(self):            
    self.window = Toplevel(self)
    self.window.title("Change plan / structure naming")

    self.structuresBefore = set()
    self.plansBefore = set()
    self.structuresAfter = set()
    self.plansAfter = set()

    for patientsInCohort in list(self.patients.values()):
        for patient in list(patientsInCohort.patients.values()):
            self.plansBefore.add(patient.getPlan())
            self.structuresBefore.add(patient.getStructure())

    self.planSubstituteList = list()
    self.structureSubstituteList = list()
    self.planEntries = list()
    self.structureEntries = list()

    self.planContainer = Frame(self.window)
    self.structureContainer = Frame(self.window)
    self.resultsContainer = Frame(self.window)
    self.buttonContainer = Frame(self.window)

    self.planContainer.pack(anchor=W, fill=X, expand=True)
    self.structureContainer.pack(anchor=W, fill=X, expand=True)
    self.resultsContainer.pack(anchor=W, fill=X, expand=True)
    self.buttonContainer.pack(fill=X, expand=True)

    Label(self.planContainer, text="Change plan names (* for wildcard)\n FROM \t  TO", font=("Helvetica", 11) ).pack(pady=(15,0))
    #Frame(self.planContainer, bg="grey", relief=SUNKEN).pack(fill=X,expand=1,anchor=W)

    self.planContainers = list()

    for _ in range(1):
        self.planContainers.append(Frame(self.planContainer))
        self.planContainers[-1].pack(expand=True)
        self.planSubstituteList.append({'from' : StringVar(), 'to' : StringVar()})
        self.planEntries.append([Entry(self.planContainers[-1], textvariable=self.planSubstituteList[0]['from'], width=15),
                                 Entry(self.planContainers[-1], textvariable=self.planSubstituteList[0]['to'], width=15)])

        self.planEntries[-1][0].pack(side=LEFT, expand=True)
        self.planEntries[-1][1].pack(side=RIGHT, expand=True)

    Label(self.structureContainer, text="Change structure names (* for wildcard)\n FROM\t  TO", font=("Helvetica", 11) ).pack(pady=(15,0), expand=True)
    #Frame(self.structureContainer, bg="grey", relief=SUNKEN).pack(fill=X,expand=1,anchor=W)

    self.structureContainers = list()

    for _ in range(1):
        self.structureContainers.append(Frame(self.structureContainer))
        self.structureContainers[-1].pack(expand=True)
        self.structureSubstituteList.append({'from' : StringVar(), 'to' : StringVar()})
        self.structureEntries.append([Entry(self.structureContainers[-1], textvariable=self.structureSubstituteList[0]['from'], width=15),
                                      Entry(self.structureContainers[-1], textvariable=self.structureSubstituteList[0]['to'], width=15)])
        self.structureEntries[-1][0].pack(side=LEFT, expand=True)
        self.structureEntries[-1][1].pack(side=RIGHT, expand=True)

    Label(self.resultsContainer, text="Plan and structure names", font=("Helvetica", 12) ).pack(pady=(15,0), expand=True)
    Frame(self.resultsContainer, bg="grey", relief=SUNKEN).pack(fill=X,expand=1,anchor=W)

    self.planAndStructureNameVariable = StringVar("")
    self.planAndStructureNameList = Label(self.resultsContainer, textvariable=self.planAndStructureNameVariable, wraplength=500)
    self.planAndStructureNameList.pack(fill=X, expand=True)

    self.calculateNewNamesCommand()
    self.drawPlanAndStructureNames()

    self.calculateNewNamesButton = Button(self.buttonContainer, text="Calculate new names (Enter)", command=self.calculateNewNamesCommand, width=self.button_width)
    self.calculateNewNamesButton.pack(side=LEFT)
    self.changeNamingQuitAndSaveButton = Button(self.buttonContainer, text="Save and close", command=self.changeNamingQuitAndSaveCommand, width=self.button_width)
    self.changeNamingQuitAndSaveButton.pack(side=LEFT)
    self.changeNamingQuitButton = Button(self.buttonContainer, text="Cancel (Esc)", command=self.changeNamingQuitCommand, width=self.button_width)
    self.changeNamingQuitButton.pack(side=LEFT)
    
    self.window.bind("<Escape>", lambda event=None: self.changeNamingQuitButton.invoke())
    self.window.bind("<Return>", lambda event=None: self.calculateNewNamesButton.invoke())
    # self.window.bind("s", lambda event=None: self.changeNamingQuitAndSaveButton.invoke())


def calculateNewNamesCommand(self):
    def sub(a,b,txt):
        a = a.replace("*", ".*")
        b = b.replace("*", ".*")
        return re.sub(a,b,txt)
    
    self.structuresAfter.clear()
    self.plansAfter.clear()
    numberOfSubstitutions = {'plan' : 0, 'structure' : 0}

    for d in self.planSubstituteList:
        if len(d['from'].get()):
            numberOfSubstitutions['plan'] += 1
    for d in self.structureSubstituteList:
        if len(d['from'].get()):
            numberOfSubstitutions['structure'] += 1

    for plan in self.plansBefore:
        for d in self.planSubstituteList:
            if len(d['from'].get()):
                plan = sub(d['from'].get(), d['to'].get(), plan)
        self.plansAfter.add(plan)

    for structure in self.structuresBefore:
        for d in self.structureSubstituteList:
            if len(d['from'].get()):
                structure = sub(d['from'].get(), d['to'].get(), structure)
        self.structuresAfter.add(structure)

    self.drawPlanAndStructureNames()

    # Add entries if list is full
    if numberOfSubstitutions['plan'] == len(self.planSubstituteList):
        self.planContainers.append(Frame(self.planContainer))
        self.planContainers[-1].pack(expand=True)
        self.planSubstituteList.append({'from' : StringVar(), 'to' : StringVar()})        
        self.planEntries.append([Entry(self.planContainers[-1], textvariable=self.planSubstituteList[-1]['from'], width=15),
                                 Entry(self.planContainers[-1], textvariable=self.planSubstituteList[-1]['to'], width=15)])

        self.planEntries[-1][0].pack(side=LEFT, expand=True)
        self.planEntries[-1][1].pack(side=RIGHT, expand=True)

    if numberOfSubstitutions['structure'] == len(self.structureSubstituteList):
        self.structureContainers.append(Frame(self.structureContainer))
        self.structureContainers[-1].pack(expand=True)
        self.structureSubstituteList.append({'from' : StringVar(), 'to' : StringVar()})
        self.structureEntries.append([Entry(self.structureContainers[-1], textvariable=self.structureSubstituteList[-1]['from'], width=15),
                                      Entry(self.structureContainers[-1], textvariable=self.structureSubstituteList[-1]['to'], width=15)])
        self.structureEntries[-1][0].pack(side=LEFT, expand=True)
        self.structureEntries[-1][1].pack(side=RIGHT, expand=True)


def changeNamingQuitCommand(self):
    self.window.destroy()


def changeNamingQuitAndSaveCommand(self):
    def sub(a,b,txt):
        a = a.replace("*", ".*")
        b = b.replace("*", ".*")
        return re.sub(a,b,txt)
    
    for patientsInCohort in list(self.patients.values()):
        for patient in list(patientsInCohort.patients.values()):
            plan = patient.getPlan()
            if not plan in self.plansAfter:
                for d in self.planSubstituteList:
                    if len(d['from'].get()):
                        plan = sub(d['from'].get(), d['to'].get(), plan)
                patient.setPlan(plan)

            structure = patient.getStructure()
            if not structure in self.structuresAfter:
                for d in self.structureSubstituteList:
                    if len(d['from'].get()):
                        structure = sub(d['from'].get(), d['to'].get(), structure)
                patient.setStructure(structure)
    self.window.destroy()


def autodetectDVHHeaderCommand(self):
    if self.options.autodetectDVHHeader.get() == 0:
        self.customDVHHeaderEntry['state'] = 'normal'
    else:
        self.customDVHHeaderEntry['state'] = 'disabled'


def drawPlanAndStructureNames(self):
    planNamesBefore = ", ".join(sorted(list(self.plansBefore)))
    structureNamesBefore = ", ".join(sorted(list(self.structuresBefore)))
    planNamesAfter = ", ".join(sorted(list(self.plansAfter)))
    structureNamesAfter = ", ".join(sorted(list(self.structuresAfter)))

    self.planAndStructureNameVariable.set(f"Plan names before:\n{planNamesBefore}\n\nPlan names after:\n{planNamesAfter}"
                                    f"\n\nStructure names before:\n{structureNamesBefore}\n\nStructure names after:\n{structureNamesAfter}")


def showGEUDvsN(self):
    for patientsInCohort in list(self.patients.values()):
        for patient in list(patientsInCohort.patients.values()):
            x = np.arange(self.options.nFrom.get(), self.options.nTo.get(), 0.02)
            y = [patient.getGEUD(k) for k in x]
            plt.plot(x,y, "-", color=patient.getTox() >= self.options.toxLimit.get() and "red" or "black")
    plt.xlabel("n values")
    plt.ylabel("gEUD [Gy]")
    plt.title("Pre-calculated gEUD vs n (organ seriality parameter)")
    if not self.options.nIsLinear.get():
        plt.xscale('log')

    custom_lines = dict()
    custom_lines["Toxicity"] = Line2D([0], [0], color="r", ls="-", lw=2)
    custom_lines["No toxicity"] = Line2D([0], [0], color="k", ls="-", lw=2)
    plt.legend(custom_lines.values(), custom_lines.keys())
    
    plt.show()


def switchNto(self):
    if self.options.NTCPcalculation.get() == "LKB":
        if self.options.fixN.get(): 
            self.paramNEntry[-1]['state'] = 'disabled'
        else:
            self.paramNEntry[-1]['state'] = 'normal'
    else:
        if self.options.fixA.get(): 
            self.paramNEntry[-1]['state'] = 'disabled'
        else:
            self.paramNEntry[-1]['state'] = 'normal'

    # IN ANY CASE RESET BESTPARAMETER LIST
    for patientsInCohort in self.patients.values():
        patientsInCohort.bestParameters = list()
        patientsInCohort.resetIdx()

def switchMto(self):
    if self.options.NTCPcalculation.get() == "LKB":
        if self.options.fixM.get(): 
            self.paramMEntry[-1]['state'] = 'disabled'
        else:
            self.paramMEntry[-1]['state'] = 'normal'
    else:
        if self.options.fixB.get(): 
            self.paramMEntry[-1]['state'] = 'disabled'
        else:
            self.paramMEntry[-1]['state'] = 'normal'

    # IN ANY CASE RESET BESTPARAMETER LIST
    for patientsInCohort in self.patients.values():
        patientsInCohort.bestParameters = list()
        patientsInCohort.resetIdx()


def switchTD50to(self):
    if self.options.fixTD50.get(): 
        self.paramTD50Entry[-1]['state'] = 'disabled'
    else:
        self.paramTD50Entry[-1]['state'] = 'normal'
    
    # IN ANY CASE RESET BESTPARAMETER LIST
    for patientsInCohort in self.patients.values():
        patientsInCohort.bestParameters = list()
        patientsInCohort.resetIdx()

def resetIdx(self):
    self.idx = {'n': 0, 'm': 1, 'TD50': 2, 'lambda': 3, 'gamma': 4}
    if self.options.fixN.get():
        self.idx['m'] -= 1
        self.idx['TD50'] -= 1
        self.idx['lambda'] -= 1
        self.idx['gamma'] -= 1

    if self.options.fixM.get():
        self.idx['TD50'] -= 1
        self.idx['lambda'] -= 1
        self.idx['gamma'] -= 1

    if self.options.fixTD50.get():
        self.idx['lambda'] -= 1
        self.idx['gamma'] -= 1

    if self.options.fixLambda.get():
        self.idx['gamma'] -= 1

    if self.options.NTCPcalculation.get() == "Logit":  # Reset lambda / gamma paramters if so
        self.idx = {'a': 0, 'b': 1, 'lambda': 2, 'gamma': 3}
        if self.options.fixA.get():
            self.idx['b'] -= 1
            self.idx['lambda'] -= 1
            self.idx['gamma'] -= 1

        if self.options.fixB.get():
            self.idx['lambda'] -= 1
            self.idx['gamma'] -= 1

        if self.options.fixLambda.get():
            self.idx['gamma'] -= 1