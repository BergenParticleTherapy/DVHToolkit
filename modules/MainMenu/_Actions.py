import os, re, random, time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from functools import partial

from tkinter import *
from tkinter import ttk
from tkinter import filedialog

from ..Patients import *

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
    
    if self.options.DVHFileType.get() == "ECLIPSE":
        plans = self.patients[cohortName].findPlans()
        structures = self.patients[cohortName].findStructures()

        self.options.structureToUse.set("")
        self.options.planToUse.set("")
        if len(structures)>1:
            self.window = Toplevel(self)
            sList = ", ".join(sorted(structures))
            Label(self.window, text=f"The following structures are available, choose one (* for wildcard, | for multiple):\n\n{sList}\n", wraplength=750).pack()
            e1 = Entry(self.window, textvariable=self.options.structureToUse, width=20)
            e1.pack()
            e1.focus()
        else:
            self.options.StructureToUse.set(structures[0])
        
        if len(plans)>1:
            if not self.window:
                self.window = Toplevel(self)
            pList = ", ".join(sorted(plans))
            Label(self.window, text=f"\n\nThe following plans are available, choose one (* for wildcard, | for multiple):\n\n{pList}\n", wraplength=750).pack()
            e2 = Entry(self.window, textvariable=self.options.planToUse, width=20)
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
        log = self.patients[cohortName].loadPatientsSIMPLE(self.progress)
        self.log("\n".join(log))
        self.progress['maximum'] = np.sum([len(k.patients) for k in self.patients.values()])
        self.addPatientCohort(cohortName, self.options.structureToUse.get(), self.options.planToUse.get(), 
                              self.patients[cohortName].getNPatients(), self.patients[cohortName].getNPatientsWithTox())

def addPatientCohort(self, cohortName, structureName, planName, nPatients, nPatientsWithTox):     
    self.cohortList[cohortName] = [Frame(self.middleLeftLowerContainer)]
    self.cohortList[cohortName][0].pack(anchor=N)
    self.cohortList[cohortName].append(Button(self.cohortList[cohortName][0], 
                                        text="DEL", command=partial(self.removePatientCohort, cohortName), width=7, height=3, bg='white'))
    self.cohortList[cohortName].append(Label(self.cohortList[cohortName][0], 
                                        text=f"{planName}\n{structureName}\n{nPatients} patients\n({nPatientsWithTox} with tox)", bg='white', width=15))
    for i in range(1,3):
        self.cohortList[cohortName][i].pack(side=LEFT, pady=i and 0 or 5)
    
    self.buttonShowDVH['state'] = 'normal'
    self.buttonCalculateGEUD['state'] = 'normal'
    
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
    
    elif not hasGEUDs and self.options.NTCPcalculation.get() == "Logit":
        self.buttonCalculateNTCP['state'] = 'normal'
        self.buttonCalculateAUROC['state'] = 'normal'
        self.buttonCalculateDVH['state'] = 'normal'
        self.buttonAggregateDVH['state'] = 'normal'
        self.buttonCalculateGEUD['state'] = 'normal'
    else:
        self.buttonCalculateNTCP['state'] = 'disabled'
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
    if not self.patients:
        self.buttonShowDVH['state'] = 'disabled'
        self.buttonCalculateGEUD['state'] = 'disabled'

def chooseStructureCommand(self,cohortName):
    self.log(f"Choosing structure: {self.options.structureToUse.get()}\n")
    self.window.destroy()
    res = self.patients[cohortName].loadPatientsECLIPSE(self.progress)
    self.addPatientCohort(cohortName, self.options.structureToUse.get(), self.options.planToUse.get(), self.patients[cohortName].getNPatients(), self.patients[cohortName].getNPatientsWithTox())
    self.log("\n".join(res))

def calculateGEUDCommand(self):
    for k,v in list(self.patients.items()):
        self.log(f"Calculating gEUD values for cohort {k}...")
        v.createGEUDs(self.progress)
    self.log("Done.")

    self.buttonCalculateNTCP['state'] = 'normal'
    self.buttonCalculateAUROC['state'] = 'normal'
    self.buttonShowGEUDvsN['state'] = 'normal'
    self.buttonCalculateDVH['state'] = 'normal'
    
def toxLimitChange(self):
    if self.cohortList:
        for cohortName, patientCohort in list(self.cohortList.items()):
            patientCohort[-1].config(text=f"{cohortName}\n{self.patients[cohortName].getNPatients()} patients\n({self.patients[cohortName].getNPatientsWithTox()}) with tox)")
    
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
    if self.options.DVHFileType.get() == "ECLIPSE":
        self.customDVHHeaderEntry['state'] = 'disabled'
    else:
        self.customDVHHeaderEntry['state'] = 'normal'
    
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
    for text,mode in [["Plan", "plan"], ["Structure", "structure"], ["None", "none"]]:
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
    styles = ["-", "--", "-.", ":", "loosely dashed", "dashdotted", "dashdotdotted"]

    plans = set()
    structures = set()

    for patientsInCohort in list(self.patients.values()):
        for patient in list(patientsInCohort.patients.values()):
            plans.add(patient.getPlan())
            structures.add(patient.getStructure())

    plotDict = dict()
    idx = 0
    if self.options.dvhPlotSeparatePlots.get() == "plan":
        for plan in plans:
            plotDict[plan] = idx
            plt.figure(idx, figsize=(15,10))
            plt.title(f"DVHs for {plans}")
            plt.xlabel("Dose [Gy]")
            plt.ylabel("Volume [%]")
            plt.ylim([0, 105])
            idx += 1
    elif self.options.dvhPlotSeparatePlots.get() == "structure":
        for structure in structures:
            plotDict[structure] = idx
            plt.figure(idx, figsize=(15,10))
            plt.title(f"DVHs for {structure}")
            plt.xlabel("Dose [Gy]")
            plt.ylabel("Volume [%]")
            plt.ylim([0, 105])
            idx += 1
    else:
        plotDict[0] = 1
        plt.figure(1, figsize=(15,10))
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
    
    for patientsInCohort in list(self.patients.values()):
        for patient in list(patientsInCohort.patients.values()):
            if self.options.dvhPlotLineStyleGrouping.get() == "plan":
                style = styleDict[patient.getPlan()]
            elif self.options.dvhPlotLineStyleGrouping.get() == "structure":
                style = styleDict[patient.getStructure()]
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
        if self.options.dvhPlotLineStyleGrouping.get() in ["plan", "structure"]:
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
    
def calculateNTCPCommand(self):        
    cohortList = list(self.patients.values())
    for cohort in cohortList:
        cohort.options = self.options
    
    if self.options.NTCPcalculation.get() == "Logit":
        for patients in cohortList:
            patients.calculateDpercent(self.options.NTCPcalculationDpercent.get())
            patients.bestParameters = self.bestParameters
    
    primaryCohort = cohortList[0]
    secondaryCohorts = len(cohortList) > 1 and cohortList[1:] or {}

    if self.options.optimizationScheme.get() == "GradientDescent":
        cStr = secondaryCohorts and "s" or ""
        pList = ",".join(list(self.patients.keys()))
        self.log(f"Performing {self.options.optimizationScheme.get()} ({self.options.optimizationMetric.get()}) optimization on cohort{cStr}: {pList}.")
        
        res = primaryCohort.doGradientOptimization(secondaryCohorts, self.progress)
        
        self.bestParameters = res.x
        self.log("\n".join([f"{k}={v}" for k,v in list(res.items())]))
        
    elif self.options.optimizationScheme.get() == "MatrixMinimization":
        self.log("Calculating toxicity array with matrix size {self.options.matrixSize.get()}")
        res = primaryCohort.doMatrixMinimization(secondaryCohorts, self.progress)
        
        self.bestParameters = res.x
        self.log("\n".join([f"{k}={v}" for k,v in list(res.items())]))
    
    self.buttonLKBuncert['state'] = 'normal'
    primaryCohort.drawSigmoid(secondaryCohorts, self.confidenceInterval, self.correlationLogit, self.aHist, self.bHist, self.TD50Hist, self.LLHhist, self.log)

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
        
        elif not hasGEUDs and self.options.NTCPcalculation.get() == "Logit":
            self.buttonCalculateNTCP['state'] = 'normal'
            self.buttonCalculateAUROC['state'] = 'normal'
            self.buttonCalculateDVH['state'] = 'normal'
            self.buttonCalculateGEUD['state'] = 'normal'
        else:
            self.buttonCalculateNTCP['state'] = 'disabled'
            self.buttonCalculateAUROC['state'] = 'disabled'
            self.buttonCalculateDVH['state'] = 'normal'
            self.buttonCalculateGEUD['state'] = 'normal'

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
    self.radioContainer = Frame(self.inputContainer)
    self.radioContainer.pack(side=LEFT, anchor=W)
    self.entryContainer = Frame(self.inputContainer)
    self.entryContainer.pack(side=LEFT, anchor=W)
    self.calculateMeanDoseContainer = Frame(self.window)
    self.calculateMeanDoseContainer.pack(anchor=W, fill=X, expand=1)
    self.buttonContainer = Frame(self.window)
    self.buttonContainer.pack(anchor=W, fill=X, expand=1)        

    self.dvhCheckVarDoseAtVolume = IntVar(value=1)
    self.dvhCheckVarVolumeAtDose = IntVar(value=1)
     
    self.dvhEntryVar1 = StringVar(value=0)
    self.dvhEntryVar2 = StringVar(value=0)
    self.outputFileNameVar = StringVar(value="Output/dvhValues.csv")
    
    Label(self.fileContainer, text="Output file: ").pack(side=LEFT)
    Entry(self.fileContainer, textvariable=self.outputFileNameVar, width=30).pack(side=LEFT, anchor=W)
    Checkbutton(self.radioContainer, text="Find dose [Gy] at volume [%] (a,b,...) ", variable=self.dvhCheckVarDoseAtVolume).pack(anchor=W)
    Checkbutton(self.radioContainer, text="Find volume [%] at dose [Gy] (a,b,...) ", variable=self.dvhCheckVarVolumeAtDose).pack(anchor=W)
    Entry(self.entryContainer, textvariable=self.dvhEntryVar1).pack(anchor=W)
    Entry(self.entryContainer, textvariable=self.dvhEntryVar2).pack(anchor=W)
    Checkbutton(self.calculateMeanDoseContainer, text="Mean dose from Eclipse? ", variable=self.calculateMeanDose).pack(anchor=W)
    
    b = Button(self.buttonContainer, text="Calculate", command=self.calculateDVHvalues, width=self.button_width)
    b.pack(side=LEFT, anchor=W)
    b2 = Button(self.buttonContainer, text="Cancel", command=self.cancelCalculateDVHvalues, width=self.button_width)
    b2.pack(side=LEFT, anchor=W)

    self.window.bind('<Return>', lambda event=None: b.invoke())
    self.window.bind('<Escape>', lambda event=None: b2.invoke())     
    
def cancelCalculateDVHvalues(self):
    plt.close("all")
    self.window.destroy()

def aggregateDVHCommand(self):
    # Mean / median
    # Subtract or show all?
    
    self.window = Toplevel(self)
    self.window.title("Aggregate DVH values")
    self.window.focus()
    self.styleContainer = Frame(self.window)
    self.styleContainer.pack(anchor=W)
    self.styleContainer2 = Frame(self.window)
    self.styleContainer2.pack(anchor=W)
    self.styleContainer3 = Frame(self.window)
    self.styleContainer3.pack(anchor=W)
    self.buttonContainer = Frame(self.window)
    self.buttonContainer.pack(anchor=W, fill=X, expand=1)
    
    self.dvhStyleVar1 = StringVar(value="mean")
    self.dvhStyleVar2 = StringVar(value="compare")
    self.dvhStyleVar3 = IntVar(value=1)
    
    Label(self.styleContainer, text="Cohort aggregation style: ").pack(side=LEFT)
    for text, mode in [["Median", "median"], ["Mean", "mean"]]:
        Radiobutton(self.styleContainer, text=text, value=mode, variable=self.dvhStyleVar1).pack(anchor=W)
        
    for text, mode in [["Compare cohorts", "compare"], ["Show all", "showAll"], ["Subtract [mean/median of all px]", "subtract"],
                        ["Mean/median [subtract per patient]", "subtractPerPatient"]]:
        Radiobutton(self.styleContainer2, text=text, value=mode, variable=self.dvhStyleVar2).pack(anchor=W)

    Label(self.styleContainer3, text="Draw 83% confidence interval: ").pack(side=LEFT)
    for text, mode in [["Yes", 1], ["No",0]]:
        Radiobutton(self.styleContainer3, text=text, value=mode, variable=self.dvhStyleVar3).pack(side=LEFT, anchor=W)
        
    b = Button(self.buttonContainer, text="Show aggregated DVH", command=self.calculateAggregatedDVH, width=self.button_width)
    b.pack(side=LEFT, anchor=W)
    self.window.bind('<Return>', lambda event=None: b.invoke())
    Button(self.buttonContainer, text="Cancel", command=self.cancelCalculateDVHvalues, width=self.button_width).pack(side=LEFT, anchor=W)

def showGEUDvsN(self):
    for patientsInCohort in list(self.patients.values()):
        for patient in list(patientsInCohort.patients.values()):
            x = np.arange(self.options.nFrom.get(), self.options.nTo.get(), 0.02)
            y = [patient.getGEUD(k) for k in x]
            plt.plot(x,y, "-", color=patient.getTox() >= self.options.toxLimit.get() and "red" or "black")
    plt.xlabel("n values")
    plt.ylabel("gEUD [Gy]")
    plt.title("Pre-calculated gEUD vs n (organ seriality parameter)")

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
    
def switchTD50to(self):
    if self.options.fixTD50.get(): 
        self.paramTD50Entry[-1]['state'] = 'disabled'
    else:
        self.paramTD50Entry[-1]['state'] = 'normal'

