# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scipy.stats as st
from scipy.optimize import basinhopping
import os, re, random, time
from math import *
from functools import partial
from bisect import bisect_left

try: # Python 3
    from tkinter import *
    from tkinter import ttk
    from tkinter import filedialog
except: # Python 2.x
    import tkFileDialog as filedialog
    import Tkinter
    import ttk
    import cProfile, pstats, StringIO

Gy = 1
dGy = 0.1
cGy = 0.01
mGy = 0.001
cc = 1

PROGRAM_VERSION = 1.35

# Changelog version 1.1 (Released 2019-02-15)
# - Added LLH + parameter bound check for ParametricBootstrap
# - Added change in basin hopping parameter names for D% logit model
# - Fixed patient.calculateNTCP (wrong formula!!), used for ParametricBooststrap
# - For plotting of individual BSd models in the sigmoid, now using CI limits on parameters
#        for model selection (draw/nodraw) instead of LLH limit -- the latter cannot be
#        trusted since its limit depend on the cohort in question (and the cohort is BSd)
#   -> More trust in CI method now. Best to use for tox asymmetry is profile likelihood.
# - Fixed ECLIPSE loading
# - Added "plan" options for choosing Eclipse patients
# - Switched to "python" engine for reading CSV files, a bit more stable...


# Changelog version 1.11 (Released 2019-02-15)
# - Added support for multiple structures in DVH set for the same patient
# - Fixed extrapolated DVH file I/O

# Changelog version 1.12 (Released 2019-02-18)
# - Autoselect ECLIPSE input with support for two columns (Dose, Volume)
# - Changed DVH calculation to simple interpolation for increased stability
# - Added Mean Dose readout from ECLIPSE files

# Changelog version 1.2
# - Finalized and verified confidence interval calculations, using bootstrap subtractions
#   - Removed "pivotal CI" option, substituted with "bootstrap subtraction correction" 
#     which acts on the optimal identified parameter set as well

# Changelog version 1.21
# - Added DVH plotting options

# Changelog version 1.3
# - Full support for python 3
# - Added 2-sample t-test CI tooltip
# - The CI calculation now acts directly on TD50 which is the variable of interest
#       (in logit this is -a/b)
# - More DVH plotting options (compare median values of tox / non tox vs cohorts)

# Changelog version 1.31
# - Fixed ECLIPSE input in python 3
# - More DVH plotting options

# Changelog version 1.32
# - Cleaned the output a bit
# - Still needs a bit of work. Shows conflicting TD50 values in Logit!!

# Changelog version 1.33
# - Added lots of DVH plotting options: Mean value of all cohort, compare between struct + plan combos

# Changelog version 1.34
# - Added specific NTCP output (adjust parameters manually in lines 1378 / 1397

# Changelog version 1.35
# - Added more flexibility in the DVH calculations (both ways + comma separated list)
# - Big update in the flexibility of the DVH plotting (grouping, multiple windows, storing plots)

# TODO:
# Clean up unneccessary code for CI calculation (remove non TD50 & non-pivotal data? less output)

def HPM(x):
    """An accurate but quick approximation of the gaussian cumulative density function. 
        See https://www.hindawi.com/journals/mpe/2012/124029/"""
    
    return 1/(exp(-358*x/23 + 111*atan(37*x/294))+1)

class Tooltip:
    '''
    It creates a tooltip for a given widget as the mouse goes on it.

    see:

    http://stackoverflow.com/questions/3221956/what-is-the-simplest-way-to-make-tooltips-
           in-tkinter/36221216#36221216

    http://www.daniweb.com/programming/software-development/
           code/484591/a-tooltip-class-for-tkinter

    - Originally written by vegaseat on 2014.09.09.

    - Modified to include a delay time by Victor Zaccardo on 2016.03.25.

    - Modified
        - to correct extreme right and extreme bottom behavior,
        - to stay inside the screen whenever the tooltip might go out on
          the top but still the screen is higher than the tooltip,
        - to use the more flexible mouse positioning,
        - to add customizable background color, padding, waittime and
          wraplength on creation
      by Alberto Vassena on 2016.11.05.
    '''

    def __init__(self, widget,
                 bg='#FFFFEA',
                 pad=(5, 3, 5, 3),
                 text='widget info',
                 waittime=400,
                 wraplength=250):

        self.waittime = waittime  # in milliseconds, originally 500
        self.wraplength = wraplength  # in pixels, originally 180
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.onEnter)
        self.widget.bind("<Leave>", self.onLeave)
        self.widget.bind("<ButtonPress>", self.onLeave)
        self.bg = bg
        self.pad = pad
        self.id = None
        self.tw = None

    def onEnter(self, event=None):
        self.schedule()

    def onLeave(self, event=None):
        self.unschedule()
        self.hide()

    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(self.waittime, self.show)

    def unschedule(self):
        id_ = self.id
        self.id = None
        if id_:
            self.widget.after_cancel(id_)

    def show(self):
        def tip_pos_calculator(widget, label, 
                    tip_delta=(10, 5), pad=(5, 3, 5, 3)):

            w = widget

            s_width, s_height = w.winfo_screenwidth(), w.winfo_screenheight()

            width, height = (pad[0] + label.winfo_reqwidth() + pad[2],
                             pad[1] + label.winfo_reqheight() + pad[3])

            mouse_x, mouse_y = w.winfo_pointerxy()

            x1, y1 = mouse_x + tip_delta[0], mouse_y + tip_delta[1]
            x2, y2 = x1 + width, y1 + height

            x_delta = x2 - s_width
            if x_delta < 0:
                x_delta = 0
            y_delta = y2 - s_height
            if y_delta < 0:
                y_delta = 0

            offscreen = (x_delta, y_delta) != (0, 0)

            if offscreen:
                if x_delta:
                    x1 = mouse_x - tip_delta[0] - width

                if y_delta:
                    y1 = mouse_y - tip_delta[1] - height

            offscreen_again = y1 < 0  # out on the top
            if offscreen_again: y1 = 0

            return x1, y1

        bg = self.bg
        pad = self.pad
        widget = self.widget

        # creates a toplevel window
        self.tw = Toplevel(widget)

        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)

        win = Frame(self.tw,
                       background=bg,
                       borderwidth=0)
        label = Label(win,
                          text=self.text,
                          justify=LEFT,
                          background=bg,
                          relief=SOLID,
                          borderwidth=0,
                          wraplength=self.wraplength)

        label.grid(padx=(pad[0], pad[2]),
                   pady=(pad[1], pad[3]),
                   sticky=NSEW)
        win.grid()

        x, y = tip_pos_calculator(widget, label)

        self.tw.wm_geometry("+%d+%d" % (x, y))

    def hide(self):
        tw = self.tw
        if tw:
            tw.destroy()
        self.tw = None

class Options():
    def __init__(self):
        self.DVHFileType = StringVar(value = "simple") # ['simple', 'ECLIPSE']
        self.NTCPcalculation = StringVar(value = "Dpercent") # ['LKB', 'Dpercent']
        self.NTCPcalculationDpercent = DoubleVar(value = 20) # [0 -> 100]
        self.dataFolder = StringVar(value = "")
        self.structureToUse = StringVar(value = "")
        self.planToUse = StringVar(value = "")
        self.customDVHHeader = StringVar(value = "Dose,Volume")
        self.CSVStyle = StringVar(value = "autodetect") # ['periodComma', 'commaSemicolon', 'autodetect']
        self.skipRows = IntVar(value=0)
        self.doseUnit = StringVar(value="cGy") # incl. autodetect
        self.optimizationScheme = StringVar(value = "GradientDescent") # ['GradientDescent', 'MatrixMinimization']
        self.optimizationMetric = StringVar(value = "LLH") # ['LLH', 'LS']
        self.confidenceIntervalMethod = StringVar(value = "NonParametricBootstrapping") # ['ProfileLikelihood', 'ParametricBootstrapping', 'NonParametricBootstrapping']
        self.confidenceIntervalIterations = IntVar(value = 1500)
        self.confidenceIntervalPercent = DoubleVar(value = 95)
        self.confidenceIntervalLikelihoodLimit = DoubleVar(value = -1.0)
        self.bootstrapCorrectionMethod = StringVar(value = "median") # ["none", "mean", "median"]
        self.confidenceIntervalShowModels = IntVar(value = 0)
        self.toxLimit = IntVar(value = 2)
        self.matrixSize = IntVar(value = 50)
        self.loadToxFromFilename = IntVar(value = 1)
        
        self.fixN = IntVar(value = 0)
        self.nFrom = DoubleVar(value = 0.02)
        self.nTo = DoubleVar(value = 1)
        
        self.fixM = IntVar(value = 0)
        self.mFrom = DoubleVar(value = 0.03)
        self.mTo = DoubleVar(value = 1)
        
        self.fixTD50 = IntVar(value = 0)
        self.TD50From = DoubleVar(value = 10)
        self.TD50To = DoubleVar(value = 85)
        
        self.fixA = IntVar(value = 0)
        self.aFrom = DoubleVar(value = -100)
        self.aTo = DoubleVar(value = 0)
        
        self.fixB = IntVar(value = 0)
        self.bFrom = DoubleVar(value = 0)
        self.bTo = DoubleVar(value = 2)
        
        self.basinHoppingIterations = IntVar(value = 10)
        self.basinHoppingTemperature = DoubleVar(value = 0.1)
        self.basinHoppingNsize = DoubleVar(value = 0.3)
        self.basinHoppingMsize = DoubleVar(value = 0.3)
        self.basinHoppingTD50size = DoubleVar(value = 10)
        self.basinHoppingAsize = DoubleVar(value = 40)
        self.basinHoppingBsize = DoubleVar(value = 0.5)
        
        self.dvhPlotUseToxAsColor = IntVar(value = 0)
        self.dvhPlotLegendMarker = StringVar(value = "structureName") # ["structureName", "folderName", "planName"]
        self.dvhPlotLegendSize = DoubleVar(value = 14)
        self.dvhPlotLineStyleGrouping = StringVar(value="plan")
        self.dvhPlotSeparatePlots = StringVar(value="structure")
        self.dvhPlotsSaveFigs = IntVar(value = 1)
        
        self.vars = {"DVHFileTypeS" : self.DVHFileType, "NTCPcalculationS" : self.NTCPcalculation, 
               "NTCPcalculationDpercentI" : self.NTCPcalculationDpercent, "dataFolderS" : self.dataFolder, 
               "structureToUseS" : self.structureToUse, "CSVStyleS" : self.CSVStyle,
               "optimizationSchemeS" : self.optimizationScheme, "optimizationMetricS" : self.optimizationMetric,
               "toxLimitI" : self.toxLimit, "loadToxFromFilenameI" : self.loadToxFromFilename, "matrixSizeI" : self.matrixSize,
               "fixNI" : self.fixN, "nFromD" : self.nFrom, "nToD" : self.nTo, 
               "fixMI" : self.fixM, "mFromD" : self.mFrom, "mToD" : self.mTo,
               "fixTD50I" : self.fixTD50, "TD50FromD" : self.TD50From, "Td50ToD" : self.TD50To,
               "fixA" : self.fixA, "aFrom" : self.aFrom, "aTo" : self.aTo,
               "fixB" : self.fixB, "bFrom" : self.bFrom, "bTo" : self.bTo,
               "basinHoppingIterationsI" : self.basinHoppingIterations, "basinHoppingTemperatureD" : self.basinHoppingTemperature,
               "basinHoppingNsizeD" : self.basinHoppingNsize, "basinHoppingMsizeD" : self.basinHoppingMsize, 
               "basinHoppingTD50sizeD" : self.basinHoppingTD50size, "confidenceIntervalMethod": self.confidenceIntervalMethod,
               "confidenceIntervalIterations" : self.confidenceIntervalIterations, "confidenceIntervalPercent" : self.confidenceIntervalPercent,
               "customDVHHeader" : self.customDVHHeader, "skipRows" : self.skipRows, "doseUnit" : self.doseUnit,
               "confidenceIntervalShowModels" : self.confidenceIntervalShowModels,
               "bootstrapCorrectionMethod" : self.bootstrapCorrectionMethod, "basinHoppingAsizeD": self.basinHoppingAsize,
               "basinHoppingBsizeD" : self.basinHoppingBsize, "dvhPlotUseToxAsColor" : self.dvhPlotUseToxAsColor,
               "dvhPlotLegendMarker" : self.dvhPlotLegendMarker, "dvhPlotLineStyleGrouping" : self.dvhPlotLineStyleGrouping,
               "dvhPlotSeparatePlots" : self.dvhPlotSeparatePlots, "dvhPlotsSaveFigs" : self.dvhPlotsSaveFigs,
               "dvhPlotLegendSize" : self.dvhPlotLegendSize}

    def loadOptions(self):
        read = False
        if os.path.exists('config.cfg'):
            with open("config.cfg", "r") as configFile:
                for line in configFile.readlines():
                    linesplit = line.rstrip().split(",")
                    var = linesplit[0]
                    value = linesplit[1]
                    if value:
                        read = True
                        if var == "customDVHHeader":
                            value = ",".join(linesplit[1:])
                            
                        if var in list(self.vars.keys()): 
                            self.vars[var].set(value)
        return read

    def saveOptions(self):
        with open("config.cfg","w") as configFile:
            for key, var in list(self.vars.items()):
                configFile.write("{},{}\n".format(key, var.get()))
                
    def getCustomHeaderIdx(self):
        try:
            header = self.customDVHHeader.get().split(",")
        except:
            header = ["Dose","Volume"]
            
        doseIdx = header.index("Dose")
        volumeIdx = header.index("Volume")
        return "{}{}".format(doseIdx, volumeIdx)

class MainMenu(Frame):
    def __init__(self, parent):
        Frame.__init__(self, parent)
        self.parent = parent
        self.patients = {}
        self.cohortList = {}
        self.parent.protocol("WM_DELETE_WINDOW", self.myQuit)
        self.parent.title("DVH tool with NTCP calculator {:.2f} - Helge Pettersen".format(PROGRAM_VERSION))
        self.window = None
        
        self.wraplength = 250
        
        self.button_width=25
        self.bestParameters = []
        self.bestParametersNone = []
        self.bestParametersMean = []
        self.bestParametersMedian = []
        self.confidenceInterval = [[0,0],[0,0], [0,0]]
        self.confidenceIntervalNone = [[0,0],[0,0], [0,0]]
        self.confidenceIntervalMean = [[0,0],[0,0], [0,0]]
        self.confidenceIntervalMedian = [[0,0],[0,0], [0,0]]
        self.correlationLogit = -0.016
        self.calculateMeanDose = IntVar(value = 1)
        self.aHist = []
        self.bHist = []
        self.TD50Hist = []
        self.LLHhist = []
        
        self.paramNEntry = []
        self.paramMEntry = []
        self.paramTD50Entry = []
        self.paramBHsizeLabel = []
        self.paramBHsizeEntry = []
        
        self.options = Options()
        res = self.options.loadOptions()
        
        self.upperContainer = Frame(self, bd=5, relief=RIDGE, height=40)  # Title
        self.middleContainer = Frame(self, bd=5)
        self.bottomContainer = Frame(self, bd=20)

        self.middleLeftContainer = Frame(self.middleContainer, bd=5) # Patient cohorts
        self.middleLeftUpperContainer = Frame(self.middleLeftContainer, bd=5) # Load patient cohort
        self.middleLeftLine = Frame(self.middleLeftContainer, bg="grey", relief=SUNKEN)
        self.middleLeftLowerContainer = Frame(self.middleLeftContainer, bd=5) # List patient cohorts                            
        self.middleMiddleLine = Frame(self.middleContainer, bg="grey", relief=SUNKEN)
        self.middleMiddleContainer = Frame(self.middleContainer, bd=5) # Options
        self.middleRightLine = Frame(self.middleContainer, bg="grey", relief=SUNKEN)
        self.middleRightContainer = Frame(self.middleContainer, bd=5) # Log
        self.middleRightContainerLine = Frame(self.middleRightContainer, bg="grey", relief=SUNKEN)
        self.middleRightLowerContainer = Frame(self.middleRightContainer)
        self.bottomLine = Frame(self.bottomContainer, bg="grey", relief=SUNKEN)
        self.bottomContainer1 = Frame(self.bottomContainer)
        self.bottomContainer2 = Frame(self.bottomContainer)
        
        self.dvhFileContainer = Frame(self.middleMiddleContainer)
        self.toxLimitContainer = Frame(self.middleMiddleContainer)
        self.toxFromFilenameContainer = Frame(self.middleMiddleContainer)
        self.NTCPcalculationContainer = Frame(self.middleMiddleContainer)
        self.NTCPpercentContainer = Frame(self.middleMiddleContainer)
        self.CSVStyleContainer = Frame(self.middleMiddleContainer)
        self.customDVHHeaderContainer = Frame(self.middleMiddleContainer)
        self.doseUnitContainer = Frame(self.middleMiddleContainer)
        self.skipRowsContainer = Frame(self.middleMiddleContainer)
        self.optimizationSchemeContainer = Frame(self.middleMiddleContainer)
        self.optimizationMetricContainer = Frame(self.middleMiddleContainer)
        self.matrixSizeContainer = Frame(self.middleMiddleContainer)
        self.basinHoppingsIterationsContainer = Frame(self.middleMiddleContainer)
        self.basinHoppingsTemperatureContainer = Frame(self.middleMiddleContainer)
        self.basinHoppingsNsizeContainer = Frame(self.middleMiddleContainer)
        self.basinHoppingsMsizeContainer = Frame(self.middleMiddleContainer)
        self.basinHoppingsTD50sizeContainer = Frame(self.middleMiddleContainer)
        self.basinHoppingsJumpLenghtsContainer = Frame(self.middleMiddleContainer)
        self.confidenceIntervalPercentContainer = Frame(self.middleMiddleContainer)
        self.confidenceIntervalSchemeContainer = Frame(self.middleMiddleContainer)
        self.confidenceIntervalIterationsContainer = Frame(self.middleMiddleContainer)
        self.bootstrapCorrectionMethodContainer = Frame(self.middleMiddleContainer)
        self.confidenceIntervalLikelihoodLimitContainer = Frame(self.middleMiddleContainer)
        self.confidenceIntervalShowModelsContainer = Frame(self.middleMiddleContainer)
        
        self.upperContainer.pack(fill=X)
        self.middleContainer.pack(fill=Y)
        self.middleLeftContainer.pack(side=LEFT,fill="both", expand=1, anchor=N)
        self.middleLeftUpperContainer.pack(fill=X)
        self.middleLeftLine.pack(fill=X, padx=5, pady=5)
        self.middleLeftLowerContainer.pack(fill=X)
        self.middleMiddleLine.pack(side=LEFT, fill=Y, padx=5, pady=5, expand=1)
        self.middleMiddleContainer.pack(side=LEFT, fill=Y, padx=5, pady=5)
        self.middleRightLine.pack(side=LEFT, fill=Y, padx=5, pady=5, expand=1)
        self.middleRightContainer.pack(side=LEFT,fill=Y)
        self.middleRightContainerLine.pack(fill=X, expand=1)
        self.middleRightLowerContainer.pack(fill=X, expand=1)
        self.bottomLine.pack(fill=X, padx=5, pady=5, expand=1)
        self.bottomContainer.pack(fill=X, anchor=N, expand=1)
        self.bottomContainer1.pack(anchor=N, expand=1)
        self.bottomContainer2.pack(anchor=N, expand=1)
        
        Label(self.upperContainer, text="DVH tool with NTCP calculator {:.2f} - Helge Pettersen".format(PROGRAM_VERSION)).pack(anchor=N)
        
        self.loadPatientsButton = Button(self.middleLeftUpperContainer, text="Load patient folder (F)", command=self.loadPatientsCommand, width=self.button_width)
        self.loadPatientsButton.pack(anchor=N, pady=3)
        Tooltip(self.loadPatientsButton, text="Load patient cohort. The file should contain DVH files in a CSV format. "
                "Upon loading, additional gEUD files will searched for in a gEUD/ subfolder. If they are not found, they should "
                "be created using the \"Calculate gEUD splines\" action prior to any LKB calculations.", wraplength=self.wraplength)

        self.parent.bind("f", lambda event=None: self.loadPatientsButton.invoke())
        
        self.progress = ttk.Progressbar(self.middleLeftUpperContainer, orient=HORIZONTAL, maximum=100, mode='determinate')
        self.progress.pack(fill=X, pady=3)
                
        Label(self.middleMiddleContainer, text="OPTIONS", font=("Helvetica", 12)).pack(anchor=N)
        Frame(self.middleMiddleContainer, bg="grey", relief=SUNKEN).pack(fill=X,expand=1,anchor=W)

        
        self.dvhFileContainer.pack(anchor=W)
        Label(self.dvhFileContainer, text="DVH File type: ").pack(side=LEFT, anchor=W)
        for mode in ["simple", "ECLIPSE"]:
            Radiobutton(self.dvhFileContainer, text=mode, variable=self.options.DVHFileType, value=mode, command=self.selectDVHFileTypeCommand).pack(side=LEFT, anchor=W)
        Tooltip(self.dvhFileContainer, text="Outline of the per-patient CSV files. The \"simple\" style assumes a flat CSV file. "
                "The \"ECLIPSE\" style includes several metadata lines in the beginning, incl. descriptions of the DVH organ / structure.", wraplength=self.wraplength)
        
        self.CSVStyleContainer.pack(anchor=W)
        Label(self.CSVStyleContainer, text="Decimal / Column separator: ").pack(side=LEFT, anchor=W)
        for text, mode in [[", / ;", "commaSemicolon"], [". / ,", "periodComma"], ["Autodetect", "autodetect"]]:
            Radiobutton(self.CSVStyleContainer, text=text, variable=self.options.CSVStyle, value=mode).pack(side=LEFT, anchor=W)
        Tooltip(self.CSVStyleContainer, text="CSV style: Comma and Semicolon, or Period and Comma for decimal and column separation.", wraplength=self.wraplength)
                
        self.customDVHHeaderContainer.pack(anchor=W)
        Label(self.customDVHHeaderContainer, text="Columns in DVH files: ").pack(side=LEFT, anchor=W)
        self.customDVHHeaderEntry = Entry(self.customDVHHeaderContainer, textvariable=self.options.customDVHHeader, width=50)
        self.customDVHHeaderEntry.pack(side=LEFT, anchor=W)
        Tooltip(self.customDVHHeaderContainer, 
                text="Define headers as a comma-separated list, indicating the wanted Dose and Volume columns without whitespace. Example: "
                "\"Volume,dummy,dummy,Dose\" or \"Dose,Volume\".", wraplength=self.wraplength)
    
        self.doseUnitContainer.pack(anchor=W)
        Label(self.doseUnitContainer, text="Dose unit: ").pack(side=LEFT, anchor=W)
        for text, mode in [ ["cGy", "cGy"], ["Gy", "Gy"], ["Autodetect", "autodetect"]]:
            Radiobutton(self.doseUnitContainer, text=text, variable=self.options.doseUnit, value=mode).pack(side=LEFT, anchor=W)
    
        self.skipRowsContainer.pack(anchor=W)
        Label(self.skipRowsContainer, text="Number of rows to skip in simple csv: ").pack(side=LEFT, anchor=W)
        Entry(self.skipRowsContainer, textvariable=self.options.skipRows, width=5).pack(side=LEFT)        
        
        self.toxFromFilenameContainer.pack(anchor=W)
        Label(self.toxFromFilenameContainer, text="Load toxicity from \'tox\' in filename: ").pack(side=LEFT, anchor=W)
        for text, mode in [('yes', 1), ('no', 0)]:
            Radiobutton(self.toxFromFilenameContainer, text=text, variable=self.options.loadToxFromFilename, 
                        value=mode, command=self.toxFromFilenameCommand).pack(side=LEFT, anchor=W)        
        Tooltip(self.toxFromFilenameContainer, text="Where to get the toxicity information from. Yes: use the presence of \'tox\' in the filename "
                "to define a yes/no toxicity. No: provide a CSV file \"Data/{folderName}_tox.csv\" with lines patientID,toxGrade. The patientID "
                "is the filename.", wraplength=self.wraplength)
        
        self.toxLimitContainer.pack(anchor=W)
        Label(self.toxLimitContainer, text="Toxicity limit (grade): ").pack(side=LEFT, anchor=W)
        for mode in range(5):
            Radiobutton(self.toxLimitContainer, text=mode, variable=self.options.toxLimit, value=mode, command=self.toxLimitChange).pack(side=LEFT, anchor=W)
        Tooltip(self.toxLimitContainer, 
                text="The toxicity threshold. Patients with the chosen grade or higher will be classified with complication.", wraplength=self.wraplength)
        
        self.NTCPcalculationContainer.pack(anchor=W)
        Label(self.NTCPcalculationContainer, text="NTCP calculation: ").pack(side=LEFT, anchor=W)
        for text, mode in [("D% + logit", "Dpercent"), ("LKB", "LKB")]:
            Radiobutton(self.NTCPcalculationContainer, text=text, variable=self.options.NTCPcalculation, 
                        value=mode, command=self.NTCPcalculationCommand).pack(side=LEFT, anchor=W)
        Tooltip(self.NTCPcalculationContainer, text="The parameter limits for the LKB model are set below. "
                "For the D% + logit model, only the percentage value has to be set.", wraplength=self.wraplength)
        
        self.NTCPpercentContainer.pack(anchor=W)
        Label(self.NTCPpercentContainer, text="NTCP calculation Dose value: ").pack(side=LEFT, anchor=W)
        self.NTCPpercentLabel = Entry(self.NTCPpercentContainer, textvariable=self.options.NTCPcalculationDpercent, width=5)
        self.NTCPpercentLabel.pack(side=LEFT, anchor=W)
        Label(self.NTCPpercentContainer, text="%").pack(side=LEFT, anchor=W)
        Tooltip(self.NTCPcalculationContainer, text="The D% value from the DVHs to be used in the logit model.", wraplength=self.wraplength)
                
        Label(self.middleMiddleContainer, text="ALGORITHM OPTIONS", font=("Helvetica", 12) ).pack(pady=(15,0))
        Frame(self.middleMiddleContainer, bg="grey", relief=SUNKEN).pack(fill=X,expand=1,anchor=W)

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
        self.paramTD50Entry.append(Checkbutton(self.basinHoppingsTD50sizeContainer, text="Fixed", variable=self.options.fixTD50, command=self.switchTD50to))
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
        Tooltip(self.basinHoppingsTemperatureContainer, text="The expected difference in error between the different local minima.", wraplength=self.wraplength)
        
        self.basinHoppingsJumpLenghtsContainer.pack(anchor=W)
        Label(self.basinHoppingsJumpLenghtsContainer, text="G.D. \"basin hopping\" jump lengths: ").pack(side=LEFT, anchor=W)
        for text, mode in [("n: ", self.options.basinHoppingNsize), ("m: ",self.options.basinHoppingMsize), ("TD50: ", self.options.basinHoppingTD50size)]:
            self.paramBHsizeLabel.append(Label(self.basinHoppingsJumpLenghtsContainer, text=text))
            self.paramBHsizeLabel[-1].pack(side=LEFT, anchor=W)
            self.paramBHsizeEntry.append(Entry(self.basinHoppingsJumpLenghtsContainer, textvariable=mode, width=5))
            self.paramBHsizeEntry[-1].pack(side=LEFT, anchor=W)

        Tooltip(self.basinHoppingsJumpLenghtsContainer, 
                text="For each \"basin hopping\", each parameter is perturbed as a random number within ± these values.", wraplength=self.wraplength)

        self.confidenceIntervalPercentContainer.pack(anchor=W)
        Label(self.confidenceIntervalPercentContainer, text="Confidence Interval percentage: ").pack(side=LEFT, anchor=W)
        Entry(self.confidenceIntervalPercentContainer, textvariable=self.options.confidenceIntervalPercent, width=7).pack(side=LEFT, anchor=W)
        Label(self.confidenceIntervalPercentContainer, text="%").pack(side=LEFT, anchor=W)
        Tooltip(self.confidenceIntervalPercentContainer, text="1 sigma: 68.75, 2 sigma: 95.45, 3 sigma: 99.73. 2-way p<0.05 @ = 83.4", wraplength=self.wraplength)
        
        self.confidenceIntervalLikelihoodLimitContainer.pack(anchor=W)
        Label(self.confidenceIntervalLikelihoodLimitContainer, text="Log Likelihood limit: ").pack(side=LEFT, anchor=W)
        Entry(self.confidenceIntervalLikelihoodLimitContainer, textvariable=self.options.confidenceIntervalLikelihoodLimit, width=7).pack(side=LEFT, anchor=W)
        Tooltip(self.confidenceIntervalLikelihoodLimitContainer, text="This is to filter out \"bad\" or at-limit parameter sets. Check LLH histogram to find this value, "
                    "typically -1 to 0 to remove parameters where the LLH is >2 too high and close to 0.", wraplength=self.wraplength)
        
        self.confidenceIntervalSchemeContainer.pack(anchor=W)
        Label(self.confidenceIntervalSchemeContainer, text="CI technique: ").pack(side=LEFT)
        for text, mode in [("Non-parametric BS", "NonParametricBootstrapping"), ("Parametric BS", "ParametricBootstrapping"), ("Profile likelihood","ProfileLikelihood")]:
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
                "Usually 1000-2000 for good statistics, this step might be a bit time demanding.", wraplength=self.wraplength)
        
        if self.options.DVHFileType.get() == "ECLIPSE":
            self.customDVHHeaderEntry['state'] = 'disabled'
        
        if self.options.NTCPcalculation.get() == "Dpercent":
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
        
        # Label(self.middleRightContainer, text="LOG", width=15).pack()
        self.logtext = Text(self.middleRightLowerContainer, height=38, width=75)
        self.scrollbar = Scrollbar(self.middleRightLowerContainer)
        self.scrollbar.config(command=self.logtext.yview)
        self.logtext.config(yscrollcommand=self.scrollbar.set)
        self.scrollbar.pack(side=RIGHT, fill="y", expand=False)            
        self.logtext.pack(side=LEFT, fill="both", expand=True)
        self.log("----- LOG -----")
        
        self.buttonCalculateGEUD = Button(self.bottomContainer1, text="Calculate gEUD splines (G)", command=self.calculateGEUDCommand, width=self.button_width, state=DISABLED)
        self.buttonShowGEUDvsN = Button(self.bottomContainer1, text="Show gEUD", command=self.showGEUDvsN, width=self.button_width, state=DISABLED)
        self.buttonShowDVH = Button(self.bottomContainer1, text="Show DVHs (D)", command=self.showDVHCommand, width=self.button_width, state=DISABLED)
        self.buttonCalculateDVH = Button(self.bottomContainer1, text="Calculate DVH values (C)", command=self.calculateDVHCommand, width=self.button_width, state=DISABLED)
        self.buttonAggregateDVH = Button(self.bottomContainer1, text="Aggregate DVH values (A)", command=self.aggregateDVHCommand, width=self.button_width, state=DISABLED)
        self.buttonCalculateNTCP = Button(self.bottomContainer2, text="Calculate NTCP model (N)", command=self.calculateNTCPCommand, width=self.button_width, state=DISABLED)
        self.buttonLKBuncert = Button(self.bottomContainer2, text="Calculate confidence intervals (B)", command=self.calculateLKBuncert, width=self.button_width, state=DISABLED)        
        self.buttonCalculateAUROC = Button(self.bottomContainer2, text="Calculate AUROC", command=self.calculateAUROCCommand, width=self.button_width, state=DISABLED)
        self.buttonQuit = Button(self.bottomContainer2, text="Exit (Esc)", command=self.myQuit, width=self.button_width)

        self.parent.bind("g", lambda event=None: self.buttonCalculateGEUD.invoke())
        self.parent.bind("d", lambda event=None: self.buttonShowDVH.invoke())
        self.parent.bind("c", lambda event=None: self.buttonCalculateDVH.invoke())
        self.parent.bind("a", lambda event=None: self.buttonAggregateDVH.invoke())
        self.parent.bind("n", lambda event=None: self.buttonCalculateNTCP.invoke())
        self.parent.bind("b", lambda event=None: self.buttonLKBuncert.invoke())
        self.parent.bind("<Escape>", lambda event=None: self.buttonQuit.invoke())
        
        for button in [self.buttonCalculateGEUD, self.buttonShowGEUDvsN, self.buttonShowDVH, self.buttonCalculateDVH, self.buttonAggregateDVH,
                  self.buttonCalculateNTCP, self.buttonLKBuncert, self.buttonCalculateAUROC, self.buttonQuit]:
            button.pack(side=LEFT, anchor=N, padx=5, pady=5)

        self.pack()
        
        if not res: self.log("Could not load options.cfg, using default values.")

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
                if cohortName == patientCohort[:-3] and patientCohort[-3] == "_" and patientCohort[-2:] in ["{:02d}".format(k) for k in range(10)]:
                    maxCohort = max(maxCohort, int(patientCohort[-2:]))
            
            if maxCohort < 0:
                cohortName += "_01"
            else:
                cohortName += "_{:02d}".format(maxCohort + 1)

            self.log("New cohort already in list of patient cohorts, renaming to %s." % (cohortName))
        
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
                Label(self.window, text="The following structures are available, choose one (* for wildcard):\n\n%s\n" % (", ".join(sorted(structures))), wraplength=750).pack()
                e1 = Entry(self.window, textvariable=self.options.structureToUse, width=20)
                e1.pack()
                e1.focus()
            else:
                self.options.StructureToUse.set(structures[0])
            
            if len(plans)>1:
                if not self.window:
                    self.window = Toplevel(self)
                    
                Label(self.window, text="\n\nThe following plans are available, choose one (* for wildcard):\n\n%s\n" % (", ".join(sorted(plans))), wraplength=750).pack()
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
                                            text=("%s\n%s\n%d patients\n(%d with tox)" % (planName, structureName, nPatients, nPatientsWithTox)), bg='white', width=15))
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
        
        elif not hasGEUDs and self.options.NTCPcalculation.get() == "Dpercent":
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
        self.log("Deleted cohort %s." % cohortName)
        if not self.patients:
            self.buttonShowDVH['state'] = 'disabled'
            self.buttonCalculateGEUD['state'] = 'disabled'
    
    def chooseStructureCommand(self,cohortName):
        self.log("Choosing structure: %s\n" % (self.options.structureToUse.get()))
        self.window.destroy()
        res = self.patients[cohortName].loadPatientsECLIPSE(self.progress)
        self.addPatientCohort(cohortName, self.options.structureToUse.get(), self.options.planToUse.get(), self.patients[cohortName].getNPatients(), self.patients[cohortName].getNPatientsWithTox())
        self.log("\n".join(res))
    
    def calculateGEUDCommand(self):
        for k,v in list(self.patients.items()):
            self.log("Calculating gEUD values for cohort %s..." % k)
            v.createGEUDs(self.progress)
        self.log("Done.")

        self.buttonCalculateNTCP['state'] = 'normal'
        self.buttonCalculateAUROC['state'] = 'normal'
        self.buttonShowGEUDvsN['state'] = 'normal'
        self.buttonCalculateDVH['state'] = 'normal'
        
    def toxLimitChange(self):
        if self.cohortList:
            for cohortName, patientCohort in list(self.cohortList.items()):
                patientCohort[-1].config(text=("%s\n%d patients\n(%d with tox)" % (cohortName, 
                                self.patients[cohortName].getNPatients(), self.patients[cohortName].getNPatientsWithTox())))
        
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
        for text, mode in [["Red (tox) / Black (noTox)", True], ["Rainbow", False]]:
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
                self.log("Best parameters: {} -> {}".format([float("{:.2f}".format(k)) for k in self.bestParameters], 
                                                             [float("{:.2f}".format(k)) for k in self.bestParametersNone]))
                self.log("Confidence Interval: {} -> {}".format(["{:.2f} - {:.2f}".format(k1,k2) for k1,k2 in self.confidenceInterval],
                                                                ["{:.2f} - {:.2f}".format(k1,k2) for k1,k2 in self.confidenceIntervalNone]))
                self.bestParameters = self.bestParametersNone
                self.confidenceInterval = self.confidenceIntervalNone
                
            elif self.options.bootstrapCorrectionMethod.get() == "mean":
                self.log("Changing bootstrap correction method to Mean.")
                self.log("Best parameters: {} -> {}".format([float("{:.2f}".format(k)) for k in self.bestParameters], 
                                                             [float("{:.2f}".format(k)) for k in self.bestParametersMean]))
                self.log("Confidence Interval: {} -> {}".format(["{:.2f} - {:.2f}".format(k1,k2) for k1,k2 in self.confidenceInterval],
                                                                ["{:.2f} - {:.2f}".format(k1,k2) for k1,k2 in self.confidenceIntervalMean]))
                self.bestParameters = self.bestParametersMean
                self.confidenceInterval = self.confidenceIntervalMean
                
            elif self.options.bootstrapCorrectionMethod.get() == "median":
                self.log("Changing bootstrap correction method to Median.")
                self.log("Best parameters: {} -> {}".format([float("{:.2f}".format(k)) for k in self.bestParameters], 
                                                             [float("{:.2f}".format(k)) for k in self.bestParametersMedian]))
                self.log("Confidence Interval: {} -> {}".format(["{:.2f} - {:.2f}".format(k1,k2) for k1,k2 in self.confidenceInterval],
                                                                ["{:.2f} - {:.2f}".format(k1,k2) for k1,k2 in self.confidenceIntervalMedian]))
                self.bestParameters = self.bestParametersMedian
                self.confidenceInterval = self.confidenceIntervalMedian
        
    def calculateNTCPCommand(self):        
        cohortList = list(self.patients.values())
        for cohort in cohortList:
            cohort.options = self.options
        
        if self.options.NTCPcalculation.get() == "Dpercent":
            for patients in cohortList:
                patients.calculateDpercent(self.options.NTCPcalculationDpercent.get())
                patients.bestParameters = self.bestParameters
        
        primaryCohort = cohortList[0]
        secondaryCohorts = len(cohortList) > 1 and cohortList[1:] or {}

        # if len(self.bestParameters) == 0: {
        if self.options.optimizationScheme.get() == "GradientDescent":
            self.log("Performing %s (%s) optimization on cohort%s: %s." % (self.options.optimizationScheme.get(), 
                    self.options.optimizationMetric.get(), secondaryCohorts and "s" or "", ",".join(list(self.patients.keys()))))
            
            if self.options.NTCPcalculation.get() == "LKB":
                res = primaryCohort.doOptimization(secondaryCohorts, self.progress)
            else:
                res = primaryCohort.doOptimizationLogistic(secondaryCohorts, self.progress)
                self.log("TD50 (-a/b) = {:.2f} Gy.".format(-res.x[0]/res.x[1]))
                
            self.bestParameters = res.x
            self.log("\n".join(["%s: %s"%(k,v) for k,v in list(res.items())]))
            
        else:
            if self.options.NTCPcalculation.get() == "LKB":
                self.log("Calculating toxicity array for {}^3 (m n TD50) parameters".format(self.options.matrixSize.get()))
                res = primaryCohort.doMatrixMinimization(secondaryCohorts, self.progress)
            else:
                self.log("Calculating toxicity array for {}^2 (a b) parameters".format(self.options.matrixSize.get()))
                res = primaryCohort.doMatrixMinimizationLogistic(secondaryCohorts, self.progress)
            self.bestParameters = res
            self.log(res)
        # }
            
        self.buttonLKBuncert['state'] = 'normal'
        primaryCohort.drawSigmoid(secondaryCohorts, self.confidenceInterval, self.correlationLogit, self.aHist, self.bHist, self.TD50Hist, self.LLHhist, self.log)
        
    def NTCPcalculationCommand(self):
        if self.options.NTCPcalculation.get() == "Dpercent":
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
            
            elif not hasGEUDs and self.options.NTCPcalculation.get() == "Dpercent":
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

        
    def calculateDVHvalues(self):
        nValues = {'Bladder' : 1/8, 'Rectum' : 1/12, 'Intestine' : 1/4}
        x1 = []; x2 = []; x3 = []; x4 = []; y1 = []; y2 = []; y3 = []; y4 = []; y5 = []; y6 = []

        if self.dvhCheckVarDoseAtVolume.get() and self.dvhCheckVarVolumeAtDose.get():
            split = self.dvhEntryVar1.get().split(",")
            volumeFromList = [float(k) for k in split]
            volumesFromList = ", ".join([f"{k:.3f}" for k in volumeFromList])
            doseFromList = [float(k) for k in self.dvhEntryVar2.get().split(",")]
            dosesFromList = ", ".join([f"{k:.3f}" for k in doseFromList])
         
            self.log(f"Calculating dose in all patients where the volume fractions are {volumesFromList} + doses are {dosesFromList}.")
            for patientsInCohort in list(self.patients.values()):
                for name, patient in list(patientsInCohort.patients.items()):
                    try:
                        meanDose = patient.getMeanDoseFromEclipse()
                    except:
                        meanDose = 0
                    doseList = [patient.getDoseAtVolume(k) for k in volumeFromList]
                    doses = ", ".join([f"{k:.3f}" for k in doseList])
                    volumeList = [patient.getVolumeAtDose(k) for k in doseFromList]
                    volumes = ", ".join([f"{k:.3f}" for k in volumeList])

                    structure = patient.getStructure()
                    plan = patient.getPlan()
#                    NTCP = HPM((patient.getGEUD(self.options.nFrom) - self.options.TD50From)/(self.options.mFrom * self.options.TD50From))
                    geud = 0
                    if structure in nValues:
                        geud = patient.getGEUD(nValues[structure])
                    x1.append(patientsInCohort.cohort)
                    x2.append(name.split("_")[0])
                    x3.append(structure)
                    x4.append(plan)
                    y1.append(doses)
                    y2.append(volumes)
                    y3.append(meanDose)
                    y4.append(geud)
                    y5.append(patient.getMaxDoseFromEclipse())
                    y6.append(patient.getMinDoseFromEclipse())

        elif self.dvhCheckVarDoseAtVolume.get():
            split = self.dvhEntryVar1.get().split(",")
            volumeFromList = [float(k) for k in split]
            volumes = ", ".join([f"{k:.3f}" for k in volumeFromList])
            volumes = ", ".join(["-" for k in volumeFromList])
         
            self.log(f"Calculating dose in all patients where the volume fraction is {volumes}.")
            for patientsInCohort in list(self.patients.values()):
                for name, patient in list(patientsInCohort.patients.items()):
                    try:
                        meanDose = patient.getMeanDoseFromEclipse()
                    except:
                        meanDose = 0
                    doseList = [patient.getDoseAtVolume(k) for k in volumeFromList]
                    doses = ", ".join([f"{k:.3f}" for k in doseList])
                    structure = patient.getStructure()
                    plan = patient.getPlan()
#                    NTCP = HPM((patient.getGEUD(self.options.nFrom) - self.options.TD50From)/(self.options.mFrom * self.options.TD50From))
                    geud = 0
                    if structure in nValues:
                        geud = patient.getGEUD(nValues[structure])
                    x1.append(patientsInCohort.cohort)
                    x2.append(name.split("_")[0])
                    x3.append(structure)
                    x4.append(plan)
                    y1.append(doses)
                    y2.append(volumes)
                    y3.append(meanDose)
                    y4.append(geud)
                    y5.append(patient.getMaxDoseFromEclipse())
                    y6.append(patient.getMinDoseFromEclipse())
        
        elif self.dvhCheckVarVolumeAtDose.get():
            doseFromList = [float(k) for k in self.dvhEntryVar2.get().split(",")]
            doses = ", ".join([f"{k:.3f}" for k in doseFromList])
            doses = ", ".join(["-" for k in doseFromList])
            
            self.log(f"Calculating volume in all patients where the dose is {doses} Gy.")
            for patientsInCohort in list(self.patients.values()):
                for name, patient in list(patientsInCohort.patients.items()):
                    volumeList = [patient.getVolumeAtDose(k) for k in doseFromList]
                    volumes = ", ".join([f"{k:.3f}" for k in volumeList])
                    structure = patient.getStructure()
                    plan = patient.getPlan()
                    try:
                        meanDose = patient.getMeanDoseFromEclipse()
                    except:
                        meanDose = 0
#                    NTCP = HPM((patient.getGEUD(self.options.nFrom) - self.options.TD50From)/(self.options.mFrom * self.options.TD50From))
                    geud = 0
                    if structure in nValues:
                        geud = patient.getGEUD(nValues[structure])
                    x1.append(patientsInCohort.cohort)
                    x2.append(name.split("_")[0])
                    x3.append(structure)
                    x4.append(plan)
                    y1.append(doses)
                    y2.append(volumes)
                    y3.append(meanDose)
                    y4.append(geud)
                    y5.append(patient.getMaxDoseFromEclipse())
                    y6.append(patient.getMinDoseFromEclipse())
        
        with open(self.outputFileNameVar.get(), "w") as dvhFile:
            dosesFromList = ", ".join([f"D{k}%" for k in self.dvhEntryVar1.get().split(",")])
            volumesFromList = ", ".join([f"V{k}Gy" for k in self.dvhEntryVar2.get().split(",")])
            meanDoseText = self.calculateMeanDose and ",Mean Dose (ECLIPSE) [Gy],Min Dose (ECLIPSE),Max Dose (ECLIPSE)" or ""
            dvhFile.write(f"Cohort,Structure,Plan,Patient ID,{dosesFromList},{volumesFromList},gEUD [Gy]{meanDoseText}\n")
            for cohort, name, structure, plan, dose,volume, meanDose, geud, maxDose, minDose, in zip(x1,x2, x3, x4, y1, y2, y3, y4, y5, y6):
                meanDoseValue = self.calculateMeanDose and f",{meanDose}"
                minDoseValue = self.calculateMeanDose and f",{minDose}"
                maxDoseValue = self.calculateMeanDose and f",{maxDose}"
                dvhFile.write(f"{cohort},{structure},{plan},{name},{dose},{volume},{geud}{meanDoseValue}{minDoseValue}{maxDoseValue}\n")

        self.window.destroy()
        
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

    
    def calculateAggregatedDVH(self):
        self.window.destroy()
        cohortDVH = {}
        tox = {}
        notox = {}
        for cohort_, patientsInCohort in list(self.patients.items()):
            plan = list(patientsInCohort.patients.values())[0].getPlan()
            structure = list(patientsInCohort.patients.values())[0].getStructure()
            cohort = f"{structure}/{plan}"
            cohortDVH[cohort] = None
            tox[cohort] = []
            notox[cohort] = []
            
            first = True
            if self.dvhStyleVar2.get() == "showAll":
                for name, patient in list(patientsInCohort.patients.items()):
                    if patient.getTox() >= self.options.toxLimit.get():
                        tox[cohort].append("Volume_{}".format(name))
                    else:
                        notox[cohort].append("Volume_{}".format(name))
                    
                    namePx = name.split("_")[0]
                    if "vmat" in namePx:
                        namePx = namePx[:-4]
                    else:
                        namePx = namePx[:-1]
                        
                    if first:
                        cohortDVH[cohort] = pd.DataFrame({"Dose": patient.dvh["Dose"], "Volume_{}".format(name) : patient.dvh["Volume"]})
                        cohortDVH[cohort].set_index("Dose", inplace=True)
                        first = False
                    else:
                        newDVH = pd.DataFrame({"Dose": patient.dvh["Dose"], "Volume_{}".format(name) : patient.dvh["Volume"]})
                        newDVH.set_index("Dose", inplace=True)
                        cohortDVH[cohort] = cohortDVH[cohort].merge(newDVH, how="outer", right_index=True, left_index=True)

                cohortDVH[cohort] = cohortDVH[cohort].interpolate(method='index', limit_direction='backward', limit = 100).fillna(0)

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
            else: # COMPARE / SUBTRACT
                plan_structures = []
                for name, patient in list(patientsInCohort.patients.items()):
                    plan = patient.getPlan()
                    structure = patient.getStructure()
                    if plan == "Clinical(1)":
                        plan = "Clinical"
                    cohort = f"{structure}/{plan}"
                    namePx = name.split("_")[0]
                        
                    if not f"{plan}_{structure}" in plan_structures:
                        plan_structures.append(f"{plan}_{structure}")
                        cohortDVH[cohort] = pd.DataFrame({"Dose": patient.dvh["Dose"], "Volume_{}".format(name) : patient.dvh["Volume"]})
                        cohortDVH[cohort].set_index("Dose", inplace=True)
                    else:
                        newDVH = pd.DataFrame({"Dose": patient.dvh["Dose"], "Volume_{}".format(name) : patient.dvh["Volume"]})
                        newDVH.set_index("Dose", inplace=True)
                        cohortDVH[cohort] = cohortDVH[cohort].merge(newDVH, how="outer", right_index=True, left_index=True)

                cohortDVH[cohort] = cohortDVH[cohort].interpolate(method='index', limit_direction='backward', limit = 100).fillna(0)

                for cohort in cohortDVH.keys():
                    if self.dvhStyleVar1.get() == "mean":
                        cohortDVH[cohort]["Volume agg"] = cohortDVH[cohort].mean(axis=1) # tox and notox
                    else:
                        cohortDVH[cohort]["Volume agg"] = cohortDVH[cohort].median(axis=1)

                    # 83.4% CI (two-sided p-test < 0.05): [ 8.3, 91.7 ]
                    cohortDVH[cohort]["Volume agg 8.3"] = cohortDVH[cohort].quantile(0.083, axis=1)
                    cohortDVH[cohort]["Volume agg 91.7"] = cohortDVH[cohort].quantile(0.917, axis=1)
                    cohortDVH[cohort]["Volume agg 5"] = cohortDVH[cohort].quantile(0.05, axis=1)
                    cohortDVH[cohort]["Volume agg 95"] = cohortDVH[cohort].quantile(0.95, axis=1)
                    
        if self.dvhStyleVar2.get() == "showAll": # Show all aggregated cohorts
            colors = ['limegreen', 'violet', 'gold', 'orangered', 'crimson', 'lightcoral', 'firebrick']
            style = ['-', '--']
            idx=0
            for k,v in cohortDVH.items():
                v["Volume agg tox"].plot(use_index=True, linestyle=style[idx%2], color=colors[idx//2], label="{} tox".format(k))
                idx += 1
                v["Volume agg notox"].plot(use_index=True, linestyle=style[idx%2], color=colors[idx//2], label="{} no tox".format(k))
                idx += 1
            plt.xlabel("Dose [Gy]")
            plt.ylabel("Volume [%]")
            plt.xlim([0, 75])
            plt.legend()
            plt.show()

        elif self.dvhStyleVar2.get() == "compare": # Compare plans
            # colors = ['limegreen', 'violet', 'gold', 'orangered', 'crimson', 'lightcoral', 'firebrick']
            
            plans = set([k.split("/")[1] for k in cohortDVH.keys()])
            structures = set([k.split("/")[0] for k in cohortDVH.keys()])
            
            styleIdx = {k:idx for idx,k in enumerate(plans)}
            colorIdx = {k:idx for idx,k in enumerate(structures)}

            style = ['-', '--']
            colorSet = {'PTV72.5' : 'darkred', 'PTV67.5' : 'indianred', 'PTV50' : 'red', 'PTV60' : 'salmon',
                        'Rectum' : 'seagreen', 'Intestine' : 'magenta', 'Bowel' : 'magenta', 'Bladder' : 'goldenrod' }
            
            custom_lines = {k:Line2D([0], [0], color="k", ls=k, lw=2) for k in style}
            custom_lines2 = {k:Line2D([0], [0], color=k, ls='-', lw=2) for k in colorSet.values()}

            custom_line_bladder = Line2D([0],[0], color=colorSet['Bladder'], ls='-', lw=2)
            custom_line_rectum = Line2D([0],[0], color=colorSet['Rectum'], ls='-', lw=2)
            custom_line_intestine = Line2D([0],[0], color=colorSet['Intestine'], ls='-', lw=2)

            fig = plt.figure(figsize=(6*1.5,8*1.5))
            plt1 = fig.add_subplot(2,1,1)
            plt2 = fig.add_subplot(2,2,3)
            plt3 = fig.add_subplot(2,2,4)
            plotsStructure = list()
            plotsPlan = list()
            
            for k,v in cohortDVH.items():
                plan = k.split("/")[1]
                structure = k.split("/")[0]

                ls = style[styleIdx[plan]]
                #c = colors[colorIdx[k.split("/")[0]]]
                c = colorSet[structure]
                thisplt=plt1
                if self.dvhStyleVar3.get():
                    """
                    if structure == "Bladder":
                        thisplt = plt1
                    else:
                        thisplt = plt2
                        
                    thisplt.fill_between(v.index, v["Volume agg 5"], v["Volume agg 95"],color=structure=="Bladder" and "gold" or c, alpha=structure=="Bladder" and 0.5 or 0.3)
                    thisplt.plot(v["Volume agg 5"].index, v["Volume agg 5"], linestyle=ls, color=c, linewidth=1, alpha=structure=="Bladder" and 0.9 or 0.7)
                    thisplt.plot(v["Volume agg 95"].index, v["Volume agg 95"], linestyle=ls, color=c, linewidth=1, alpha=structure=="Bladder" and 0.9 or 0.7)
                    """
                    
                p = thisplt.plot(v["Volume agg"].index, v["Volume agg"], linestyle=ls, color=c)
                
                if not [structure, c] in plotsStructure:
                    plotsStructure.append([structure, c])
                if not [plan, ls] in plotsPlan:
                    plotsPlan.append([plan,ls])


            # ONLY FOR CURRENT FIGURE
            for k,v in cohortDVH.items():
                plan = k.split("/")[1]
                structure = k.split("/")[0]

                if not structure in ['Bladder', 'Rectum']: continue

                ls = style[styleIdx[plan]]
                #c = colors[colorIdx[k.split("/")[0]]]
                c = colorSet[structure]
                if self.dvhStyleVar3.get():
                    if structure == "Bladder":
                        thisplt = plt2
                    else:
                        thisplt = plt3
                        
                    thisplt.fill_between(v.index, v["Volume agg 5"], v["Volume agg 95"],color=structure=="Bladder" and "gold" or c, alpha=structure=="Bladder" and 0.5 or 0.3)
                    thisplt.plot(v["Volume agg 5"].index, v["Volume agg 5"], linestyle=ls, color=c, linewidth=1, alpha=structure=="Bladder" and 0.9 or 0.7)
                    thisplt.plot(v["Volume agg 95"].index, v["Volume agg 95"], linestyle=ls, color=c, linewidth=1, alpha=structure=="Bladder" and 0.9 or 0.7)
                    
                p = thisplt.plot(v["Volume agg"].index, v["Volume agg"], linestyle=ls, color=c)
                
                if not [structure, c] in plotsStructure:
                    plotsStructure.append([structure, c])
                if not [plan, ls] in plotsPlan:
                    plotsPlan.append([plan,ls])

            plt.subplots_adjust(.08,.12,.95,.97, .19, .09)
            plt1.set_xlabel("Dose [Gy]", fontsize=12)
            plt2.set_xlabel("Dose [Gy]", fontsize=12)
            plt3.set_xlabel("Dose [Gy]", fontsize=12)
            plt1.set_ylabel("Volume [%]", fontsize=12)
            plt2.set_ylabel("Volume [%]", fontsize=12)
#            plt3.set_ylabel("Volume [%]", fontsize=12)
            plt1.set_xlim([0, 75])
            plt2.set_xlim([0, 75])
            plt3.set_xlim([0, 75])
            all_legends0 = [custom_lines[k[1]] for k in plotsPlan] + [Line2D([],[],linestyle='')] + [custom_lines2[k[1]] for k in plotsStructure]
            all_legends1 = [custom_lines[k[1]] for k in plotsPlan] + [Line2D([],[],linestyle='')] + [custom_line_bladder]
            all_legends2 = [custom_lines[k[1]] for k in plotsPlan] + [Line2D([],[],linestyle='')] + [custom_line_rectum]
            all_labels1 = [k[0] == "Clinical" and "Clinical Plan" or k[0] for k in plotsPlan] + [''] + ['Bladder']
            all_labels2 = [k[0] == "Clinical" and "Clinical Plan" or k[0] for k in plotsPlan] + [''] + ['Rectum']
            all_labels0 = [k[0] == "Clinical" and "Clinical Plan" or k[0] for k in plotsPlan] + [''] + [k[0] == "Intestine" and "Bowel" or k[0] for k in plotsStructure]
            plt1.legend(all_legends0, all_labels0, loc='lower left')
            plt2.legend(all_legends1, all_labels1, loc='lower left')
            plt3.legend(all_legends2, all_labels2, loc='lower left')

            def get_axis_limits(ax, scale=.87):
                return ax.get_xlim()[1]*scale, ax.get_ylim()[1]*scale
            
            plt1.annotate('A', xy=get_axis_limits(plt1,0.95), fontsize=16, weight='bold')
            plt2.annotate('B', xy=get_axis_limits(plt2), fontsize=16, weight='bold')
            plt3.annotate('C', xy=get_axis_limits(plt3), fontsize=16, weight='bold')

            plt.rcParams['font.size'] = '12'            
            plt1.tick_params(labelsize=12)
            plt2.tick_params(labelsize=12)
            plt3.tick_params(labelsize=12)
            txt = "Figure 1: Population mean DVH for a) all OARs and PTVs, b) bladder with 90% CI and c) rectum with 90% CI. " \
                      "Solid lines are used for RP, dashed lines for CP and shaded areas with thin border lines for CI."
            plt.figtext(0.5,0.01,txt,wrap=True,horizontalalignment='center',fontsize=12)
            
            plt.show()
                
        else: # Subtract two cohorts
            style = ['-', '--']
            fig = plt.figure(figsize=(6*1.5,8*1.5))
            colorSet = {'PTV72.5' : 'darkred', 'PTV67.5' : 'indianred', 'PTV50' : 'red', 'PTV60' : 'salmon',
            'Rectum' : 'seagreen', 'Intestine' : 'magenta', 'Bowel' : 'magenta', 'Bladder' : 'goldenrod' }
            
            custom_lines = {k:Line2D([0], [0], color="k", ls=k, lw=2) for k in style}
            custom_lines2 = {k:Line2D([0], [0], color=k, ls='-', lw=2) for k in colorSet.values()}

            plans = sorted(set([k.split("/")[1] for k in cohortDVH.keys()]))
            structures = sorted(set([k.split("/")[0] for k in cohortDVH.keys()]))

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
                cohortDiffPerPatient[structure].drop(["Volume agg", "Volume agg 8.3", "Volume agg 91.7"], axis=1, inplace=True)
                if self.dvhStyleVar1.get() == "mean":
                    cohortDiffPerPatient[structure]["Volume agg"] = cohortDiffPerPatient[structure].mean(axis=1)
                else:
                    cohortDiffPerPatient[structure]["Volume agg"] = cohortDiffPerPatient[structure].median(axis=1)

                if self.dvhStyleVar3.get():
                    cohortDiffPerPatient[structure]["Volume agg 5"] = cohortDiffPerPatient[structure].quantile(0.05, axis=1)
                    cohortDiffPerPatient[structure]["Volume agg 95"] = cohortDiffPerPatient[structure].quantile(0.95, axis=1)

            plotsStructure = list()
            for structure in structures:
                c = colorSet[structure]
                ls = '-'
                
                if self.dvhStyleVar2.get() == "subtract":
                    cohortDiff[structure].plot(use_index=True, color=c, label=structure)
                    plt.title(f"{self.dvhStyleVar1.get().capitalize()} DVH {kLarge} - self.dvhStyleVar1.get().capitalize() DVH {kSmall}")
                elif self.dvhStyleVar2.get() == "subtractPerPatient":
                    cohortDiffPerPatient[structure]["Volume agg"].plot(use_index=True, color=c, label=structure)
                    if self.dvhStyleVar3.get():
                        plt.fill_between(cohortDiffPerPatient[structure].index,
                                         cohortDiffPerPatient[structure]["Volume agg 5"], cohortDiffPerPatient[structure]["Volume agg 95"],
                                         color=structure=="Bladder" and "gold" or c, alpha=structure=="Bladder" and 0.5 or 0.3)
                        
                        plt.plot(cohortDiffPerPatient[structure]["Volume agg 5"].index,
                                 cohortDiffPerPatient[structure]["Volume agg 5"],
                                 linestyle=ls, color=c, linewidth=1, alpha=structure=="Bladder" and 0.9 or 0.7)
                        
                        plt.plot(cohortDiffPerPatient[structure]["Volume agg 95"].index,
                                 cohortDiffPerPatient[structure]["Volume agg 95"], linestyle=ls,
                                 color=c, linewidth=1, alpha=structure=="Bladder" and 0.9 or 0.7)
                        
                        #plt.fill_between(cohortDiffPerPatient[structure].index, cohortDiffPerPatient[structure]["Volume agg 5"],
                                         #cohortDiffPerPatient[structure]["Volume agg 95"], color=c, linewidth=3, alpha=0.25)
                        
                
            txt = "Figure 2: Population mean and 90% CI of the absolute difference in " \
                  "relative volume as a function of dose between RP and CP for each patient."
            plt.figtext(0.5,0.01,txt,wrap=True,horizontalalignment='center',fontsize=12)


            plt.rcParams['font.size'] = '12'
            plt.tick_params(labelsize=12)
            plt.xlabel("Dose [Gy]", fontsize=12)
            plt.ylabel("RapidPlan Volume [%] - Clinical Plan Volume [%]", fontsize=12)
            plt.legend()

            """
            plt.subplot(132)
            dropColumns = ["Volume agg", "Volume agg 2.5", "Volume agg 97.5"]
            melted = [cohortDiffPerPatient[s].drop(dropColumns, axis=1).melt()['value'] for s in structures]
            aggCohortDiff = pd.DataFrame(melted, index=structures).T
            aggCohortDiff.boxplot(whis=[5,95],showfliers=False)
            plt.xlabel("Structures")
            plt.ylabel(f"DVH difference {kLarge} - {kSmall} [%]")

            plt.subplot(143)
            # k = 4,4,8 for tarm,rektum,blære
            k_values = {'Bladder' : 8, 'Intestine' : 4, 'Rectum' : 12}
            geud_value_kLarge = {s:list() for s in structures}
            geud_value_kSmall = {s:list() for s in structures}
            V50_value_kLarge = {s:list() for s in structures}
            V50_value_kSmall = {s:list() for s in structures}
            for cohort_, patientsInCohort in list(self.patients.items()):
                for name, patient in list(patientsInCohort.patients.items()):
                    plan = patient.getPlan()
                    structure = patient.getStructure()
                    #geud = patient.getGEUD(1/k_values[structure])
                    if kLarge == plan:
                        #geud_value_kLarge[structure].append( geud )
                        print("large", name)
                        #V50_value_kLarge[structure].append( patient.getVolumeAtDose(35) )
                        V50_value_kLarge[structure].append( patient.getDoseAtVolume(98) )
                    else:
                        #geud_value_kSmall[structure].append( geud )
                        print("small", name)
                        #V50_value_kSmall[structure].append( patient.getVolumeAtDose(35) )
                        V50_value_kSmall[structure].append( patient.getDoseAtVolume(98) )
            #geud_value_kLarge_pd = pd.DataFrame(geud_value_kLarge)
            #geud_value_kSmall_pd = pd.DataFrame(geud_value_kSmall)
            #geud_difference = geud_value_kLarge_pd - geud_value_kSmall_pd
            #pd.DataFrame(geud_difference).boxplot()
            #plt.ylabel("GEUD")
            
            plt.subplot(133)
            V50_value_kLarge_pd = pd.DataFrame(V50_value_kLarge)
            V50_value_kSmall_pd = pd.DataFrame(V50_value_kSmall)
            V50_difference = V50_value_kLarge_pd - V50_value_kSmall_pd
            V50_difference.boxplot()
            plt.ylabel("V50%")

            """
                    
            plt.show()
                
    def showGEUDvsN(self):
        for patientsInCohort in list(self.patients.values()):
            for patient in list(patientsInCohort.patients.values()):
                x = np.arange(self.options.nFrom.get(), self.options.nTo.get(), 0.02)
                y = [patient.getGEUD(k) for k in x]
                plt.plot(x,y, "-", color=patient.getTox() >= self.options.toxLimit.get() and "red" or "black")
        plt.xlabel("n values")
        plt.ylabel("gEUD [Gy]")
        plt.show()
        
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
        
        nIterations = self.options.confidenceIntervalIterations.get()        
        
        origIt = self.options.basinHoppingIterations.get()
        self.options.basinHoppingIterations.set(2)
        self.progress['maximum'] = nIterations 
        
        nHist = np.zeros(nIterations)
        mHist = np.zeros(nIterations)
        TD50Hist = np.zeros(nIterations)
        LLHhist = np.zeros(nIterations)

        nInit = self.bestParameters[0]
        mInit = self.bestParameters[1]
        if self.options.NTCPcalculation.get() == "LKB":
            TD50Init = self.bestParameters[2]
        else:
            TD50Init = -nInit/mInit

        grade = self.options.toxLimit.get()
        
        for patientsInCohort in list(self.patients.values()):
            patientsInCohort.calculateNTCP()
        
        cohortList = list(self.patients.values())
        primaryCohort = cohortList[0]
        secondaryCohorts = len(cohortList) > 1 and cohortList[1:] or {}
        time1 = time.time()
        if self.options.confidenceIntervalMethod.get() == "ParametricBootstrapping":
            idx = 0
            nTot = 0
            nCorrect = 0
            
            for patientsInCohort in list(self.patients.values()):
                patientsInCohort.saveTox()
            
            for k in range(nIterations):
                self.progress.step(1)
                self.progress.update_idletasks()
                                
                for patientsInCohort in list(self.patients.values()):
                    for patient in list(patientsInCohort.patients.values()):
                        nTot += 1
                        rn = random.random()
                        ntcp = patient.getNTCP()
                        patient.setTox(rn < ntcp and grade or 0)
                        
                        if patient.tox == patient.storedTox:
                            nCorrect += 1
        
                if self.options.NTCPcalculation.get() == "LKB":
                    if self.options.optimizationScheme.get() == "GradientDescent":
                        res = primaryCohort.doOptimization(secondaryCohorts, None)
                    else:
                        res = primaryCohort.doMatrixMinimization(secondaryCohorts, None)
                    
                    if self.options.fixN.get() and res.x[0] == self.options.aTo.get(): # Many models incorrectly? put n=nTo, don't use these for CI calculation
                        continue
                        
                elif self.options.NTCPcalculation.get() == "Dpercent":
                    if self.options.optimizationScheme.get() == "GradientDescent":
                        res = primaryCohort.doOptimizationLogistic(secondaryCohorts, None)
                    else:
                        res = primaryCohort.doMatrixMinimizationLogistic(secondaryCohorts, None)
                    
                    if res.x[0] < self.options.aFrom.get() + 0.5 or res.x[0] > self.options.aTo.get() - 0.5:
                        continue
                    
                    if res.x[1] < self.options.bFrom.get() + 0.01 or res.x[1] > self.options.bTo.get() - 0.01:
                        continue
                
                if res.fun < -self.options.confidenceIntervalLikelihoodLimit.get():
                    continue
                
                nHist[idx] = res.x[0]
                mHist[idx] = res.x[1]
                if self.options.NTCPcalculation.get() == "LKB":
                    TD50Hist[idx] = res.x[2]
                LLHhist[idx] = -res.fun
                idx += 1
            
            for patientsInCohort in list(self.patients.values()):
                patientsInCohort.restoreTox()
            
            print("Of {} patients, {} had correctly guessed tox.".format(nTot, nCorrect))
            
        elif self.options.confidenceIntervalMethod.get() == "NonParametricBootstrapping":
            idx = 0
            patientZip = []
            for cohort, patientsInCohort in list(self.patients.items()):
                for patientName in list(patientsInCohort.patients.keys()):
                    patientZip.append((cohort, patientName))
            nPatients = len(patientZip)
                    
            for k in range(nIterations):
                self.progress.step(1)
                self.progress.update_idletasks()
                
                newPatientCohort = Patients(self.options)
                newPatientCohort.bestParameters = [x for x in self.bestParameters]
                
                nTox = 0
                for n in range(nPatients):
                    thisPatient = Patient(None)
                    randomPatientID = random.randint(0, nPatients-1)
                    
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
                
                if self.options.NTCPcalculation.get() == "LKB":
                    if self.options.optimizationScheme.get() == "GradientDescent":
                        res = newPatientCohort.doOptimization({}, None)
                    elif self.options.optimizationScheme.get() == "MatrixMinimization":
                        res = newPatientCohort.doMatrixMinimization({}, None)
                        
                    if self.options.fixN.get() and res.x[0] == self.options.aTo.get(): # Many models incorrectly? put n=nTo, don't use these for CI calculation
                        continue
                        
                elif self.options.NTCPcalculation.get() == "Dpercent":
                    if self.options.optimizationScheme.get() == "GradientDescent":
                        res = newPatientCohort.doOptimizationLogistic({}, None)
                    elif self.options.optimizationScheme.get() == "MatrixMinimization":
                        res = newPatientCohort.doMatrixMinimizationLogistic({}, None)
                    
                    if res.x[0] < self.options.aFrom.get() + 0.5 or res.x[0] > self.options.aTo.get() - 0.5:
                        continue
                    
                    if res.x[1] < self.options.bFrom.get() + 0.01 or res.x[1] > self.options.bTo.get() - 0.01:
                        continue
                
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
            if self.options.NTCPcalculation.get() == "Dpercent":
                self.log("a = {:.3f} ({:.3f} - {:.3f})".format(self.bestParameters[0], res[0][0], res[0][1]))
                self.log("b = {:.3f} ({:.3f} - {:.3f})".format(self.bestParameters[1], res[1][0], res[1][1]))
            else:
                self.log("n = {:.2f} ({:.2f}-{:.2f})".format(self.bestParameters[0], res[0][0], res[0][1]))
                self.log("m = {:.2f} ({:.2f}-{:.2f})".format(self.bestParameters[1], res[1][0], res[1][1]))
                self.log("TD50 = {:.2f} ({:.2f}-{:.2f})".format(self.bestParameters[2], res[2][0], res[2][1]))
            
            self.confidenceInterval = res
            return
        time2 = time.time()
        self.options.basinHoppingIterations.set(origIt)      

        nHist = np.trim_zeros(nHist)
        mHist = np.trim_zeros(mHist)
        TD50Hist = np.trim_zeros(TD50Hist)
        LLHhist = np.trim_zeros(LLHhist)
        
        if self.options.NTCPcalculation.get() == "LKB":
            with open("Output/LKBuncert{}.csv".format(self.options.confidenceIntervalMethod.get()), "w") as out:
                out.write("n,m,TD50\n")
                for k in range(len(nHist)):
                    out.write("{},{},{}\n".format(nHist[k], mHist[k], TD50Hist[k]))
            
            lowerPercent = (100 - self.options.confidenceIntervalPercent.get()) / 2
            upperPercent = 100 - lowerPercent        
            
            nMedian = np.median(nHist)
            mMedian = np.median(mHist)
            TD50Median = np.median(TD50Hist)
            
            nMean = np.mean(nHist)
            mMean = np.mean(mHist)
            TD50Mean = np.mean(TD50Hist)
            
            print("nMedian = {:.2f}, mMedian = {:.2f}, TD50Median = {:.2f}".format(nMedian, mMedian, TD50Median))
            print("nMean = {:.2f}, mMean = {:.2f}, TD50Mean = {:.2f}".format(nMean, mMean, TD50Mean))
            
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
                
            self.log("\nFinished Confidence Interval tests ({:.1f} minutes).".format((time2-time1)/60))
            self.log("{}% CI calculated as the percentiles of a {} procedure.".format(self.options.confidenceIntervalPercent.get(), 
                             self.options.confidenceIntervalMethod.get()))
            self.log("The original parameters were n = {:.2f}, m = {:.2f} and TD50 = {:.1f} Gy.".format(nInit, mInit, TD50Init))
            if self.options.bootstrapCorrectionMethod.get() == 'mean':
                self.log("Using the MEAN bias correction method, the corrected best fits are:")
                self.log("n = {:.2f};\tm = {:.2f};TD50 = {:.2f}\t".format(2*nInit-nMean, 2*mInit-mMean, 2*TD50Init-TD50Mean))
            elif self.options.bootstrapCorrectionMethod.get() == 'median':
                self.log("Using the MEDIAN bias correction method, the corrected best fits are:")
                self.log("n = {:.2f};\tm = {:.2f};TD50 = {:.2f}\t".format(2*nInit-nMedian, 2*mInit-mMedian, 2*TD50Init-TD50Median))
            if self.options.bootstrapCorrectionMethod.get() == 'none':
                self.log("Bootstrapped confidence intervals:")
                self.log("n\t= [{:.2f},  {:.2f}]".format(nCI[0], nCI[1]))
                self.log("m\t= [{:.2f},  {:.2f}]".format(mCI[0], mCI[1]))
                self.log("TD50\t= [{:.1f},  {:.1f}]".format(TD50CI[0], TD50CI[1]))
            else:
                self.log("\"Pivot\" bootstrapped confidence intervals:")
                self.log("n\t= [{:.2f},  {:.2f}]".format(nCIPivot[0], nCIPivot[1]))
                self.log("m\t= [{:.2f},  {:.2f}]".format(mCIPivot[0], mCIPivot[1]))
                self.log("TD50\t= [{:.1f},  {:.1f}]".format(TD50CIPivot[0], TD50CIPivot[1]))
            
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
                
            fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(1, 6, sharex=False, sharey=False, figsize=(30,5))
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
            
        elif self.options.NTCPcalculation.get() == "Dpercent":
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

            with open("Output/LKBuncert{}.csv".format(self.options.confidenceIntervalMethod.get()), "w") as out:
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
            
            # In case of overshooting the distributions!
            #nCIPivot = [max(nCIPivot[0], nCI[0]), min(nCIPivot[1], nCI[1])]
            #mCIPivot = [max(mCIPivot[0], mCI[0]), min(mCIPivot[1], mCI[1])]
            
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
                print("PUSH UP")
                TD50CIPivot = [k + TD50Lim[0] - TD50CIPivot[0] for k in TD50CIPivot]
            if TD50CIPivot[1] > TD50Lim[1]:
                print("PUSH DOWN")
                TD50CIPivot = [k - (TD50CIPivot[1] - TD50Lim[1]) for k in TD50CIPivot]
            
            self.log("\nFinished Confidence Interval tests ({:.1f} minutes).".format((time2-time1)/60))
            self.log("{}% CI calculated as the percentiles of a {} procedure.".format(self.options.confidenceIntervalPercent.get(), 
                             self.options.confidenceIntervalMethod.get()))
            self.log("The original parameters were a = {:.2f} and b = {:.2f} (-> TD50 = {:.2f} Gy).".format(nInit, mInit, TD50Init))
            if self.options.bootstrapCorrectionMethod.get() == 'mean':
                self.log("Using the MEAN bias correction method, the corrected best fits are:")
                self.log("a = {:.2f};\tb = {:.2f}".format(2*nInit-nMean, 2*mInit-mMean, 2*TD50Init-TD50Mean))
            elif self.options.bootstrapCorrectionMethod.get() == 'median':
                self.log("Using the MEDIAN bias correction method, the corrected best fits are:")
                self.log("a = {:.2f};\tb = {:.2f}\tTD50 = {:.2f}".format(2*nInit-nMedian, 2*mInit-mMedian, 2*TD50Init-TD50Median))
            if self.options.bootstrapCorrectionMethod.get() == 'none':
                self.log("The values of the bootstrapped distributions are shown below as median (CI):")
                self.log("a\t= [{:.2f},  {:.2f}]".format(nCI[0], nCI[1]))
                self.log("b\t= [{:.2f},  {:.2f}]".format(mCI[0], mCI[1]))
                self.log("TD50\t= {:.2f} - {:.2f}.".format(TD50CI[0], TD50CI[1]))
            else:
                self.log("If the \"pivot\" method is applied, the following CIs are obtained:")
                self.log("a\t= {:.2f} - {:.2f}".format(nCIPivot[0], nCIPivot[1]))
                self.log("b\t= {:.2f} - {:.2f}".format(mCIPivot[0], mCIPivot[1]))
                self.log("TD50\t: {:.2f} - {:.2f}.".format(TD50CIPivot[0], TD50CIPivot[1]))
              
            #fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharex=False, sharey=False, figsize=(20,8))
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
            if len(self.TD50Hist):
                self.log("Found prior TD50 distribution, making differential histogram:")
                TD50DiffHist = [TD50A - TD50B for TD50A, TD50B in zip(self.TD50Hist, TD50Hist)]
                ax6.hist(TD50DiffHist, bins=50)
                ax6.set_xlabel("Differential TD50 distribution")
                diffMean = np.mean(TD50DiffHist)
                diffMedian = np.mean(TD50DiffHist)                

                oldInit = oldBestParametersNone[2]
                thisInit = self.bestParametersNone[2]
                diffInit = oldInit - thisInit
                
                print("oldInit = {:.2f}, thisInit = {:.2f}".format(oldInit, thisInit))

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

                self.log("The differential TD50 distribution (A - B) is {:.2f} ({:.2f} - {:.2f})".format(bestParameterDiff, diffCI[0], diffCI[1]))
                
            self.TD50Hist = TD50Hist
            plt.show()

        # Write to output files
        bootstrapOutput = open("Output/bs.csv")
        for idx in range(len(self.mHist)):
            bootstrapOutput.write(f"{mHist[idx]:.4f},{nHist[idx]:.4f},{TD50Hist[idx]:.4f}\n")
        bootstrapOutput.close()
            
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
        self.lastN = None
        self.lastGEUD = None
        self.GEUDlist = None
        self.nList = None

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
        
    def setMeanDoseFromEclipse(self, d):
        self.meanDoseFromEclipse = d

    def setMinDoseFromEclipse(self, d):
        self.minDoseFromEclipse = d
        
    def setMaxDoseFromEclipse(self, d):
        self.maxDoseFromEclipse = d
        
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
        nrows = self.dvh.shape[0]
        avgDoseList = np.zeros(nrows)
        diffVolumeList = np.zeros(nrows)
        doses = self.dvh["Dose"]
        volumes = self.dvh["Volume"]
        
        maxVolume = volumes.at[0]
        volumes = volumes / maxVolume
    
        for idx in range(nrows):
            try:    avgDose = np.mean([doses[idx], doses[idx+1]])
            except: avgDose = doses[idx]
            try:    diffVolume = volumes[idx] - volumes[idx+1]
            except: diffVolume = 0
            
            avgDoseList[idx] = avgDose
            diffVolumeList[idx] = diffVolume

        # Add the two new columns (avg dose, differential volume) to the dataframe
        self.dvh = self.dvh.assign(avgDose = avgDoseList, diffVolume = diffVolumeList)
        self.dvh.columns = ["Dose", "Volume", "Avg. Dose", "Diff. Volume"]
        
        # Remove rows where the differential volume is zero, not needed for calculation, reduce data set by 10%-20%
        self.dvh = self.dvh[self.dvh["Diff. Volume"] > 0]
        self.dvh = self.dvh.reset_index()
    
    def checkGEUDsplines(self, options):
        try:
            with open("%s/gEUD/gEUD_%s_%s_%s_%s.csv" % (self.dataFolder, self.cohort, self.structure, options.getCustomHeaderIdx(), self.ID), "r") as fh:
                GEUDlist = []               
                nList = []
                for line in fh.readlines():
                    linesplit = line.split(",")
                    nList.append(float(linesplit[0]))
                    GEUDlist.append(float(linesplit[1]))
                self.nList = np.array(nList)
                self.GEUDlist = np.array(GEUDlist)
                options.nFrom.set(nList[0])
                options.nTo.set(nList[-1]-0.02)
                
            return True
        except:
            return False
    
    def createGEUDspline(self, options):
        def integrate_eud(row):
            return row["Avg. Dose"] ** ninv * row["Diff. Volume"]
            
        if not self.checkGEUDsplines(options):
            #nList = np.arange(options.nFrom.get(), options.nTo.get()+0.04, 0.02)
            nList = np.arange(0.02, 1.2, 0.02)
            GEUDlist = []
            for n in nList:
                ninv = 1/n
                EUD = self.dvh.apply(integrate_eud, axis=1)
                GEUD = EUD.sum() ** n
                GEUDlist.append(GEUD)
            try:
                self.nList = np.array(nList)
                self.GEUDlist = np.array(GEUDlist)
            except:
                print("Could not create gEUD list for patient {}.".format(self.ID))
                self.nList = None
                self.GEUDlist = None
            
            # Save result to CSV files
            if not os.path.exists("%s/gEUD/" % (self.dataFolder)): 
                os.makedirs("%s/gEUD/" % (self.dataFolder))
            with open("%s/gEUD/gEUD_%s_%s_%s_%s.csv" % (self.dataFolder, self.cohort, self.structure, options.getCustomHeaderIdx(), self.ID), "wb") as fh:
                for n, GEUD in zip(nList, GEUDlist): 
                    fh.write(b"%.3f,%.3f\n" % (n, GEUD))
                    
    def getGEUD(self, n):
        if self.lastN == n:
            return self.lastGEUD
            
        else:
            index = bisect_left(self.nList, n)
            _xrange = self.nList[index] - self.nList[index-1]
            xdiff = n - self.nList[index-1]
            modolo = xdiff/_xrange
            ydiff = self.GEUDlist[index] - self.GEUDlist[index-1]
            self.lastN = n
            self.lastGEUD = self.GEUDlist[index-1] + modolo*ydiff
        
        return self.lastGEUD
        
    def getDoseAtVolume(self, volume):
        dvh_sort = self.dvh.ix[self.dvh["Volume"].argsort()]
        return np.interp(volume, dvh_sort["Volume"], dvh_sort["Dose"])
        
    def getVolumeAtDose(self, dose):
        if self.dvh.ix[self.dvh["Dose"] > dose]["Volume"].sum() == 0:
            return 0        
        return np.interp(dose, self.dvh["Dose"], self.dvh["Volume"])
        
    def calculateDpercent(self, percent):
        Dpercent = self.getDoseAtVolume(percent)
        self.Dpercent = Dpercent

    def getDpercent(self):
        return self.Dpercent

class Patients:
    def __init__(self, options):
        self.patients = {}
        self.bestParameters = []
        self.dataFolder = None
        self.structure = None
        self.cohort = None
        self.doseUnit = None
        self.options = options

        self.plot = None
        self.ax1 = None
        self.style1 = ["darkred", "darkblue", "k"]
        self.style2 = ["r", "b", "k"]

    def setDataFolder(self, folder):
        self.dataFolder = folder

    def getFilePath(self, fileName):
        return self.dataFolder + "/" + fileName

    def findStructures(self):
        structureNames = []
        for root, d, f, in os.walk(self.dataFolder):
            for filename in f:
                if not ".txt" in filename[-4:]: continue
                if "gEUD" in filename: continue                
                with open(self.getFilePath(filename), "r") as textin:
                    for line in textin:
                        if "Structure:" in line:
                            structureNames.append(re.sub("\s+", "", line.split(": ")[-1]))
        return list(set(structureNames))
    
    def findPlans(self):
        planNames = []
        for root, d, f, in os.walk(self.dataFolder):
            for filename in f:
                if not ".txt" in filename[-4:]: continue
                if "gEUD" in filename: continue
                try:
                    with open(self.getFilePath(filename), "r", encoding="utf8") as textin:
                        for line in textin:
                            if "Plan:" in line:
                                planNames.append(re.sub("\s+", "", line.split(": ")[-1]))
                except:
                    print(f"Skipping filename {filename}")
                    continue
        return list(set(planNames))

    def findNPatients(self):
        nPatients = 0
        for root, d, f in os.walk(self.dataFolder):
            for filename in f:
                if not ".txt" in filename[-4:]: continue
                if "gEUD" in filename: continue
                nPatients +=1
        return nPatients

    def loadPatientsECLIPSE(self, progress):
        def match(a,b):
            a += "$" # Don't match end-of-lines
            a = a.replace("*", ".*") # Use regex type wildcard
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
                log.append("Cannot open toxicity file %s_tox.csv, please include in format (filename,tox)." % (self.dataFolder))
                return log
                
        n=0
        for root, d, f, in os.walk(self.dataFolder):
            for filename in f:
                if not ".txt" in filename[-4:]: continue
                if "gEUD" in filename: continue
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
                skipPlan = 3 # one extra "Plan: " per file

                with open(self.getFilePath(filename), "r") as textin:
                    for line in textin:
                        if "Plan: " in line and skipPlan:
                            skipPlan -= 1
                            idx += 1
                            continue
                        
                        if "Structure: " in line:
                            structureNames.append(re.sub("\s+", "", line.split(": ")[-1]))
                            if structureStarts:
                                structureLength.append(idx-structureStarts[-1]-1)
                        
                        if "Plan: " in line:
                            planNames.append(re.sub("\s+", "", line.split(": ")[-1]))
                            if planStarts:
                                planLength.append(idx-planStarts[-1]-3)

                        if "Min Dose [" in line:
                            structureMinDose.append(float(line.split(": ")[-1]))

                        if "Max Dose [" in line:
                            structureMaxDose.append(float(line.split(": ")[-1]))
                        
                        if "Mean Dose [" in line:
                            structureMeanDose.append(float(line.split(": ")[-1]))
                                
                        if "Dose [" in line and "Volume [" in line:
                            structureStarts.append(idx+1)
                            planStarts.append(idx+1)
                        
                        idx += 1

                # To get to end of file
                structureLength.append(1000000)
                planLength.append(1000000)
                
                structures = zip(structureNames, structureStarts, structureLength, structureMeanDose, structureMinDose, structureMaxDose)
                plans = zip(planNames, planStarts, planLength)
                
                for structure, plan in zip(structures, plans):
                    if match(self.options.structureToUse.get(), structure[0]) and match(self.options.planToUse.get(), plan[0]):
                        headers = self.options.customDVHHeader.get().split(",")
                        headers = ["Dose", "Relative dose", "Volume"]

                        dvh = pd.read_csv(self.getFilePath(filename), header=None, names = headers,
                                                usecols = ["Dose", "Volume"], decimal=".", sep="\s+",
                                                skiprows=structure[1], nrows = structure[2], engine="python")

                        if np.sum(dvh.isnull().values):
                            headers = ["Dose", "Volume"]
                            dvh = pd.read_csv(self.getFilePath(filename), header=None, names = headers,
                                            usecols = ["Dose", "Volume"], decimal=".", sep="\s+",
                                            skiprows=structure[1], nrows = structure[2], engine="python")

                        # Create Patient object with DVH data
                        if self.options.doseUnit.get() == 'autodetect' and not self.doseUnit:
                            maxDose = dvh["Dose"].max()
                            if maxDose > 1000:
                                doseUnit = cGy
                            else:
                                doseUnit = Gy
                            self.doseUnit = doseUnit
                        elif self.doseUnit:
                            doseUnit = self.doseUnit
                        else:
                            doseUnit = self.options.doseUnit.get() == "cGy" and cGy or Gy

                        # SORT BY DESCENDING VOLUME
                        #dvh = dvh.sort_values(by=["Volume"], ascending=False).reset_index(drop=True)
                        
                        progress.step(1)
                        progress.update_idletasks()
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
                        patient.dvh["Dose"] = patient.dvh["Dose"] * doseUnit
                        patient.setID("%s_%s_%s" % (patientName, patient.getPlan(), patient.getStructure()))
                        patient.setMeanDoseFromEclipse(structure[3] * doseUnit)
                        patient.setMinDoseFromEclipse(structure[4] * doseUnit)
                        patient.setMaxDoseFromEclipse(structure[5] * doseUnit)

                        # Add object to dictionary
                        self.cohort = self.dataFolder.split("/")[-1]
                        self.structure = structure[0]
                        self.patients[patient.getID()] = patient
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
                if not "gEUD" in filename: nFiles += 1
        progress['maximum'] = nFiles            
            
        for root, d, f, in os.walk(self.dataFolder): # GENERALIZE THIS IF MORE 'simple' FILES ARE TO BE USED ...
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
                            dec = ","; sep = ";"
                        else:
                            dec = "."; sep = ","
                    
                elif self.options.CSVStyle.get() == "periodComma":
                    dec = "."; sep = ","
                elif self.options.CSVStyle.get() == "commaSemicolon":
                    dec = ","; sep = ";"
                                              
                headers = self.options.customDVHHeader.get().split(",")
                dvh = pd.read_csv(self.getFilePath(filename), decimal=dec, sep=sep, names = headers, 
                                  skiprows=self.options.skipRows.get(), usecols = ["Dose", "Volume"], dtype=np.float64) # skiprows = 2 if name

                if self.options.doseUnit.get() == 'autodetect' and not self.doseUnit:
                    maxDose = dvh["Dose"].max()
                    if maxDose > 1000:
                        doseUnit = cGy
                    else:
                        doseUnit = Gy
                    self.doseUnit = doseUnit
                elif self.doseUnit:
                    doseUnit = self.doseUnit
                else:
                    doseUnit = self.options.doseUnit.get() == "cGy" and cGy or Gy
                    
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
            self.log("Could not find best parameters...")
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

    def calculateDpercent(self, percent):
        for patient in list(self.patients.values()):
            patient.calculateDpercent(percent)
            
    def calculateTDxFromLogit(self, percent, a = None, b = None):
        X_test  = np.linspace(0, 100, 300)
        if not a:
            a = self.bestParameters[0]
        if not b:
            b = self.bestParameters[1]
        NTCP = np.array([1 - 1 / (1 + exp(a + b*k)) for k in X_test])
        TDx = X_test[np.argmin(abs(NTCP-percent/100))]
        return TDx
        
    def doMatrixMinimization(self, extraPatients, progress):
        """Original method from Piero Fossati."""
        
        nStep = (self.options.nTo.get()-self.options.nFrom.get())/self.options.matrixSize.get()
        mStep = (self.options.mTo.get()-self.options.mFrom.get())/self.options.matrixSize.get()
        TD50Step = (self.options.TD50To.get()-self.options.TD50From.get())/self.options.matrixSize.get()

        nList = self.options.fixN.get() and [self.options.nFrom.get()] or np.arange(self.options.nFrom.get(), self.options.nTo.get()+nStep, nStep)
        mList = self.options.fixM.get() and [self.options.mFrom.get()] or np.arange(self.options.mFrom.get(), self.options.mTo.get()+mStep, mStep)
        TD50List = self.options.fixTD50.get() and [self.options.TD50From.get()] or np.arange(self.options.TD50From.get(), self.options.TD50To.get()+TD50Step, TD50Step)

        toxArray = np.zeros((len(nList), len(mList), len(TD50List)), dtype=np.longdouble)

        if progress:
            progress['maximum'] = len(nList) * len(mList) * len(TD50List) / 5
        for nidx, n in enumerate(nList):
            if progress:
                progress.step(1)
                progress.update_idletasks()
                
            for midx, m in enumerate(mList):
                for TD50idx, TD50 in enumerate(TD50List):
                    error = 0
                    for patient in list(self.patients.values()):
                        gEUD = patient.getGEUD(n)
                        LKB = HPM((gEUD - TD50) / (m * TD50))
                        error += (patient.getTox() >= self.options.toxLimit.get() - LKB) ** 2
                    for cohort in extraPatients:
                        for patient in list(cohort.patients.values()):
                            gEUD = patient.getGEUD(n)
                            LKB = HPM((gEUD - TD50) / (m * TD50))
                            error += (patient.getTox() >= self.options.toxLimit.get() - LKB) ** 2
                    toxArray[nidx, midx, TD50idx] = error
                        
        minIdx = np.unravel_index(np.argmin(toxArray), np.shape(toxArray))
        
        self.bestParameters = [nList[minIdx[0]], mList[minIdx[1]], TD50List[minIdx[2]]]
        for cohort in extraPatients:
            cohort.bestParameters = self.bestParameters
            
        if progress: 
            progress['value'] = 0
        
        return "The best fit value is {:.2e}: n = {:.3f}{}, m = {:.3f}{} and TD50 = {:.1f} Gy.".format(
            toxArray[minIdx], nList[minIdx[0]], self.options.fixN.get() and " (fixed)" or "",
            mList[minIdx[1]], self.options.fixM.get() and " (fixed)" or "", TD50List[minIdx[2]])

    def doMatrixMinimizationLogistic(self, extraPatients, progress):
        
        aStep = (self.options.aTo.get()-self.options.aFrom.get())/self.options.matrixSize.get()
        bStep = (self.options.bTo.get()-self.options.bFrom.get())/self.options.matrixSize.get()
                
        aList = self.options.fixA.get() and [self.options.aFrom.get()] or np.arange(self.options.aFrom.get(), self.options.aTo.get()+aStep, aStep)        
        bList = self.options.fixB.get() and [self.options.bFrom.get()] or np.arange(self.options.bFrom.get(), self.options.bTo.get()+bStep, bStep)        

        toxArray = np.zeros((len(aList), len(bList)), dtype=np.longdouble)

        if progress:
            progress['maximum'] = len(aList)
            
        for aidx, a in enumerate(aList):
            if progress:
                progress.step(1)
                progress.update_idletasks()
                
            for bidx, b in enumerate(bList):
                error = 0
                for patient in list(self.patients.values()):
                    tox = patient.getTox() >= self.options.toxLimit.get()
                    NTCP = 1 - 1 / (1 + exp(a + b*patient.getDpercent()))
                    
                    if self.options.optimizationMetric.get() == "LS":
                        error += (tox - NTCP) ** 2
                    else:
                        if tox:
                            if NTCP>0:  error -= log(NTCP)
                            else:       error += 2500 # assume minimum probability of ~10^1000
                        else:
                            if NTCP<1:  error -= log(1-NTCP)
                            else:       error += 2500 # assume minimum probability of ~10^1000            
                            
                for cohort in extraPatients:
                    for patient in list(cohort.patients.values()):
                        NTCP = 1 - 1 / (1 + exp(a + b*patient.getDpercent()))
                        if self.options.optimizationMetric.get() == "LS":
                            error += (tox - NTCP) ** 2
                        else:
                            if tox:
                                if NTCP>0:  error -= log(NTCP)
                                else:       error += 2500 # assume minimum probability of ~10^1000
                            else:
                                if NTCP<1:  error -= log(1-NTCP)
                                else:       error += 2500 # assume minimum probability of ~10^1000                
                toxArray[aidx, bidx] = error

        minIdx = np.unravel_index(np.argmin(toxArray), np.shape(toxArray))
        
        self.bestParameters = [aList[minIdx[0]], bList[minIdx[1]]]
        for cohort in extraPatients:
            cohort.bestParameters = self.bestParameters
            
        if progress: 
            progress['value'] = 0
        
        return "The best fit value is {:.2e}: n = {:.3f}{}, m = {:.3f}{}.".format(
            toxArray[minIdx], aList[minIdx[0]], self.options.fixA.get() and " (fixed)" or "",
            bList[minIdx[1]], self.options.fixB.get() and " (fixed)" or "")

    def doOptimizationLogistic(self, extraPatients, progress):
        class MyTakeStep(object):
            def __init__(self, options, stepsize=1):
                self.stepsize = stepsize
                self.options = options
    
            def __call__(self, x):
                s = self.stepsize
                x[0] += np.random.uniform(-s*self.options.basinHoppingAsize.get(), s*self.options.basinHoppingAsize.get())
                x[1] += np.random.uniform(-s*self.options.basinHoppingBsize.get(), s*self.options.basinHoppingBsize.get())
                return x
            
        def funLS(x, *args):
            error = 0
            for tox, Dpercent in args:
                NTCP = 1 - 1 / (1 + exp(x[0] + x[1]*Dpercent))
                error += (tox - NTCP) ** 2
            return error
            
        def funLLH(x, *args):
            error = 0
            for tox, Dpercent in args:
                NTCP = 1 - 1 / (1 + exp(x[0] + x[1]*Dpercent))

                if tox:
                    if NTCP>0:  error -= log(NTCP)
                    else:       error += 2500 # assume minimum probability of ~10^1000
                else:
                    if NTCP<1:  error -= log(1-NTCP)
                    else:       error += 2500 # assume minimum probability of ~10^1000
            return error
        
        def print_fun(x,f,accepted):
            if progress:
                progress.step(1)
                progress.update_idletasks()
                
        if progress:
            progress['maximum'] = self.options.basinHoppingIterations.get()
                
        argTuple = ()
        for name, patient in list(self.patients.items()):
            argTuple += ((patient.getTox() >= self.options.toxLimit.get(), patient.getDpercent()),)
        for cohort in extraPatients:
            for name, patient in list(cohort.patients.items()):
                argTuple += ((patient.getTox() >= self.options.toxLimit.get(), patient.getDpercent()),)
            
        mytakestep = MyTakeStep(self.options)
        bounds = ((self.options.aFrom.get(), self.options.fixA.get() and self.options.aFrom.get() or self.options.aTo.get()),
                  (self.options.bFrom.get(), self.options.fixB.get() and self.options.bFrom.get() or self.options.bTo.get()))

        if len(self.bestParameters):
            x0 = np.array(self.bestParameters[:2])
        else:
            x0 = np.array([-10, 0.2])
        
        if self.options.optimizationMetric.get() == "LLH": fun = funLLH
        elif self.options.optimizationMetric.get() == "LS": fun = funLS
        
        res = basinhopping(fun, x0, niter=self.options.basinHoppingIterations.get(), T=self.options.basinHoppingTemperature.get(), 
                           minimizer_kwargs={'args':argTuple, 'method':'TNC', 'bounds':bounds}, take_step=mytakestep, callback=print_fun)
        
        self.bestParameters = [res.x[0], res.x[1]]
        
        for cohort in extraPatients:
            cohort.bestParameters = self.bestParameters
            
        if progress:
            progress['value'] = 0
            
        res["TD5%"]  = self.calculateTDxFromLogit(5)
        res["TD50%"] = self.calculateTDxFromLogit(50)
        
        return res
        
    def doOptimization(self, extraPatients, progress):
        """See http://stacks.iop.org/1742-6596/489/i=1/a=012087?key=crossref.181b59106e0d253de74e704220e16c36 for method."""
        
        class MyTakeStep(object):
            def __init__(self, options, stepsize=1):
                self.stepsize = stepsize
                self.options = options
    
            def __call__(self, x):
                s = self.stepsize
                x[0] += np.random.uniform(-s*self.options.basinHoppingNsize.get(), s*self.options.basinHoppingNsize.get())
                x[1] += np.random.uniform(-s*self.options.basinHoppingMsize.get(), s*self.options.basinHoppingMsize.get())
                x[2] += np.random.uniform(-s*self.options.basinHoppingTD50size.get(), s*self.options.basinHoppingTD50size.get())
                return x
        
        def funLS(x, *args):
            error = 0
            n = x[0]
            m = x[1]
            TD50 = x[2]
            for tox, GEUDspline in args:
                LKB = HPM((GEUDspline(n) - TD50) / (m*TD50))
                error += (tox - LKB) ** 2
            return error
    
        def funLLH(x, *args):
            error = 0
            for tox, GEUDspline in args:
                LKB = HPM((GEUDspline(x[0]) - x[2]) / (x[1]*x[2]))
                if tox:
                    if LKB>0: error -= log(LKB)
                    else:     error += 2500 # assume minimum probability of 10^-1000
                elif not tox:
                    if LKB<1: error -= log(1-LKB)
                    else:     error += 2500 # assume minimum probability of 10^-1000
            return error
        
        def print_fun(x, f, accepted):
            if progress: 
                progress.step(1)
                progress.update_idletasks()
        
        if progress:
            progress['maximum'] = self.options.basinHoppingIterations.get()
        
        argTuple = ()
        for name, patient in list(self.patients.items()):
            argTuple += ((patient.getTox() >= self.options.toxLimit.get(), patient.getGEUD),)
        for cohort in extraPatients:
            for name, patient in list(cohort.patients.items()):
                argTuple += ((patient.getTox() >= self.options.toxLimit.get(), patient.getGEUD),)
            
        mytakestep = MyTakeStep(self.options)
        bounds = ((self.options.nFrom.get(), self.options.fixN.get() and self.options.nFrom.get() or self.options.nTo.get()), 
                  (self.options.mFrom.get(), self.options.fixM.get() and self.options.mFrom.get() or self.options.mTo.get()),
                  (self.options.TD50From.get(), self.options.fixTD50.get() and self.options.TD50From.get() or self.options.TD50To.get()))
                  
        if len(self.bestParameters):
            x0 = np.array(self.bestParameters)
            if (len(x0)) == 3: x0[2] = 50
        else:
            x0 = np.array([0.2, 0.5, 50])
        
        if self.options.optimizationMetric.get() == "LLH": fun = funLLH
        elif self.options.optimizationMetric.get() == "LS": fun = funLS
        
        res = basinhopping(fun, x0, niter=self.options.basinHoppingIterations.get(), T=self.options.basinHoppingTemperature.get(), 
                           minimizer_kwargs={'args':argTuple, 'method':'TNC', 'bounds':bounds}, take_step=mytakestep, callback=print_fun)
        
        self.bestParameters = [res.x[0], res.x[1], res.x[2]]
        
        for cohort in extraPatients:
            cohort.bestParameters = self.bestParameters
            
        if progress:
            progress['value'] = 0
        
        return res
    
    def drawSigmoid(self, extraPatients, confidenceInterval, correlationLogit, aHist, bHist, TD50Hist, LLHhist, log):
        """Plot the patient outcomes vs sigmoid probability from the best parameter set."""

        D5_lower = None
        D50_lower = None
        D5_upper = None
        D50_upper = None

        if self.options.NTCPcalculation.get() == "LKB":
            n = self.bestParameters[0]
            m = self.bestParameters[1]
            TD50 = self.bestParameters[2]
            print("n = {}, m = {}, TD50 = {}".format(n,m,TD50))
            if np.sum(np.ravel(confidenceInterval)):
                nMin = confidenceInterval[0][0]
                nMax = confidenceInterval[0][1]
                mMin = confidenceInterval[1][0]
                mMax = confidenceInterval[1][1]
                TD50Min = confidenceInterval[2][0]
                TD50Max = confidenceInterval[2][1]
        
        else:
            a = self.bestParameters[0]
            b = self.bestParameters[1]
            print("a = {}, b = {}".format(a,b))
            if np.sum(np.ravel(confidenceInterval)):
                aMin = confidenceInterval[0][0]
                aMax = confidenceInterval[0][1]
                bMin = confidenceInterval[1][0]
                bMax = confidenceInterval[1][1]
                TD50Min = confidenceInterval[2][0]
                TD50Max = confidenceInterval[2][1]
        
        x = np.arange(0, 120, 0.1)
        yExtra = []
        
        if self.options.NTCPcalculation.get() == "LKB":
            y = np.fromiter([HPM(( k - TD50 ) / (m * TD50)) for k in x], float)
            if np.sum(np.ravel(confidenceInterval)):
                for n_, m_, TD50_ in zip(aHist, bHist, TD50Hist):
                    if (self.options.fixN.get() or nMin < n_ < nMax) and \
                        (self.options.fixM.get() or mMin < m_ < mMax) and \
                        (self.options.fixTD50.get() or TD50Min < TD50_ < TD50Max):
                        yExtra.append(np.fromiter([HPM(( k - TD50_ ) / (m_ * TD50_)) for k in x], float))
            
                # EMPIRICAL CONFIDENCE LIMITS (EXTREMITIES OF ALL MODELS WITHIN PARAMETER CIs)
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
            y = np.fromiter([1 - 1 / (1 + exp(a + b*k)) for k in x], float)
            if np.sum(np.ravel(confidenceInterval)):
                for a_,b_,TD50_ in zip(aHist, bHist,TD50Hist):
                    #if (self.options.fixA.get() or aMin < a_ < aMax) and \
                        #(self.options.fixB.get() or bMin < b_ < bMax):
                    if TD50Min < TD50_ < TD50Max:
                        yExtra.append(np.fromiter([1 - 1/(1 + exp(a_+b_*k)) for k in x], float))
                
                # EMPIRICAL CONFIDENCE LIMITS (EXTREMITIES OF ALL MODELS WITHIN PARAMETER CIs)
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
        for cohort in extraPatients:
            size += len(cohort.patients)
        
        px = np.zeros(size)
        py = np.zeros(size)

        idx = 0
        
        for patient in list(self.patients.values()):
            if self.options.NTCPcalculation.get() == "LKB":
                px[idx] = patient.getGEUD(n)
            else:
                px[idx] = patient.getDpercent()
            py[idx] = patient.getTox() >= self.options.toxLimit.get()
            idx += 1
            
        for cohort in extraPatients:
            for patient in list(cohort.patients.values()):
                if self.options.NTCPcalculation.get() == "LKB":
                    px[idx] = patient.getGEUD(n)
                else:
                    px[idx] = patient.getDpercent()
                py[idx] = patient.getTox() >= self.options.toxLimit.get()
                idx += 1

        # Don't make new if already... But this is secondary usage
        self.plot = plt.figure(figsize=(15, 10))
        self.ax1 = self.plot.add_subplot(111)
            
        style1 = self.style1.pop(0)
        style2 = self.style2.pop(0)

        self.ax1.plot(x,y, "-", color=style1, zorder=15, linewidth=3)
        self.ax1.plot(px,py, "o", color=style2, zorder=0)
        if self.options.confidenceIntervalShowModels.get():
            for each in yExtra:
                self.ax1.plot(x, each, '-', color="black", linewidth = 0.25, zorder=0)
                
        if np.sum(np.ravel(confidenceInterval)):
            self.ax1.fill_between(x, yMinEmpirical, yMaxEmpirical, color=style2, alpha=0.3, label="{:.0f}% confidence interval".format(self.options.confidenceIntervalPercent.get()), zorder=10)
            
        plt.ylim([-0.1,1.1])
        if self.options.NTCPcalculation.get() == "LKB":
            plt.xlabel("gEUD [Gy]")
            if D50_lower:
                plt.title("LKB for {}{}; n = {:.3f} ({:.3f}-{:.3f}), m = {:.3f} ({:.3f}-{:.3f}), TD50 = {:.2f} ({:.2f}-{:.2f}) Gy.".format(
                    self.cohort,np.sum(np.ravel(confidenceInterval)) and " with {}% CI".format(self.options.confidenceIntervalPercent.get()) or "", 
                    n, nMin, nMax, m, mMin, mMax, TD50, D50_lower, D50_upper))
            else:
                 plt.title("LKB for {}{}; n = {:.3f}, m = {:.3f}, TD50 = {:.2f} Gy.".format(
                    self.cohort,np.sum(np.ravel(confidenceInterval)) and " with {}% CI".format(self.options.confidenceIntervalPercent.get()) or "", 
                    n, m, TD50))
                
        else:
            plt.xlabel("D{}% [Gy]".format(self.options.NTCPcalculationDpercent.get()))
            if not D50_lower:
                plt.title("Logit for {}{}, using D{:.0f}%. TD5 = {:.1f} Gy, TD50 = {:.1f} Gy.".format(
                            self.cohort, np.sum(np.ravel(confidenceInterval)) and " with {}% CI".format(self.options.confidenceIntervalPercent.get()) or "",
                            self.options.NTCPcalculationDpercent.get(), self.calculateTDxFromLogit(5), self.calculateTDxFromLogit(50)))
            else:
                plt.title("Logit for {}{}, using D{:.0f}%. TD5 = {:.1f} ({:.1f}-{:.1f}) Gy, TD50 = {:.1f} ({:.1f}-{:.1f}) Gy.".format(
                        self.cohort, np.sum(np.ravel(confidenceInterval)) and " with {}% CI".format(self.options.confidenceIntervalPercent.get()) or "",
                        self.options.NTCPcalculationDpercent.get(), self.calculateTDxFromLogit(5), D5_lower, D5_upper, self.calculateTDxFromLogit(50), D50_lower, D50_upper))
                log("TD5 = {:.1f} ({:.1f}-{:.1f}) Gy ({}% CI)".format(
                        self.calculateTDxFromLogit(5), D5_lower, D5_upper, self.options.confidenceIntervalPercent.get()))
                log("TD50 = {:.1f} ({:.1f}-{:.1f}) Gy ({}% CI)".format(
                        self.calculateTDxFromLogit(50), D50_lower, D50_upper, self.options.confidenceIntervalPercent.get()))
        plt.ylabel("Toxicity / probability")
        
        plt.show()
        plt.savefig("Output/LKBgraph_%s_%s.png" % (self.cohort, self.structure)) # VIRKER IKKE !!! 
        
        if self.options.NTCPcalculation.get() == "LKB":
            of = open("Output/LKBmodel_%s_%s.csv" % (self.cohort, self.structure), "w") 
            of.write("n = {:.3f}, m = {:.3f}, TD50 = {:.2f} Gy.\n\n".format(n, m, TD50))
            of.write("gEUD tox from patients\n")
            of.write("\n".join([" ".join([str(k[0]), str(k[1])]) for k in zip(px, py)]))
            of.write("\n\ngEUD LKB from model\n")
            of.write("\n".join([" ".join([str(k[0]), str(k[1])]) for k in zip(x,y)]))
            of.close()
        else:
            of = open("Output/LogitModel_%s_%s.csv" % (self.cohort, self.structure), "w") 
            of.write("a = {:.3f}, b = {:.3f}.\n\n".format(a, b))
            if np.sum(np.ravel(confidenceInterval)):
                of.write("{}% CI: a = {:.3f} - {:.3f}, b = {:.3f} - {:.3f}\n".format(self.options.confidenceIntervalPercent.get(),
                confidenceInterval[0][0], confidenceInterval[0][1], confidenceInterval[1][1], confidenceInterval[1][1]))
            of.write("gEUD tox from patients\n")
            of.write("\n".join([" ".join([str(k[0]), str(k[1])]) for k in zip(px, py)]))
            of.write("\n\ngEUD LKB from model\n")
            of.write("\n".join([" ".join([str(k[0]), str(k[1])]) for k in zip(x,y)]))
            of.close()
            
    def profileLikelihood(self, extraPatients):
        """For some reason this article gets cited for this method:
            Stryhn, H, and J Christensen. “Confidence Intervals by the Profile Likelihood Method, 
            with Applications in Veterinary Epidemiology,” 2003. http://www.sciquest.org.nz/elibrary/edition/5008."""
        q = 1 - self.options.confidenceIntervalPercent.get()/100
        gamma = st.chi2.isf(q, df=1) / 2
        
        if self.options.NTCPcalculation.get() == "Dpercent":
            aInit = self.bestParameters[0]
            aValues = np.arange(aInit*0.5, aInit*2, aInit/100)
            bInit = self.bestParameters[1]
            bValues = np.arange(bInit*0.5, bInit*2, bInit/100)
            
            LLH_a = np.zeros(len(aValues))
            LLH_b = np.zeros(len(bValues))

            for idx, a in enumerate(aValues): 
                error = 0
                for patient in list(self.patients.values()):
                    tox = patient.getTox() >= self.options.toxLimit.get() 
                    NTCP = 1 - 1 / (1 + exp(a + bInit*patient.getDpercent()))
                    if tox:
                        if NTCP>0: error -= log(NTCP)
                        else:     error += 2500 # assume minimum probability of 10^-1000
                    elif not tox:
                        if NTCP<1: error -= log(1-NTCP)
                        else:     error += 2500 # assume minimum probability of 10^-1000
                        
                for cohort in extraPatients:
                    for patient in list(cohort.patients.values()):
                        tox = patient.getTox() >= self.options.toxLimit.get() 
                        NTCP = 1 - 1 / (1 + exp(a + bInit*patient.getDpercent()))
                        if tox:
                            if NTCP>0: error -= log(NTCP)
                            else:     error += 2500 # assume minimum probability of 10^-1000
                        elif not tox:
                            if NTCP<1: error -= log(1-NTCP)
                            else:     error += 2500 # assume minimum probability of 10^-1000
                LLH_a[idx] = -error           
                
            aCIvalue = np.max(LLH_a) - gamma
            idx = np.argwhere(np.diff(np.sign(LLH_a - aCIvalue))).flatten()
            aLow = min(aValues[idx[0]], aValues[idx[1]])
            aHigh = max(aValues[idx[0]], aValues[idx[1]])
            
            for idx, b in enumerate(bValues): 
                error = 0
                for patient in list(self.patients.values()):
                    tox = patient.getTox() >= self.options.toxLimit.get() 
                    NTCP = 1 - 1 / (1 + exp(aInit + b*patient.getDpercent()))
                    if tox:
                        if NTCP>0: error -= log(NTCP)
                        else:     error += 2500 # assume minimum probability of 10^-1000
                    elif not tox:
                        if NTCP<1: error -= log(1-NTCP)
                        else:     error += 2500 # assume minimum probability of 10^-1000
                        
                for cohort in extraPatients:
                    for patient in list(cohort.patients.values()):
                        tox = patient.getTox() >= self.options.toxLimit.get() 
                        NTCP = 1 - 1 / (1 + exp(a + bInit*patient.getDpercent()))
                        if tox:
                            if NTCP>0: error -= log(NTCP)
                            else:     error += 2500 # assume minimum probability of 10^-1000
                        elif not tox:
                            if NTCP<1: error -= log(1-NTCP)
                            else:     error += 2500 # assume minimum probability of 10^-1000
                LLH_b[idx] = -error

            bCIvalue = np.max(LLH_b) - gamma
            idx = np.argwhere(np.diff(np.sign(LLH_b - bCIvalue))).flatten()
            bLow = min(bValues[idx[0]], bValues[idx[1]])
            bHigh = max(bValues[idx[0]], bValues[idx[1]])

            fig, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(20,8))
            fig.add_subplot(111, frameon=False)
            plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')            
            plt.title("{}% profile likelihood scan over NTCP = 1 - 1 / (1 + exp(a + b * D{}%))".format(self.options.confidenceIntervalPercent.get(), self.options.NTCPcalculationDpercent.get()))            
            ax1.plot(aValues, LLH_a, label="Log Likelihood")
            ax1.plot([aValues[0], aValues[-1]], [aCIvalue, aCIvalue], "-", label="chi2(1)/2 = {:.2f} reduction in LLH".format(gamma))            
            ax1.set_xlabel("a values")
            ax1.set_ylabel("log likelihood")
            #ax1.set_xscale('log')
            ax1.legend(loc='lower right')
            ax2.plot(bValues, LLH_b, label="Log likelihood")
            ax2.plot([bValues[0], bValues[-1]], [bCIvalue, bCIvalue], "-", label="chi2(1)/2 = {:.2f} reduction in LLH".format(gamma))
            ax2.set_xlabel("b values")
            ax2.set_ylabel("log likelihood")
            #ax2.set_xscale('log')
            ax2.legend(loc='lower right')
            
            plt.savefig("Output/profileLikelihood.png")

            return [[aLow, aHigh], [bLow, bHigh]]
            
        if self.options.NTCPcalculation.get() == "LKB":
            nInit = self.bestParameters[0]
            nValues = np.arange(self.options.nFrom.get(), self.options.nTo.get() ,nInit/50)
            mInit = self.bestParameters[1]
            mValues = np.arange(mInit*0.2, mInit*3, mInit/100)
            TD50Init = self.bestParameters[2]
            TD50Values = np.arange(TD50Init*0.5, TD50Init*3, TD50Init/100)
            
            LLH_n = np.zeros(len(nValues))
            LLH_m = np.zeros(len(mValues))
            LLH_TD50 = np.zeros(len(TD50Values))
            
            for idx, n in enumerate(nValues): # scan over n
                error = 0
                for patient in list(self.patients.values()):
                    tox = patient.getTox() >= self.options.toxLimit.get() 
                    NTCP = HPM((patient.getGEUD(n) - TD50Init) / (mInit * TD50Init))
                    if tox:
                        if NTCP>0: error -= log(NTCP)
                        else:     error += 2500 # assume minimum probability of 10^-1000
                    elif not tox:
                        if NTCP<1: error -= log(1-NTCP)
                        else:     error += 2500 # assume minimum probability of 10^-1000
                        
                for cohort in extraPatients:
                    for patient in list(cohort.patients.values()):
                        tox = patient.getTox() >= self.options.toxLimit.get() 
                        NTCP = HPM((patient.getGEUD(n) - TD50Init) / (mInit * TD50Init))
                        if tox:
                            if NTCP>0: error -= log(NTCP)
                            else:     error += 2500 # assume minimum probability of 10^-1000
                        elif not tox:
                            if NTCP<1: error -= log(1-NTCP)
                            else:     error += 2500 # assume minimum probability of 10^-1000
                LLH_n[idx] = -error           
            
            nCIvalue = np.max(LLH_n) - gamma          
            idx = np.argwhere(np.diff(np.sign(LLH_n - nCIvalue))).flatten()
            nLow = nValues[idx[0]]
            nHigh = nValues[idx[1]]

            for idx, m in enumerate(mValues): # Scan over m
                error = 0
                for patient in list(self.patients.values()):
                    tox = patient.getTox() >= self.options.toxLimit.get() 
                    NTCP = HPM((patient.getGEUD(nInit) - TD50Init) / (m * TD50Init))
                    if tox:
                        if NTCP>0: error -= log(NTCP)
                        else:     error += 2500 # assume minimum probability of 10^-1000
                    elif not tox:
                        if NTCP<1: error -= log(1-NTCP)
                        else:     error += 2500 # assume minimum probability of 10^-1000
                        
                for cohort in extraPatients:
                    for patient in list(cohort.patients.values()):
                        tox = patient.getTox() >= self.options.toxLimit.get() 
                        NTCP = HPM((patient.getGEUD(nInit) - TD50Init) / (m * TD50Init))
                        if tox:
                            if NTCP>0: error -= log(NTCP)
                            else:     error += 2500 # assume minimum probability of 10^-1000
                        elif not tox:
                            if NTCP<1: error -= log(1-NTCP)
                            else:     error += 2500 # assume minimum probability of 10^-1000
                LLH_m[idx] = -error           
                
            mCIvalue = np.max(LLH_m) - gamma
            idx = np.argwhere(np.diff(np.sign(LLH_m - mCIvalue))).flatten()
            mLow = mValues[idx[0]]
            mHigh = mValues[idx[1]]
            
            for idx, TD50 in enumerate(TD50Values): # Scan over TD50
                error = 0
                for patient in list(self.patients.values()):
                    tox = patient.getTox() >= self.options.toxLimit.get() 
                    NTCP = HPM((patient.getGEUD(nInit) - TD50) / (mInit * TD50))
                    if tox:
                        if NTCP>0: error -= log(NTCP)
                        else:     error += 2500 # assume minimum probability of 10^-1000
                    elif not tox:
                        if NTCP<1: error -= log(1-NTCP)
                        else:     error += 2500 # assume minimum probability of 10^-1000
                        
                for cohort in extraPatients:
                    for patient in list(cohort.patients.values()):
                        tox = patient.getTox() >= self.options.toxLimit.get() 
                        NTCP = HPM((patient.getGEUD(nInit) - TD50) / (mInit * TD50))
                        if tox:
                            if NTCP>0: error -= log(NTCP)
                            else:     error += 2500 # assume minimum probability of 10^-1000
                        elif not tox:
                            if NTCP<1: error -= log(1-NTCP)
                            else:     error += 2500 # assume minimum probability of 10^-1000
                LLH_TD50[idx] = -error           
                
            TD50CIvalue = np.max(LLH_TD50) - gamma
            idx = np.argwhere(np.diff(np.sign(LLH_TD50 - TD50CIvalue))).flatten()
            TD50Low = TD50Values[idx[0]]
            TD50High = TD50Values[idx[1]]
            
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=False, sharey=False, figsize=(20,8))
            fig.add_subplot(111, frameon=False)
            plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')            
            plt.title("{}% profile likelihood scan over NTCP = LKB(n,m,TD50)".format(self.options.confidenceIntervalPercent.get()))
            ax1.plot(nValues, LLH_n, label="Log Likelihood")
            ax1.plot([nValues[0], nValues[-1]], [nCIvalue, nCIvalue], "-", label="chi2(1)/2 = {:.2f} reduction in LLH".format(gamma))            
            ax1.set_xlabel("n values")
            ax1.set_ylabel("log likelihood")
            ax1.legend(loc='lower right')
            ax2.plot(mValues, LLH_m, label="Log likelihood")
            ax2.plot([mValues[0], mValues[-1]], [mCIvalue, mCIvalue], "-", label="chi2(1)/2 = {:.2f} reduction in LLH".format(gamma))
            ax2.set_xlabel("m values")
            ax2.set_ylabel("log likelihood")
            ax2.legend(loc='lower right')
            ax3.plot(TD50Values, LLH_TD50, label="Log likelihood")
            ax3.plot([TD50Values[0], TD50Values[-1]], [TD50CIvalue, TD50CIvalue], "-", label="chi2(1)/2 = {:.2f} reduction in LLH".format(gamma))
            ax3.set_xlabel("TD50 values")
            ax3.set_ylabel("log likelihood")
            ax3.legend(loc='lower right')
            
            plt.savefig("Output/profileLikelihood.png")

            return [[nLow, nHigh], [mLow, mHigh], [TD50Low, TD50High]]
        
    def drawAUROC(self, extraPatients, progress):
        """Plot the AUROC curve for the chosen patient cohorts."""
        
        def auroc_error(theta, nplus, nminus):
            """From https://www.tandfonline.com/doi/full/10.1080/02841860903078513"""
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
            doseThresholds += [(gEUDlist[k] + gEUDlist[k+1]) / 2 for k in range(len(gEUDlist) - 1)]
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

                if (nN): x[idx] = fn/nN
                if (nP): y[idx] = tp/nP

            area = np.sum([(x[k] - x[k+1]) * y[k+1] for k in range(len(doseThresholds)-1)])
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
        plt.title("AUROC for {}. Max is {:.2f} at n={:.2f} +- {:.2f}.".format(self.cohort, maxval, xmean, xwidth))
        plt.ylabel('Area Under ROC')
        plt.xlabel('gEUD n-value')
        plt.show()
        progress['value'] = 0

"""
pr = cProfile.Profile()
pr.enable()
"""

root = Tk()
mainmenu = MainMenu(root)
root.mainloop()

"""
pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby).reverse_order()
ps.print_stats()
print(s.getvalue())
"""
