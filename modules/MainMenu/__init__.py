from tkinter import *
from tkinter import ttk

from ..Tools import *
from ..Patients import *


class MainMenu(Frame):
    def __init__(self, parent):
        from .._Options import Options
        from ._Tooltip import Tooltip
        self.options = Options()

        Frame.__init__(self, parent)
        self.parent = parent
        self.parent.protocol("WM_DELETE_WINDOW", self.myQuit)
        self.parent.title(f"DVH tool with NTCP calculator {self.options.PROGRAM_VERSION:.2f} - Helge Pettersen")
        self.window = None
        self.wraplength = 250
        self.button_width = 25

        self.bestParameters = []

        """
        self.bestParametersNone = []
        self.bestParametersMean = []
        self.bestParametersMedian = []
        self.confidenceInterval = [[0, 0], [0, 0], [0, 0]]
        self.confidenceIntervalNone = [[0, 0], [0, 0], [0, 0]]
        self.confidenceIntervalMean = [[0, 0], [0, 0], [0, 0]]
        self.confidenceIntervalMedian = [[0, 0], [0, 0], [0, 0]]
        self.correlationValueLogit = -0.016
        self.calculateMeanDose = IntVar(value=1)
        self.aHist = []
        self.bHist = []
        self.TD50Hist = []
        self.LLHhist = []
        """

        self.paramNEntry = []
        self.paramMEntry = []
        self.paramTD50Entry = []
        self.paramBHsizeLabel = []
        self.paramBHsizeEntry = []
        self.patients = {}
        self.cohortList = {}
        self.NTCPAxis = None
        self.useCustomAggregateDVHPlot = False
        self.style1 = ["darkred", "darkblue", "k"] * 100
        self.style2 = ["r", "b", "k"] * 100

        res = self.options.loadOptions()

        self.upperContainer = Frame(self, bd=5, relief=RIDGE, height=40)  # Title
        self.middleContainer = Frame(self, bd=5)
        self.bottomContainer = Frame(self, bd=20)

        self.middleLeftContainer = Frame(self.middleContainer, bd=5)  # Patient cohorts
        self.middleLeftUpperContainer = Frame(self.middleLeftContainer, bd=5)  # Load patient cohort
        self.middleLeftLine = Frame(self.middleLeftContainer, bg="grey", relief=SUNKEN)
        self.middleLeftLowerContainer = Frame(self.middleLeftContainer, bd=5)  # List patient cohorts
        self.middleMiddleLine = Frame(self.middleContainer, bg="grey", relief=SUNKEN)
        self.middleMiddleContainer = Frame(self.middleContainer, bd=5)  # Options
        self.middleRightLine = Frame(self.middleContainer, bg="grey", relief=SUNKEN)
        self.middleRightContainer = Frame(self.middleContainer, bd=5)  # Log
        self.middleRightContainerLine = Frame(self.middleRightContainer, bg="grey", relief=SUNKEN)
        self.middleRightLowerContainer = Frame(self.middleRightContainer)
        self.bottomLine = Frame(self.bottomContainer, bg="grey", relief=SUNKEN)
        self.bottomContainer1 = Frame(self.bottomContainer)
        self.bottomContainer2 = Frame(self.bottomContainer)

        self.dvhFileContainer = Frame(self.middleMiddleContainer)
        self.toxLimitContainer = Frame(self.middleMiddleContainer)
        self.toxFromFilenameContainer = Frame(self.middleMiddleContainer)
        self.changeNamingContainer = Frame(self.middleMiddleContainer)

        self.CSVStyleContainer = Frame(self.middleMiddleContainer)
        self.customDVHHeaderContainer = Frame(self.middleMiddleContainer)
        self.doseUnitContainer = Frame(self.middleMiddleContainer)
        self.skipRowsContainer = Frame(self.middleMiddleContainer)

        self.upperContainer.pack(fill=X)
        self.middleContainer.pack(anchor=N, fill=Y)
        self.middleLeftContainer.pack(side=LEFT, fill="both", expand=1, anchor=N)
        self.middleLeftUpperContainer.pack(fill=X)
        self.middleLeftLine.pack(fill=X, padx=5, pady=5)
        self.middleLeftLowerContainer.pack(fill=X)
        self.middleMiddleLine.pack(anchor=N, side=LEFT, fill=Y, padx=5, pady=5, expand=1)
        self.middleMiddleContainer.pack(anchor=N, side=LEFT, padx=5, pady=5)
        self.middleRightLine.pack(side=LEFT, fill=Y, padx=5, pady=5, expand=1)
        self.middleRightContainer.pack(side=LEFT, fill=Y)
        self.middleRightContainerLine.pack(fill=X, expand=1)
        self.middleRightLowerContainer.pack(fill=X, expand=1)
        self.bottomLine.pack(fill=X, padx=5, pady=5, expand=1)
        self.bottomContainer.pack(fill=X, anchor=N, expand=1)
        self.bottomContainer1.pack(anchor=N, expand=1)
        self.bottomContainer2.pack(anchor=N, expand=1)

        Label(self.upperContainer, text=f"DVH tool with NTCP calculator {self.options.PROGRAM_VERSION:.2f} - Helge Pettersen").pack(anchor=N)

        self.loadPatientsButton = Button(self.middleLeftUpperContainer, text="Load patient folder (F)",
                                         command=self.loadPatientsCommand, width=self.button_width)
        self.loadPatientsButton.pack(anchor=N, pady=3)
        Tooltip(self.loadPatientsButton, text="Load patient cohort. The file should contain DVH files in a CSV format. "
                "Upon loading, additional gEUD files will searched for in a gEUD/ subfolder. If they are not found, they should "
                "be created using the \"Calculate gEUD LUT\" action prior to any LKB calculations.", wraplength=self.wraplength)

        self.parent.bind("f", lambda event=None: self.loadPatientsButton.invoke())

        self.progress = ttk.Progressbar(self.middleLeftUpperContainer, orient=HORIZONTAL, maximum=100, mode='determinate')
        self.progress.pack(fill=X, pady=3)

        Label(self.middleMiddleContainer, text="FILE OPTIONS", font=("Helvetica", 12)).pack(anchor=N)
        Frame(self.middleMiddleContainer, bg="grey", relief=SUNKEN).pack(fill=X, expand=1, anchor=N)

        self.dvhFileContainer.pack(anchor=W)
        Label(self.dvhFileContainer, text="DVH File type: ").pack(side=LEFT, anchor=W)
        for mode in ["simple", "ECLIPSE", "RayStation"]:
            Radiobutton(self.dvhFileContainer, text=mode, variable=self.options.DVHFileType, value=mode,
                        command=self.selectDVHFileTypeCommand).pack(side=LEFT, anchor=W)
        Tooltip(self.dvhFileContainer, text="Outline of the per-patient CSV files. The \"simple\" style assumes a flat CSV file. "
                "The \"ECLIPSE\" and \"RayStation\" style includes several metadata lines in the beginning, incl. descriptions of the DVH organ / structure.",
                wraplength=self.wraplength)

        self.CSVStyleContainer.pack(anchor=W)
        Label(self.CSVStyleContainer, text="Decimal / Column separator: ").pack(side=LEFT, anchor=W)
        for text, mode in [[", / ;", "commaSemicolon"], [". / ,", "periodComma"], ["Autodetect", "autodetect"]]:
            Radiobutton(self.CSVStyleContainer, text=text, variable=self.options.CSVStyle, value=mode).pack(side=LEFT, anchor=W)
        Tooltip(self.CSVStyleContainer, text="CSV style: Comma and Semicolon, or Period and Comma for decimal and column separation.",
                wraplength=self.wraplength)

        self.customDVHHeaderContainer.pack(anchor=W)
        Label(self.customDVHHeaderContainer, text="Columns in DVH files: ").pack(side=LEFT, anchor=W)
        self.customDVHHeaderEntry = Entry(self.customDVHHeaderContainer, textvariable=self.options.customDVHHeader, width=50)
        self.customDVHHeaderEntry.pack(side=LEFT, anchor=W)
        Tooltip(self.customDVHHeaderContainer,
                text="Define headers as a comma-separated list, indicating the wanted Dose and Volume columns without whitespace. Example: "
                "\"Volume,dummy,dummy,Dose\" or \"Dose,Volume\".", wraplength=self.wraplength)

        self.doseUnitContainer.pack(anchor=W)
        Label(self.doseUnitContainer, text="Dose unit: ").pack(side=LEFT, anchor=W)
        for text, mode in [["cGy", "cGy"], ["Gy", "Gy"], ["Autodetect", "autodetect"]]:
            Radiobutton(self.doseUnitContainer, text=text, variable=self.options.doseUnit, value=mode).pack(side=LEFT, anchor=W)

        self.skipRowsContainer.pack(anchor=W)
        Label(self.skipRowsContainer, text="Number of rows to skip in simple csv: ").pack(side=LEFT, anchor=W)
        Entry(self.skipRowsContainer, textvariable=self.options.skipRows, width=5).pack(side=LEFT)

        Label(self.middleMiddleContainer, text="TOXICITY OPTIONS", font=("Helvetica", 12)).pack(anchor=N)
        Frame(self.middleMiddleContainer, bg="grey", relief=SUNKEN).pack(fill=X, expand=1, anchor=N)

        self.toxFromFilenameContainer.pack(anchor=W)
        Label(self.toxFromFilenameContainer, text="Load toxicity from: : ").pack(side=LEFT, anchor=W)
        for text, mode in [('\'tox\' in filename', 1), ('CSV_file ', 0)]:
            Radiobutton(self.toxFromFilenameContainer, text=text, variable=self.options.loadToxFromFilename,
                        value=mode, command=self.toxFromFilenameCommand).pack(side=LEFT, anchor=W)
        Tooltip(self.toxFromFilenameContainer, text="Where to get the toxicity information from. Either use the presence of \'tox\' in the filename "
                "to define a yes/no toxicity; or provide a CSV file \"Data/<folderName>_tox.csv\" with lines patientID,toxGrade. The patientID "
                "is the filename.", wraplength=self.wraplength)

        self.toxLimitContainer.pack(anchor=W)
        Label(self.toxLimitContainer, text="Toxicity limit (grade): ").pack(side=LEFT, anchor=W)
        for mode in range(5):
            Radiobutton(self.toxLimitContainer, text=mode, variable=self.options.toxLimit, value=mode,
                        command=self.toxLimitChange).pack(side=LEFT, anchor=W)
        Tooltip(self.toxLimitContainer,
                text="The toxicity threshold. Patients with the chosen grade or higher will be classified with complication.",
                wraplength=self.wraplength)

        self.changeNamingContainer.pack()
        self.changeNamingButton = Button(self.changeNamingContainer, text="Change plan/structure names",
                                         command=self.changeNamingCommand, width=self.button_width, state=DISABLED)
        self.changeNamingButton.pack()

        if self.options.DVHFileType.get() in ["ECLIPSE", "RayStation"]:
            self.customDVHHeaderEntry['state'] = 'disabled'

        self.logtext = Text(self.middleRightLowerContainer, height=25, width=75)
        self.scrollbar = Scrollbar(self.middleRightLowerContainer)
        self.scrollbar.config(command=self.logtext.yview)
        self.logtext.config(yscrollcommand=self.scrollbar.set)
        self.scrollbar.pack(side=RIGHT, fill="y", expand=False)
        self.logtext.pack(side=LEFT, fill="both", expand=True)
        self.log("----- LOG -----")

        self.buttonCalculateGEUD = Button(self.bottomContainer1, text="Calculate gEUD LUT (G)", command=self.calculateGEUDCommand, width=self.button_width, state=DISABLED)
        self.buttonShowGEUDvsN = Button(self.bottomContainer1, text="Show gEUD/n (H)", command=self.showGEUDvsN, width=self.button_width, state=DISABLED)
        self.buttonShowDVH = Button(self.bottomContainer1, text="Show DVHs (D)", command=self.showDVHCommand, width=self.button_width, state=DISABLED)
        self.buttonCalculateDVH = Button(self.bottomContainer1, text="Save DVH values (C)", command=self.calculateDVHCommand, width=self.button_width, state=DISABLED)
        self.buttonAggregateDVH = Button(self.bottomContainer1, text="Aggregate DVH plots (A)", command=self.aggregateDVHCommand, width=self.button_width, state=DISABLED)
        self.buttonCalculateNTCP = Button(self.bottomContainer2, text="Calculate NTCP model (N)", command=self.calculateNTCPWindow, width=self.button_width, state=DISABLED)
        self.buttonLKBuncert = Button(self.bottomContainer2, text="Calculate confidence intervals (B)", command=self.calculateBootstrapWindow, width=self.button_width, state=DISABLED)
        self.buttonCalculateAUROC = Button(self.bottomContainer2, text="Calculate AUROC (S)", command=self.calculateAUROCCommand, width=self.button_width, state=DISABLED)
        self.buttonQuit = Button(self.bottomContainer2, text="Exit (Esc)", command=self.myQuit, width=self.button_width)

        self.parent.bind("g", lambda event=None: self.buttonCalculateGEUD.invoke())
        self.parent.bind("d", lambda event=None: self.buttonShowDVH.invoke())
        self.parent.bind("c", lambda event=None: self.buttonCalculateDVH.invoke())
        self.parent.bind("a", lambda event=None: self.buttonAggregateDVH.invoke())
        self.parent.bind("n", lambda event=None: self.buttonCalculateNTCP.invoke())
        self.parent.bind("b", lambda event=None: self.buttonLKBuncert.invoke())
        self.parent.bind("s", lambda event=None: self.buttonCalculateAUROC.invoke())
        self.parent.bind("h", lambda event=None: self.buttonShowGEUDvsN.invoke())
        self.parent.bind("<Escape>", lambda event=None: self.buttonQuit.invoke())

        for button in [self.buttonCalculateGEUD, self.buttonShowGEUDvsN, self.buttonShowDVH, self.buttonCalculateDVH, self.buttonAggregateDVH,
                       self.buttonCalculateNTCP, self.buttonLKBuncert, self.buttonCalculateAUROC, self.buttonQuit]:
            button.pack(side=LEFT, anchor=N, padx=5, pady=5)

        self.pack()

        if not res:
            self.log("Could not load options.cfg, using default values.")

    # Imported modules
    from ._GUIelements import myQuit, log, loadPatientsCommand, addPatientCohort, removePatientCohort
    from ._GUIelements import chooseStructureCommand, calculateGEUDCommand, toxLimitChange, toxFromFilenameCommand
    from ._GUIelements import customDVHHeaderCommand, selectDVHFileTypeCommand, showDVHCommand, showDVHPlotCommand
    from ._GUIelements import bootstrapCorrectionMethodCommand, calculateNTCPCommand, NTCPcalculationCommand
    from ._GUIelements import calculateAUROCCommand, calculateDVHCommand, cancelCalculateDVHvalues, aggregateDVHCommand
    from ._GUIelements import showGEUDvsN, switchNto, switchMto, switchTD50to, changeNamingCommand, calculateNewNamesCommand
    from ._GUIelements import changeNamingQuitCommand, changeNamingQuitAndSaveCommand, drawPlanAndStructureNames
    from ._GUIelements import customAggregatedDVHCommand, cancelCustomAggregateDVHCommand, matchCustomAggregateDVHCommand
    from ._GUIelements import saveCustomAggregateDVHCommand, packCustomAggregateDVHCommand, calculateNTCPWindow, calculateNTCPWindowCancel, switchToNTCPcc
    from ._GUIelements import calculateBootstrapWindow, calculateBootstrapWindowCancel

    from ._Analysis import calculateDVHvalues, calculateAggregatedDVH
    from ._NTCPbootstrap import calculateLKBuncert
