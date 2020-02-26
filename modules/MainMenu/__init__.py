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
        self.patients = {}
        self.cohortList = {}
        self.NTCPAxis = None
        self.style1 = ["darkred", "darkblue", "k"]*100
        self.style2 = ["r", "b", "k"]*100
        
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
        self.changeNamingContainer = Frame(self.middleMiddleContainer)
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
        self.makeDifferentialCIContainer = Frame(self.middleMiddleContainer)
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
        Frame(self.middleMiddleContainer, bg="grey", relief=SUNKEN).pack(fill=X,expand=1,anchor=W)
        
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
        for text, mode in [ ["cGy", "cGy"], ["Gy", "Gy"], ["Autodetect", "autodetect"]]:
            Radiobutton(self.doseUnitContainer, text=text, variable=self.options.doseUnit, value=mode).pack(side=LEFT, anchor=W)
    
        self.skipRowsContainer.pack(anchor=W)
        Label(self.skipRowsContainer, text="Number of rows to skip in simple csv: ").pack(side=LEFT, anchor=W)
        Entry(self.skipRowsContainer, textvariable=self.options.skipRows, width=5).pack(side=LEFT)        
        
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

        Label(self.middleMiddleContainer, text="NTCP OPTIONS", font=("Helvetica", 12) ).pack(pady=(15,0))
        Frame(self.middleMiddleContainer, bg="grey", relief=SUNKEN).pack(fill=X,expand=1,anchor=W)
        
        self.NTCPcalculationContainer.pack(anchor=W)
        Label(self.NTCPcalculationContainer, text="NTCP calculation: ").pack(side=LEFT, anchor=W)
        for text, mode in [("D% + logit", "Logit"), ("LKB", "LKB")]:
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

        Label(self.middleMiddleContainer, text="BOOTSTRAP OPTIONS", font=("Helvetica", 12) ).pack(pady=(15,0))
        Frame(self.middleMiddleContainer, bg="grey", relief=SUNKEN).pack(fill=X,expand=1,anchor=W)

        self.confidenceIntervalPercentContainer.pack(anchor=W)
        Label(self.confidenceIntervalPercentContainer, text="Confidence Interval percentage: ").pack(side=LEFT, anchor=W)
        Entry(self.confidenceIntervalPercentContainer, textvariable=self.options.confidenceIntervalPercent, width=7).pack(side=LEFT, anchor=W)
        Label(self.confidenceIntervalPercentContainer, text="%").pack(side=LEFT, anchor=W)
        Tooltip(self.confidenceIntervalPercentContainer, text="1 sigma: 68.75, 2 sigma: 95.45, 3 sigma: 99.73. 2-way p<0.05 @ = 83.4", wraplength=self.wraplength)
        
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
        
        if self.options.DVHFileType.get() in ["ECLIPSE", "RayStation"]:
            self.customDVHHeaderEntry['state'] = 'disabled'
        
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
        
        # Label(self.middleRightContainer, text="LOG", width=15).pack()
        self.logtext = Text(self.middleRightLowerContainer, height=38, width=75)
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
        self.buttonCalculateNTCP = Button(self.bottomContainer2, text="Calculate NTCP model (N)", command=self.calculateNTCPCommand, width=self.button_width, state=DISABLED)
        self.buttonLKBuncert = Button(self.bottomContainer2, text="Calculate confidence intervals (B)", command=self.calculateLKBuncert, width=self.button_width, state=DISABLED)
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
        
        if not res: self.log("Could not load options.cfg, using default values.")

    # Imported modules   
    from ._Actions import myQuit, log, loadPatientsCommand, addPatientCohort, removePatientCohort
    from ._Actions import chooseStructureCommand, calculateGEUDCommand, toxLimitChange, toxFromFilenameCommand
    from ._Actions import customDVHHeaderCommand, selectDVHFileTypeCommand, showDVHCommand, showDVHPlotCommand
    from ._Actions import bootstrapCorrectionMethodCommand, calculateNTCPCommand, NTCPcalculationCommand
    from ._Actions import calculateAUROCCommand, calculateDVHCommand, cancelCalculateDVHvalues, aggregateDVHCommand
    from ._Actions import showGEUDvsN, switchNto, switchMto, switchTD50to, changeNamingCommand, calculateNewNamesCommand
    from ._Actions import changeNamingQuitCommand, changeNamingQuitAndSaveCommand, drawPlanAndStructureNames
    from ._Actions import customAggregatedDVHCommand, cancelCustomAggregateDVHCommand, matchCustomAggregateDVHCommand
    from ._Actions import saveCustomAggregateDVHCommand, packCustomAggregateDVHCommand

    from ._Analysis import calculateDVHvalues, calculateAggregatedDVH, calculateLKBuncert
