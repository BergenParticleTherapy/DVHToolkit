import os
from tkinter import *


class Options():
    def __init__(self):
        self.PROGRAM_VERSION = 1.64
        self.Gy = 1
        self.dGy = 0.1
        self.cGy = 0.01
        self.mGy = 0.001
        self.cc = 1

        self.DVHFileType = StringVar(value="ECLIPSE")  # ['simple', 'ECLIPSE', 'RayStation']
        self.NTCPcalculation = StringVar(value="Logit")  # ['LKB', 'Logit']
        self.NTCPcalculationDpercent = DoubleVar(value=20)  # [0 -> 100]
        self.useNTCPcc = IntVar(value=0)
        self.dataFolder = StringVar(value="")
        self.structureToUse = StringVar(value="")
        self.planToUse = StringVar(value="")
        self.customDVHHeader = StringVar(value="Dose,Volume")
        self.CSVStyle = StringVar(value="autodetect")  # ['periodComma', 'commaSemicolon', 'autodetect']
        self.skipRows = IntVar(value=0)
        self.doseUnit = StringVar(value="autodetect")  # incl. autodetect
        self.optimizationScheme = StringVar(value="GradientDescent")  # ['GradientDescent', 'MatrixMinimization']
        self.optimizationMetric = StringVar(value="LLH")  # ['LLH', 'LS']
        self.confidenceIntervalMethod = StringVar(value="NonParametricBootstrapping")  # ['ProfileLikelihood', 'ParametricBootstrapping', 'NonParametricBootstrapping']
        self.confidenceIntervalIterations = IntVar(value=1500)
        self.confidenceIntervalPercent = DoubleVar(value=95)
        self.confidenceIntervalLikelihoodLimit = DoubleVar(value=-1.0)
        self.makeDifferentialCI = IntVar(value=0)
        self.bootstrapCorrectionMethod = StringVar(value="median")  # ["none", "mean", "median"]
        self.confidenceIntervalShowModels = IntVar(value=0)
        self.toxLimit = IntVar(value=2)
        self.matrixSize = IntVar(value=50)
        self.loadToxFromFilename = IntVar(value=1)
        self.autodetectDVHHeader = IntVar(value=1)

        self.fixN = IntVar(value=0)
        self.nFrom = DoubleVar(value=0.02)
        self.nTo = DoubleVar(value=1)
        self.nSet = DoubleVar(value=1)

        self.fixM = IntVar(value=0)
        self.mFrom = DoubleVar(value=0.03)
        self.mTo = DoubleVar(value=1)

        self.fixTD50 = IntVar(value=0)
        self.TD50From = DoubleVar(value=10)
        self.TD50To = DoubleVar(value=85)

        self.fixA = IntVar(value=0)
        self.aFrom = DoubleVar(value=-100)
        self.aTo = DoubleVar(value=0)

        self.fixB = IntVar(value=0)
        self.bFrom = DoubleVar(value=0)
        self.bTo = DoubleVar(value=2)

        self.basinHoppingIterations = IntVar(value=10)
        self.basinHoppingTemperature = DoubleVar(value=0.1)
        self.basinHoppingNsize = DoubleVar(value=0.3)
        self.basinHoppingMsize = DoubleVar(value=0.3)
        self.basinHoppingTD50size = DoubleVar(value=10)
        self.basinHoppingAsize = DoubleVar(value=40)
        self.basinHoppingBsize = DoubleVar(value=0.5)

        self.NTCPTimeDependent = IntVar(value=0)
        self.fixLambda = IntVar(value=1)
        self.lambdaFrom = DoubleVar(value=0.38)
        self.lambdaTo = DoubleVar(value=0.54)
        self.NTCPLambdaSize = DoubleVar(value=0.38)
        self.fixGamma = IntVar(value=1)
        self.gammaFrom = DoubleVar(value=1.37)
        self.gammaTo = DoubleVar(value=2)
        self.NTCPGammaSize = DoubleVar(value=1.37)

        self.dvhPlotUseToxAsColor = IntVar(value=0)
        self.dvhPlotLegendMarker = StringVar(value="structureName")  # ["structureName", "folderName", "planName"]
        self.dvhPlotLegendSize = DoubleVar(value=14)
        self.dvhPlotLineStyleGrouping = StringVar(value="plan")
        self.dvhPlotSeparatePlots = StringVar(value="structure")
        self.dvhPlotsSaveFigs = IntVar(value=1)

        self.vars = {"DVHFileTypeS": self.DVHFileType, "NTCPcalculationS": self.NTCPcalculation,
                     "NTCPcalculationDpercentI": self.NTCPcalculationDpercent, "dataFolderS": self.dataFolder,
                     "structureToUseS": self.structureToUse, "CSVStyleS": self.CSVStyle,
                     "optimizationSchemeS": self.optimizationScheme, "optimizationMetricS": self.optimizationMetric,
                     "toxLimitI": self.toxLimit, "loadToxFromFilenameI": self.loadToxFromFilename, "matrixSizeI": self.matrixSize,
                     "fixNI": self.fixN, "nFromD": self.nFrom, "nToD": self.nTo,
                     "fixMI": self.fixM, "mFromD": self.mFrom, "mToD": self.mTo,
                     "fixTD50I": self.fixTD50, "TD50FromD": self.TD50From, "Td50ToD": self.TD50To,
                     "fixA": self.fixA, "aFrom": self.aFrom, "aTo": self.aTo,
                     "fixB": self.fixB, "bFrom": self.bFrom, "bTo": self.bTo,
                     "basinHoppingIterationsI": self.basinHoppingIterations, "basinHoppingTemperatureD": self.basinHoppingTemperature,
                     "basinHoppingNsizeD": self.basinHoppingNsize, "basinHoppingMsizeD": self.basinHoppingMsize,
                     "basinHoppingTD50sizeD": self.basinHoppingTD50size, "confidenceIntervalMethod": self.confidenceIntervalMethod,
                     "confidenceIntervalIterations": self.confidenceIntervalIterations, "confidenceIntervalPercent": self.confidenceIntervalPercent,
                     "customDVHHeader": self.customDVHHeader, "skipRows": self.skipRows, "doseUnit": self.doseUnit,
                     "confidenceIntervalShowModels": self.confidenceIntervalShowModels,
                     "bootstrapCorrectionMethod": self.bootstrapCorrectionMethod, "basinHoppingAsizeD": self.basinHoppingAsize,
                     "basinHoppingBsizeD": self.basinHoppingBsize, "dvhPlotUseToxAsColor": self.dvhPlotUseToxAsColor,
                     "dvhPlotLegendMarker": self.dvhPlotLegendMarker, "dvhPlotLineStyleGrouping": self.dvhPlotLineStyleGrouping,
                     "dvhPlotSeparatePlots": self.dvhPlotSeparatePlots, "dvhPlotsSaveFigs": self.dvhPlotsSaveFigs,
                     "dvhPlotLegendSize": self.dvhPlotLegendSize, "makeDifferentialCI": self.makeDifferentialCI,
                     "useNTCPcc": self.useNTCPcc, "nSet": self.nSet}

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
        with open("config.cfg", "w") as configFile:
            for key, var in list(self.vars.items()):
                configFile.write(f"{key},{var.get()}\n")

    def getCustomHeaderIdx(self):
        try:
            header = self.customDVHHeader.get().split(",")
        except:
            header = ["Dose", "Volume"]

        doseIdx = header.index("Dose")
        volumeIdx = header.index("Volume")
        return "{}{}".format(doseIdx, volumeIdx)
