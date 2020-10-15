from _Options import *
import numpy as np


class ConfidenceInterval:
    def __init__(self):
        self.histogram = list()

    def getMean(self, histogram):
        return np.mean(histogram)

    def getMedian(self, histogram):
        return np.median(histogram)

    def getUpperCI(self, histogram):
        pass

    def getLowerCI(self, histogram):
        pass

    def doCorrectPivot(self, histogram):
        pass


class LKB(ConfidenceInterval):
    def __init__(self, n, m, TD50, options):
        self.options = options
        self.n = n
        self.m = m
        self.TD50 = TD50
        self.nIterations = self.options.ConfidenceIntervalIterations.get()
        self.histogram = {'n': np.zeros(nIterations), 'm': np.zeros(nIterations), 'TD50': np.zeros(nIterations)}

    def __repr__(self):
        return {'n': self.n, 'm': self.m, 'TD50': self.TD50}


class Logit(ConfidenceInterval):
    def __init__(self, a, b, options):
        self.options = options
        self.a = a
        self.b = b
        self.nIterations = self.options.ConfidenceIntervalIterations.get()
        self.histogram = {'a': np.zeros(nIterations), 'b': np.zeros(nIterations)}
        self.idx = 0

    def __repr__(self):
        return {'a': self.a, 'b': self.b}


class Results:
    def __init__(self, options):
        self.LKB = LKB(None, None, None)
        self.Logit = Logit(None, None, None)
        self.options = options
        self.currentModel = self.options.NTCPcalculation.get() == "LKB" and self.LKB or self.Logit
        self.parameterList = self.options.NTCPcalculation.get() == "LKB" and ['n', 'm', 'TD50'] or ['a', 'b']

    def __repr__(self):
        return f"Results({self.currentModel})"

    def getCurrentModel(self):
        return self.options.NTCPcalculation.get()

    def updateOptions(self, options):
        self.options = options

    def setParameters(self, parameters):
        if isinstance(parameters, dict):
            for key, value in parameters.items():
                assert hasattr(self.currentModel, key)
                setattr(self.currentModel, key, value)

        elif isinstance(parameters, list):
            for key, value in zip(self.parameterList, parameters):
                assert hasattr(self.currentModel, key)
                setattr(self.currentModel, key, value)

    def addHistogram(self, parameters):
        if isinstance(parameters, dict):
            for key, value in parameters.items():
                self.currentModel[key][self.idx] = value

        elif isinstance(parameters, list):
            for key, value in zip(self.parameterList, parameters):
                self.currentModel[key][self.idx] = value

        self.idx += 1

    def getOutputParameters(self):
        return f"Results({self.currentModel})"

    def getOutputParametersCI(self):
        output = str()
        for key, value in self.currentModel:
            output += f"{key} = {value.getMean():.2f} ({value.getLowerCI():.2f}-{value.getUpperCI():.2f})\n"
        return output
