import json
import numpy as np
from os.path import isfile
from tools import *
from dataset import Dataset
from config import Config
from recording import Recording
from analyzer import Analyzer
import pickle


# Change if you want your screen to keep quiet
# 0 = quiet
# 1 = in between
# 2 = detailed information
verbatim = 2

# Init project with configuration file
config = Config('../config/general/newsettings_10.json', verbatim=verbatim)
config.readAndCheck()

# Make or load analyzer (model+scaler)
analyzer = Analyzer(config, verbatim=verbatim)
# allData, allLabels = analyzer.learn(config, returnData=True) # If you want the data
analyzer.learn(config)
analyzedSet = Dataset(config,verbatim=verbatim)
analyzedSet.analyze(analyzer, config, save=True)
analyzedSet.makeDecision(config, save=True)
analyzedSet.display(config, onlineDisplay=False, saveDisplay=True, forChecking=False)
#analyzedSet.display(config, onlineDisplay=False, saveDisplay=True, forChecking=True, labelEncoder=analyzer.labelEncoder)
#trueLabels,predictedLabels=analyzedSet.getNumericalResults(config, analyzer.labelEncoder)
