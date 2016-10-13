__author__ = "Ben Nelsen"

import backwaterLength


s_folderPath = "C:\\Users\\bnelsen2\\Documents\\SeniorYearOne\\CEEN433"
s_fileName = "backwaterConfigFile2.txt"
s_writeFileName = "backwaterProblem2.txt"

backwater = backwaterLength.Channel()
backwater.readConfigurationfile(s_folderPath, s_fileName)
data = backwater.backwaterDetermination()
backwater.outputPrinting(data, s_folderPath, s_writeFileName)
