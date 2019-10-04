#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Wed Jun 19 09:49:52 2019

@author: chrisunderwood
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [6.0,4.0]
import matplotlib.pyplot as plt

# Load my module of functions
import CUnderwood_Functions3 as func

savePath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/"

rootExperimentFolder = "/Volumes/GoogleDrive/Shared drives/Murphy Group/GeminiRR_2019/"
fileTag = "_Probe_Interfero.tiff"

# The format of the saving at Gemini is yyyymmdd, search for begining of the year
DayFolders = func.SubFolders_in_Folder(rootExperimentFolder)
DayFolders = [x for x in DayFolders if x.startswith("20")]


outputData = {}
for day in DayFolders:
    dayPath = rootExperimentFolder + day + "/"
    runFolders = func.SubFolders_in_Folder(dayPath)
    runFolders = [x for x in runFolders if x.startswith("20")] # Search for begining of the year again

    print (day, runFolders)

    for run in runFolders:
        runPath = dayPath + run + "/"
        fileList = func.FilesInFolder(runPath, fileTag)
        if len(fileList) > 0:
# =============================================================================
# Sort the list into order by shot number
# This can be done in two ways, so check that the method chosen works for the case
# =============================================================================
            # sortIndexes = func.sortArrSplice(fileList, 13,16)
            sortIndexes = func.SplitArr(func.SplitArr(fileList, "s", 1), "_", 0)

            # Sort the files by their shot number
            fileList, sortIndexes = func.sortArrAbyB(fileList, sortIndexes)
            filePathList = []
            for f in fileList:
                filePathList.append(runPath + f)
            outputData[run] = filePathList

runFolders = list(outputData)
Date = 0.001 * np.array(func.SplitArr(runFolders, "r", 1)) +  np.array(func.SplitArr(runFolders, "r", 0))
runFolders, Date = func.sortArrAbyB(runFolders, Date)

# =============================================================================
# Write the data to a csv file
# =============================================================================
outFile = open(savePath + "2019RadiationReactionDatabase.csv", "w")
outFile.write("Day,Run,Shot,FilePath,Load (1) / Make (0),Use same as x, LoadPathReference, Fringe, , , ,Angle,CropROI, , , ,PC_height, PC_heightSize,crop_To_PlasmaChannel, , , ,FFTCrop, , , ,Laser Wavelenth, m per pxiels\n")

laserwavelength = 8e-07
sizePerPixel = 4.22e-05
for run in runFolders:
    print (run)
    day = run.split("r")[0]
    runNos =  run.split("r")[1]
    for filePath in outputData[run]:
        shotNumber = filePath.split("/")[-1].split("s")[1].split("_")[0]
        outFile.write("{},{},{},{},0, , ,,,,,,,,,,,,,,,,,,,,{},{}\n".format(day, runNos, shotNumber, filePath, laserwavelength, sizePerPixel))
        
outFile.close()


