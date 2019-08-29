#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 14:45:56 2019
XRayFilterPack
@author: Gruse
"""
import os
import scipy.io
import numpy as np


def TupelOfFiles(path="", Filetype=".tif"):
    if len(path) == 0:
        path = "."
    FileList = []
    for files in os.listdir(path):
        if files.endswith(Filetype):
            FileList.append(os.path.join(path, files))
            print(FileList[-1])
    return FileList


def ImportImageFiles(FileList):
    from PIL import Image
    im = Image.open(FileList[0])
    imarray = np.array(im)
    Images = np.zeros(shape=(imarray.shape[0], imarray.shape[1], len(FileList)))
    for i in range(0, len(FileList)):
        im = Image.open(FileList[i])
        Images[:, :, i] = np.array(im)
    return Images


def getCameraSettings(File):
    file = open(File, "r")
    file.readline()
    LineTemporary = file.readline()
    TemporaryValues = LineTemporary.split(",")
    Values = [float(TemporaryValues[i]) for i in range(0, len(TemporaryValues))]
    PixelSize, GasCell2Camera, RepRate, Alpha, sigma_Alpha = Values
    return PixelSize, GasCell2Camera, RepRate, Alpha, sigma_Alpha


def ImportXRayTransmissions(File):
    VariableNames = XRayFilters(File)
    file = open(File, "r")
    # First Line only consists of the variable names and needs to be neglected:
    file.readline()
    # Then the data starts.
    # This first part of the code just checks how many data points exist in the file
    LineTemporary = file.readline()
    cnt = 0
    while LineTemporary:
        cnt += 1
        LineTemporary = file.readline()
    energy = np.zeros(cnt)
    # This has to be adapted if there is no background filter (in this case its Tungsten)
    # '-1' instead of '-2' if there is none
    T = np.zeros(shape=(cnt, len(VariableNames) - 2))
    file.close()
    file = open(File, "r")
    # as before ignoring the line with the variable names:
    file.readline()
    LineTemporary = file.readline()
    cnt = 0
    while LineTemporary:
        TemporaryValues = LineTemporary.split(",")
        TemporaryValues = [float(TemporaryValues[i]) for i in range(0, len(TemporaryValues))]
        energy[cnt] = float(TemporaryValues[0])
        # temporaryTransmissionValues = TemporaryValues[([1, 3, 4, 5, 6, 7])]
        temporaryTransmissionValues = [float(TemporaryValues[i]) for i in [1, 3, 4, 5, 6, 7]]
        T[cnt, :] = temporaryTransmissionValues
        cnt += 1
        LineTemporary = file.readline()
    return VariableNames, energy, T


def XRayFilters(File):
    file = open(File, "r")
    LineVariableNames = file.readline()
    VariableNames = LineVariableNames.split(",")
    for i in range(0, len(VariableNames)):
        VariableNames[i] = VariableNames[i].replace('\n', '')
    file.close()
    return VariableNames


def loadImageMask(FilterName, SettingPath, MaskingName='_Mask.mat'):
    FilterName.replace('\n', '')
    FileName = FilterName.replace('\n', '') + MaskingName
    PathToLoad = os.path.join(SettingPath, FileName)
    MaskFile = scipy.io.loadmat(PathToLoad)
    Mask = MaskFile['MaterialMask']
    return Mask


def getAverageTransmission(Image, Mask):
    Tmp = Image[Mask > 0]
    MeanValue = sum(Tmp) / len(Tmp)
    return MeanValue


def getXRaySignal(XRayImage, SettingPath, SettingFile, MaskingName):
    VariableNames, energy, T = ImportXRayTransmissions(os.path.join(SettingPath, SettingFile))
    AverageValues = np.zeros(len(VariableNames) - 2)
    BackgroundSignal = 0
    j = 0
    FilterNames = []
    for i in range(0, len(VariableNames)):
        if VariableNames[i] == 'Energy':
            continue
        MaskHere = loadImageMask(VariableNames[i], SettingPath, MaskingName)
        if VariableNames[i] == 'W':
            BackgroundSignal = getAverageTransmission(XRayImage, MaskHere)
        else:
            AverageValues[j] = getAverageTransmission(XRayImage, MaskHere)
            FilterNames.append(VariableNames[i])
            j += 1
    for j in range(0, len(AverageValues)):
        # StrToPrint1 = "%s Signal %f" % (FilterNames[j], AverageValues[j])
        # print(StrToPrint1)
        AverageValues[j] = AverageValues[j] - BackgroundSignal
        # StrToPrint2 = "%s Signal without Bckgrnd %f" % (FilterNames[j], AverageValues[j])
        # print(StrToPrint2)
    PeakIntensity = AverageValues[0]
    AverageValues = normaliseArrayOnAverage(AverageValues)
    T, AverageValues, FilterNames = sortTransmissions(T, AverageValues, FilterNames)
    return energy, T, AverageValues, FilterNames, PeakIntensity


def sortTransmissions(T, AverageValues, VariableNames):
    SumT = np.sum(T, axis=0)
    Order = np.argsort(SumT)
    NewAverageValues = AverageValues[Order]
    NewVariableNames = [VariableNames[Order[0]]]
    NewT = np.zeros(T.shape)
    NewT[:, 0] = T[:, Order[0]]
    for i in range(1, len(Order)):
        NewVariableNames.append(VariableNames[Order[i]])
        NewT[:, i] = T[:, Order[i]]
    return NewT, NewAverageValues, NewVariableNames


def normaliseArrayOnAverage(InArray):
    MeanValue = np.sum(InArray)
    for i in range(0, len(InArray)):
        InArray[i] = InArray[i] / MeanValue
    return InArray
