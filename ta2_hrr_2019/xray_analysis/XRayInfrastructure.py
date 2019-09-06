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
import csv

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


def ImportXRayTransmissionsCSV(csv_name):
    LineLength = 0
    with open(csv_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            LineLength += 1
    cnt = 0
    energy = np.empty([LineLength-1, ])
    with open(csv_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            if cnt == 0:
                VariableNames = row
                T = np.empty([LineLength-1, len(VariableNames)-1])
                cnt += 1
                continue
            energy[cnt-1] = float(row[0])
            T[cnt-1, :] = np.array( row[1:] ).astype('float')
            cnt += 1
    return VariableNames, energy, T


def XRayFilters(File):
    file = open(File, "r")
    LineVariableNames = file.readline()
    VariableNames = LineVariableNames.split(",")
    for i in range(0, len(VariableNames)):
        VariableNames[i] = VariableNames[i].replace('\n', '')
    file.close()
    return VariableNames


def loadImageMask(FilterName, SettingPath):
    FilterName.replace('\n', '')
    maskingName = '_Mask.mat'
    FileName = FilterName.replace('\n', '') + maskingName
    PathToLoad = os.path.join(SettingPath, FileName)
    MaskFile = scipy.io.loadmat(PathToLoad)
    Mask = MaskFile['MaterialMask']
    return Mask


def getAverageTransmission(Image, Mask):
    Tmp = Image[Mask > 0]
    MeanValue = sum(Tmp) / len(Tmp)
    return MeanValue


def getRequiredCalibration(SettingPath):
    SettingFile = '2019TA2Transmission.csv'
    CameraSettingFile = 'CameraSettings.txt'
    BckImgPath = 'Darkfield.mat'
    VariableNames, energy, T = ImportXRayTransmissionsCSV(os.path.join(SettingPath, SettingFile))
    BackgroundSignal = 0
    FilterNames = []
    MaskHere = []
    for i in range(0, len(VariableNames)):
        if VariableNames[i] == 'Energy':
            continue
        FilterNames.append(VariableNames[i])
        MaskHere.append(loadImageMask(VariableNames[i], SettingPath))
    File = os.path.join(SettingPath, CameraSettingFile)
    T, FilterNames, MaskHere = sortTransmissions(energy, T, FilterNames, MaskHere)
    PixelSize, GasCell2Camera, RepRate, Alpha, sigma_Alpha = getCameraSettings(File)    
    MaskFile = scipy.io.loadmat(os.path.join(SettingPath, BckImgPath))
    BackgroundImage = MaskFile['BackgroundImage']
    BackgroundImageStd = MaskFile['BackgroundImageStd']
    return energy, T, FilterNames, MaskHere, PixelSize, GasCell2Camera, RepRate, Alpha, sigma_Alpha, BackgroundImage, BackgroundImageStd


def getXRaySignal(image, Masks, FilterNames, BackgroundImage = 0):
    image = image - BackgroundImage
    # one of the filters is the background Tungsten = W
    AverageValues = np.zeros(len(FilterNames) - 1)
    BackgroundSignal = 0
    j = 0
    for i in range(0, len(FilterNames)):
        MaskHere = Masks[i]
        if FilterNames[i] == 'W':
            BackgroundSignal = getAverageTransmission(image, MaskHere)
        else:
            AverageValues[j] = getAverageTransmission(image, MaskHere)
            j += 1
    for j in range(0, len(AverageValues)):
        AverageValues[j] = AverageValues[j] - BackgroundSignal # hard hits
    AverageValues = normaliseArrayOnAverage(AverageValues)
    PeakIntensity = np.amax(AverageValues)
    return AverageValues, PeakIntensity


def sortTransmissions(energy, T, VariableNames, MaskHere):
    SumT = np.empty(T.shape[1],)
    for i in range(0, T.shape[1]):
        SumT[i] = scipy.integrate.trapz(T[:, i], energy)
    Order = np.argsort(SumT)
    NewVariableNames = [VariableNames[Order[0]]]
    NewMaskHere = [MaskHere[Order[0]]]
    NewT = np.zeros(T.shape)
    NewT[:, 0] = T[:, Order[0]]
    for i in range(1, len(Order)):
        NewVariableNames.append(VariableNames[Order[i]])
        NewMaskHere.append(MaskHere[Order[i]])
        NewT[:, i] = T[:, Order[i]]
    return NewT, NewVariableNames, NewMaskHere


def sortTransmissionsValues(T, AverageValues, VariableNames):
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
