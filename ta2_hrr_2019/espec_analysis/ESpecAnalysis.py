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
import skimage.transform
import matplotlib.pyplot as plt
import time
import csv
try:
    from imageTransformation import *
except ImportError:
    from ta2_hrr_2019.espec_analysis.imageTransformation import *
import pkg_resources

import ta2_hrr_2019


def getCalibrationFromCSV(runName, DataPath=ta2_hrr_2019.utils.DATA_FOLDER):
    # csv_file = os.path.join(os.path.realpath(__file__), 'calibration.csv')
    csv_name = pkg_resources.resource_filename(__name__, 'calibration.csv')
    with open(csv_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            RunDateHere, FileNameHere = row
            if RunDateHere <= runName:
                # change it here, depending on how we do it with the calibration information
                CalibrationFile = FileNameHere
                print('Found Calibration')
            else:
                break
    CalibrationFilePath = os.path.join(*[DataPath, "Calibrations", "HighEspec", CalibrationFile])
    J, W, pts, E, dxoverdE, BckgndImage = constantsDefinitions(CalibrationFilePath)
    return J, W, pts, E, dxoverdE, BckgndImage


def TupelOfFiles(path="", Filetype='.tif'):
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


def getESpecSettings(File):
    file = open(File, "r")
    file.readline()
    LineTemporary = file.readline()
    TemporaryValues = LineTemporary.split(",")
    Values = [float(TemporaryValues[i]) for i in range(0, len(TemporaryValues))]
    PixelSize, GasCell2Camera, RepRate, Alpha, sigma_Alpha = Values
    return PixelSize, GasCell2Camera, RepRate, Alpha, sigma_Alpha


def loadMatFile(SettingPath, FileName, VariableName):
    PathToLoad = os.path.join(SettingPath, FileName)
    VabiableLibrary = scipy.io.loadmat(PathToLoad)
    Variable = []
    for i in range(0, len(VariableName)):
        Variable.append(VabiableLibrary[VariableName[i]])
    return Variable


def imageTransformation(Image, LimitView, J2, Projection, ROI1, ROI2):
    Image = Image[LimitView[0]:LimitView[1], LimitView[2]:LimitView[3]]
    Image = skimage.transform.warp(Image, Projection)
    Image = np.multiply(Image, 1 / J2)
    Image = Image.transpose()
    Image = Image[ROI1, ROI2]
    return Image


def uint16ToDoubleJ2(J2, UintConversion):
    J2 = J2.astype(float) / UintConversion[0] * UintConversion[1] + UintConversion[2]
    return J2


def constantsDefinitions(SettingPath):
    """
    Define paths and constanst here. It will later be simply loaded. The
    """
    ImageWarpingFile = 'CalibrationParamters.mat'
    ImageWarpingVariable = ['Jacobian', 'Length', 'Width', 'pts', 'BckgndImage']
    CompleteSetVar = loadMatFile(SettingPath, ImageWarpingFile, ImageWarpingVariable)
    J, L, W, pts, BckgndImage = CompleteSetVar
    ImageWarpingFile = 'PositionVsEnergy.mat'
    ImageWarpingVariable = ['Screen', 'EnergyOnAverage']
    ScreenEnergyOnAverage = loadMatFile(SettingPath, ImageWarpingFile, ImageWarpingVariable)
    Screen, EnergyOnAverage = ScreenEnergyOnAverage
    poly3 = np.poly1d(np.polyfit(Screen[0, :], EnergyOnAverage[0, :], 3))
    E = poly3(L / 1e3)
    dpoly3 = np.polyder(poly3)
    dEoverdx = dpoly3(L / 1e3)
    dxoverdE = 1 / dEoverdx
    # dEoverE = np.multiply(dE, 1/E)
    return J, W, pts, E, dxoverdE, BckgndImage


def calibrationFunction(ImagePath, BackgroundPath, SavePath):
    """
    This can be run to sort out all the calibrations needed. It requires some manual analysis:
    pts: are the points, which are taken at the corners on the edges of the spectrometer screen to flatten.
    Length: mm length of the espec screen, manually analised (= m. a.)
    Width: mm width of the espec screen (m. a.)
    CentrePoint: point where the tape is a "CentrePointDistance" away from the start of the tape (m. a.)
    CentrePointDistance: distance from start of lanex to centre point in mm (m. a.)

    Steps to do this recommended:
    1) Find the pts with matlab or imageJ
    2) Make Sure there is an image with spatial information in the folder of this script. Alternatively add the
        file path in the line FileList = TupelOfFiles(<input path>) as an input argument.
    3) run this script to save a warped image (currently saved as .mat file, but adjust as pleased)
    4) manual find the CentrePoint
    5) Length, Width, CentrePointDistance should be measured from the experimental set-up
    6) A calibration file will be save in "Settings" containing all the important parameters and the path needs
        to be change in the function "constantsDefinitions():"
    +) The path of the background images need to be defined. Otherwise the code just subtracts 0
    :return:
    """
    pts = np.array([[46, 849], [46, 1158], [2035, 871], [2033, 1136]])

    FileList = TupelOfFiles(ImagePath)
    ImageWithSpatialPattern = ImportImageFiles(FileList)

    WarpedImage, M = four_point_transform(ImageWithSpatialPattern[:, :, 0], pts)
    scipy.io.savemat(os.path.join(SavePath, 'WarpedImage.mat'), {'WarpedImage': ImageWithSpatialPattern})
    Length = 230
    Width = 30
    CentrePoint = 991
    CentrePointDistance = 255
    J, L = getJacobianAndSpatialCalibration(WarpedImage, M, Length, CentrePoint, CentrePointDistance)
    PixelWidth = WarpedImage.shape[0]
    W = np.arange(0, PixelWidth) - round(PixelWidth / 2)
    W = W / PixelWidth * Width
    # define the path of the background images here:
    BckgndImage = backgroundImages(BackgroundPath, J, pts)
    scipy.io.savemat(os.path.join(SavePath, 'CalibrationParamters.mat'), {'Jacobian': J, 'Length': L, 'Width': W,
                                                                          'pts': pts, 'BckgndImage': BckgndImage})
    return J, L, W, BckgndImage


def backgroundImages(Path, J, pts):
    if len(Path) > 0:
        FileList = TupelOfFiles(Path)
        BackgroundImages = ImportImageFiles(FileList)
        BackgroundImage = np.mean(BackgroundImages, 2)
        WarpedBackgroundImage = four_point_transform(BackgroundImage, pts, J)
    else:
        WarpedBackgroundImage = 0
    return WarpedBackgroundImage


def electronSpectrum(Image, dxoverdE):
    SumCounts = np.sum(Image, 0)
    Spectrum = np.multiply(SumCounts, dxoverdE)
    return Spectrum


"""
def analysisImages(PathOfImage, Filetype='.tif'):
    if PathOfImage[-4:] == 'tiff' or PathOfImage[-4:] == '.tif':
        FileList = PathOfImage
    else:
        FileList = TupelOfFiles(PathOfImage, Filetype)
    Images = ImportImageFiles(FileList)
    Image = np.mean(Images, 2)
    J, W, pts, E, dxoverdE, BckgndImage = constantsDefinitions()
    WarpedImage, __ = four_point_transform(Image, pts, J)
    PureWarpedImage = WarpedImage - BckgndImage
    Spectrum = electronSpectrum(PureWarpedImage, dxoverdE)
    return E, Spectrum
"""


class HighEspec:
    """
    Function to load calibration values for the electron spectrometer.
    function: __init__ takes a run name and the path of the data.
    returns after initialising:
    return: BckgndImage  Background image (isnt done coding yet)
    return:  dxoverdE change the spectrum from a spatial domain into an energy domain
    return: E energy domain
    return: pts point for warping the image
    return: W width of the camera image in [m] (at the moment)
    """
    def __init__(self, runName, DataPath=ta2_hrr_2019.utils.DATA_FOLDER):
        J, W, pts, E, dxoverdE, BckgndImage = getCalibrationFromCSV(runName, DataPath)
        self.J = J
        self.BckgndImage = BckgndImage
        self.dxoverdE = dxoverdE
        self.E = E
        self.pts = pts
        self.W = W

    def SpectrumFromImage(self, rawImageuint16):
        image = rawImageuint16.astype(float)
        WarpedImage, __ = four_point_transform(image, self.pts, self.J)
        WarpedImageWithoutBckgnd = WarpedImage - self.BckgndImage
        Spectrum = electronSpectrum(WarpedImages[i], self.dxoverdE)
        return Spectrum

    def AverageSpectrum(self, Spectra):
        if not isinstance(Spectra, list):
            Spectra = [Spectra]
        Spectrum = np.zeros(Spectra[0].shape)
        for i in range(0, len(Spectra)):
            Spectrum = Spectrum + Spectra[i]
        Spectrum = Spectrum / len(Spectra)
        return Spectrum


if __name__ == "__main__":
    """
    Things done for initialising the espec calibration files and analysing the performance:
    calibrationFunction() -> saving calibration files, see descriptions
    Checking performance about 
    b = time.time()
    for i in range(0, 1000):
        J, W, pts, E = constantsDefinitions()
    a = time.time()
    print(a - b)
        
    """
    # once the background pictures are acquired, read the description of calibrationFunction() and run it once
    # calibrationFunction()

    # to analyse espec images run:
    # 1) define PathOfImage where the images of a run are
    # 1.1) If PathOfImage is a folder it will average all images and get the spectrum
    # 1.2) If PathOfImage is a single tif or tiff file it will give the spectrum for just that.
    # 2) run: EnergyVector, Spectrum = analysisImages(PathOfImage)

    """
    old stuff:
    b = time.time()
    a = time.time()
    print(a - b)
    # me analysing some part of the code:
    b = time.time()
    FileList = TupelOfFiles()
    NumberOfImages = 50
    ImageWithSpatialPattern = ImportImageFiles(FileList)
    J, W, pts, E, dxoverdE, BckgndImage = constantsDefinitions()
    Images = np.zeros(shape=(ImageWithSpatialPattern.shape[0], ImageWithSpatialPattern.shape[1], NumberOfImages))
    Images[:, :, 0] = ImageWithSpatialPattern[:, :, 0]
    for i in range(1, NumberOfImages):
        Images[:, :, i] = ImportImageFiles(FileList)[:, :, 0]
    a = time.time()
    d = time.time()
    Image = np.mean(Images, 2)
    WarpedImage, __ = four_point_transform(Image, pts, J)
    PureWarpedImage = WarpedImage - BckgndImage
    Spectrum = electronSpectrum(PureWarpedImage, dxoverdE)
    c = time.time()
    print(a - b)
    print(c - d)
    """
