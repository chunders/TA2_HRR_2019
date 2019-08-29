#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 14:45:56 2019
XRayFilterPack
@author: Gruse
"""
from typing import Tuple, Union

import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from numpy.core.multiarray import ndarray
from scipy.special import kv
from scipy.integrate import trapz
from scipy import optimize
import time
from XRayInfrastructure import *
try:
    from XRayInfrastructure import *
except ValueError:
    from ta2_hrr_2019.xray_analysis.XRayInfrastructure import *
import pkg_resources


def synchrotronFunction(energy, ecrit):
    phi = energy / (2 * ecrit)
    K2over3 = scipy.special.kv(2 / 3, phi)
    NumberSpectrum = energy / ecrit ** 2 * K2over3 ** 2
    Intensity = NumberSpectrum * energy
    # , NumberSpectrum
    return Intensity


def numberSpectrum(energy, ecrit):
    phi = energy / (2 * ecrit)
    K2over3 = scipy.special.kv(2 / 3, phi)
    NumberSpectrum = energy / ecrit ** 2 * K2over3 ** 2
    return NumberSpectrum


"""
def getPhotonFlux(SettingPath, CameraSettingFile, energy, ecrit, PeakIntensity, TQ):
    File = os.path.join(SettingPath, CameraSettingFile)
    PixelSize, GasCell2Camera, RepRate, Alpha, sigma_Alpha = getCameraSettings(File)
    Theta2 = (PixelSize/GasCell2Camera)**2*1e6
    S = numberSpectrum(energy, ecrit)
    NormFactor = scipy.integrate.trapz(S, energy)
    Snorm = S/NormFactor
    Divisor = np.multiply(Snorm, energy)
    Divisor = np.multiply(Divisor, TQ)
    Divisor = scipy.integrate.trapz(Divisor, energy)
    NPhotons = PeakIntensity/(Alpha*Divisor)
    energy01percent = np.arange(0.9995*ecrit, 1.0005*ecrit, (1.0005*ecrit-0.9995*ecrit)/100)
    S01percent = numberSpectrum(energy01percent, ecrit)/NormFactor
    NPhotons01percent = NPhotons*scipy.integrate.trapz(S01percent, energy01percent)
    NPhotons_01percent_omega_s = NPhotons01percent*RepRate/Theta2
    return NPhotons, NPhotons01percent, NPhotons_01percent_omega_s
"""


def getPhotonFlux(energy, ecrit, PeakIntensity, TQ, PixelSize, GasCell2Camera, RepRate, Alpha):
    Theta2 = (PixelSize/GasCell2Camera)**2*1e6
    S = numberSpectrum(energy, ecrit)
    NormFactor = scipy.integrate.trapz(S, energy)
    Snorm = S/NormFactor
    Divisor = np.multiply(Snorm, energy)
    Divisor = np.multiply(Divisor, TQ)
    Divisor = scipy.integrate.trapz(Divisor, energy)
    NPhotons = PeakIntensity/(Alpha*Divisor)
    energy01percent = np.arange(0.9995*ecrit, 1.0005*ecrit, (1.0005*ecrit-0.9995*ecrit)/100)
    S01percent = numberSpectrum(energy01percent, ecrit)/NormFactor
    NPhotons01percent = NPhotons*scipy.integrate.trapz(S01percent, energy01percent)
    NPhotons_01percent_omega_s = NPhotons01percent*RepRate/Theta2
    return NPhotons, NPhotons01percent, NPhotons_01percent_omega_s


def transmissionSpectrum(ecrit, energy, T):
    TSQE = np.zeros(T.shape)
    Spectrum = synchrotronFunction(energy, ecrit)
    SpectrumMatrix = np.transpose(np.tile(Spectrum, (T.shape[1], 1)))
    TSQE = np.multiply(SpectrumMatrix, T)
    return TSQE


def showTransmissionCurves(energy, TSQE, FilterNames, ecrit):
    for i in range(0, len(FilterNames)):
        plt.semilogx(energy, TSQE[:, i])
        plt.xlabel('Photon Energy [keV]')
        plt.ylabel('Transmission')
    plt.legend(FilterNames)
    StrToPrint2 = "Filter Transmissions for a critical energy of %.0f keV" % ecrit
    plt.title(StrToPrint2)
    plt.show()


def showSignalPoints(energy, T, ecrits, FilterNames, AverageValues):
    plt.plot(range(0, len(FilterNames)), AverageValues)
    plt.xlabel('Filter')
    plt.ylabel('Norm. Transmission')
    LegendList = ["Image"]
    if hasattr(ecrits, "__len__"):
        for i in range(0, len(ecrits)):
            ecrit = ecrits[i]
            TSQE = transmissionSpectrum(ecrit, energy, T)
            TheoryValues = theoreticalTransmissionSignal(energy, TSQE)
            plt.plot(range(0, len(FilterNames)), TheoryValues)
            StrToPrint2 = "%.0f keV" % ecrit
            LegendList.append(StrToPrint2)
    else:
        ecrit = ecrits
        TSQE = transmissionSpectrum(ecrit, energy, T)
        TheoryValues = theoreticalTransmissionSignal(energy, TSQE)
        plt.plot(range(0, len(FilterNames)), TheoryValues)
        StrToPrint2 = "%.1f keV" % ecrit
        LegendList.append(StrToPrint2)
    plt.legend(LegendList)
    plt.xticks(range(0, len(FilterNames)), FilterNames, rotation=12)
    plt.show()


def dataVStheoryToMinimise(data, theory):
    Value = 0
    for i in range(0, len(data)):
        Value += (data[i]-theory[i])**2
    return Value


def fitDataToSynchrotronSpectrum(energy, AverageValues, T):
    TempFunc = lambda ecrit: dataVStheoryToMinimise(theoreticalTransmissionSignal(energy, transmissionSpectrum(ecrit, energy, T)), AverageValues)
    ecritMin = 0.1
    ecritMax = 30
    BestCriticalEnergy = scipy.optimize.fminbound(TempFunc, ecritMin, ecritMax)
    return BestCriticalEnergy


def theoreticalTransmissionSignal(energy, TSQE):
    theoryValue = np.zeros(TSQE.shape[1])
    for i in range(0, TSQE.shape[1]):
        theoryValue[i] = scipy.integrate.trapz(TSQE[:, i], energy)
    theoryValue = normaliseArrayOnAverage(theoryValue)
    return theoryValue


def getCalibrationFromCSV(runName, DataPath=ta2_hrr_2019.utils.DATA_FOLDER):
    # csv_file = os.path.join(os.path.realpath(__file__), 'calibration.csv')
    csv_file = pkg_resources.get_resource_filename(__name__, 'calibration.csv')
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        row[0].split("/")
        if row[0] <= runName:
            # change it here, depending on how we do it with the calibration information
            CalibrationFile = row[1]
        else:
            break
    SettingPath = os.path.join(*[DataPath, "Calibrations", "XRay", CalibrationFile])
    return SettingPath


class XRay:
    def __init__(self, runName, DataPath=ta2_hrr_2019.utils.DATA_FOLDER):
        self.SettingPath = getCalibrationFromCSV(runName, DataPath)
        energy, T, FilterNames, MaskHere, PixelSize, GasCell2Camera, RepRate, Alpha, sigma_Alpha = \
            getRequiredCalibration(self.SettingPath)
        self.energy = energy
        self.T = T
        self.FilterNames = FilterNames
        self.MaskHere = MaskHere
        self.PixelSize = PixelSize
        self.GasCell2Camera = GasCell2Camera
        self.RepRate = RepRate
        self.Alpha = Alpha
        self.sigma_Alpha = sigma_Alpha

    def prepImages(self, image):
        image = image.astype(float)
        AverageValues, PeakIntensity = getXRaySignal(image, self.MaskHere, self.FilterNames)
        return AverageValues, PeakIntensity

    def analyseImages(self, AverageValuesList, PeakIntensityList):
        if not isinstance(AverageValuesList, list):
            AverageValuesList = [AverageValuesList]
            PeakIntensityList = [PeakIntensityList]
        AverageValues = np.zeros(AverageValuesList[0].shape)
        PeakIntensity = 0
        # outter loop through all the filters
        for j in range(0, len(AverageValuesList[0])):
            # inner loop going through all the images
            for i in range(0, len(AverageValuesList)):
                AverageValues[j] = AverageValues[j] + AverageValuesList[i][j]
                if j == 0:
                    PeakIntensity = PeakIntensity + PeakIntensityList[i]
            if j == 0:
                PeakIntensity = PeakIntensity / len(PeakIntensityList)
            AverageValues[j] = AverageValues[j] / len(AverageValuesList)
        ecrit = fitDataToSynchrotronSpectrum(self.energy, AverageValues, self.T)
        NPhotons, NPhotons01percent, NPhotons_01percent_omega_s = getPhotonFlux(self.energy, ecrit, PeakIntensity,
                                                                                self.T[:, -1], self.PixelSize,
                                                                                self.GasCell2Camera, self.RepRate,
                                                                                self.Alpha)
        return ecrit, NPhotons, NPhotons01percent, NPhotons_01percent_omega_s


def xrayAnalysis(Images, SettingPath):
    """
    This function is to run for the anaylis. Coded in
    :param SettingPath: path where the settings can be found
    :param Images: numpy array of images
    :return: ecrit, NPhotons, NPhotons01percent, NPhotons_01percent_omega_s
    For test purposes the following code was partially used (re-use for error analysis):
    SettingPath = 'Settings'
    # ecrit = 10
    SettingFile = '2019TA2Transmission.txt'
    CameraSettingFile = 'CameraSettings.txt'
    ImageFile = '2019TA2TestImage.mat'
    TestImageFile = scipy.io.loadmat(os.path.join(SettingPath, ImageFile))
    TestImage = TestImageFile['img']
    < run the scripts with these settings, then look at the results: >
    showSignalPoints(energy, T, ecrit, FilterNames, AverageValues)
    """
    CameraSettingFile = 'CameraSettings.txt'

    energy, T, AverageValues, FilterNames, PeakIntensity = getXRaySignal(Images, SettingPath)
    ecrit = fitDataToSynchrotronSpectrum(energy, AverageValues, T)
    NPhotons, NPhotons01percent, NPhotons_01percent_omega_s = getPhotonFlux(SettingPath, CameraSettingFile, energy,
                                                                            ecrit, PeakIntensity, T[:, -1])
    return ecrit, NPhotons, NPhotons01percent, NPhotons_01percent_omega_s


if __name__ == "__main__":
    """
    Here are the main parameter, which has to be defined to run the scripts. Also it can be chosen which value is
        returned by the script. For additional information on the parameters look at the description of "xrayAnalysis".
    A script for the 5 Hz run or any other things should be easy to code with this information. 
    In case of question: 
    jg1916@ic.ac.uk
    """
    # Define a path
    ImagePath = ''
    # Whatever filetype
    Filetype = '.tiff'
    # '2019AugustTA2.txt' is the set-up with 700 nm aluminium as light tightning
    SettingFile = '2019AugustTA2.txt'
    MaskingName = '_Mask.mat'

    ecrit, NPhotons, NPhotons01percent, NPhotons_01percent_omega_s = xrayAnalysis(
        ImagePath, Filetype, SettingFile, MaskingName)
