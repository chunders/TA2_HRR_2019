#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 14:45:56 2019
XRayFilterPack
@author: Gruse
"""
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import kv
from scipy.integrate import trapz
from scipy import optimize
import time
from XRayInfrastructure import *


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


def xrayAnalysis(ImagePath, Filetype, SettingFile, MaskingName):
    """
    This function is to run for the anaylis. Coded in
    :param MaskingName: MaskingName is the end of the masking files. These have the names of the filters, plus
        some additional strings to make them unique. In the TA2 2019 run this is either
        a) '_Mask.mat' for the full 1024 x 1024 X-Ray images
        b) '_MaskBinned.mat' for the faster binned 256 x 256 X-Ray images, which are used in the 5 Hz run
    :param Filetype: '.tiff' or '.tif'
    :param SettingFile: '2019JulyTA2.txt' or '2019AugustTA2.txt'
    :param ImagePath: Can be either an individual image or a directory with images
    :return: ecrit, NPhotons, NPhotons01percent, NPhotons_01percent_omega_s

    For test purposes the following code was partially used (re-use for error analysis):
    SettingPath = 'Settings'
    # ecrit = 10
    SettingFile = '2019AugustTA2.txt'
    CameraSettingFile = 'CameraSettings.txt'
    ImageFile = '2019TA2TestImage.mat'
    TestImageFile = scipy.io.loadmat(os.path.join(SettingPath, ImageFile))
    TestImage = TestImageFile['img']
    < run the scripts with these settings, then look at the results: >
    showSignalPoints(energy, T, ecrit, FilterNames, AverageValues)
    """
    SettingPath = 'Settings'
    CameraSettingFile = 'CameraSettings.txt'

    if ImagePath[-4:] == 'tiff' or ImagePath[-4:] == '.tif':
        FileList = ImagePath
    else:
        FileList = TupelOfFiles(ImagePath, Filetype)
    Images = ImportImageFiles(FileList)

    energy, T, AverageValues, FilterNames, PeakIntensity = getXRaySignal(Images, SettingPath, SettingFile, MaskingName)
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
