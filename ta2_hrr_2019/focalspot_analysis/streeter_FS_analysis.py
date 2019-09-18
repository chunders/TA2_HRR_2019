#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Thu Jul 18 16:46:25 2019

@author: chrisunderwood
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [8.0,8.0]
import matplotlib.pyplot as plt

# Load my module of functions
import sys
sys.path.insert(0, '/Users/chrisunderwood/Documents/Python/')
import CUnderwood_Functions3 as func

sys.path.insert(0, '/Users/chrisunderwood/Documents/Experimental_Tools/FocalSpot_Analysis/')
from FocalSpot import focal_spot

from skimage import io


def spatial_cal_with_grating():
    """
    Spatial calibration
    
    17.5 cm with 25 slits
    d = 0.007 m
        eq:
            d \sin \theta_{m} = m\lambda
    
    np.deg2rad(np.arcsin( (1 * lambda_l) / d ))        
    """
    d = 0.007
    fl = 1
    lambda_l = 800e-9
    Space_of_diffraction_peaks = lambda_l * fl / d
    # = 0.00011428571428571428
    return (Space_of_diffraction_peaks)

def calibration_with_grating():
    
    peakSpacing = spatial_cal_with_grating()
    print ("Predicted spacing of spots")
    print (peakSpacing)
    
    # Positions of the spots from the image, extracted with ImageJ
    yPos = np.array([55.97, 170.83, 290])
    xPos = np.array([34.12, 147.83, 260.47])
    yGap = yPos[1:] - yPos[:-1]
    xGap = xPos[1:] - xPos[:-1]
    
    aveSpace = [np.average(yGap), np.average( xGap)]
    print ("Measured spacing in pixesl")
    print (aveSpace )
    print ()
    
    pixelSize = peakSpacing /  np.average(aveSpace)
    print ("pixels_per_m ", pixelSize)
    return pixelSize

folderPath = "/Volumes/GoogleDrive/My Drive/2019_Streeter_TA2/FocalSpot/"
# laserShots = "20190718_LP/"
# backgroundShots = "20190718_bg/"

# laserShots = "20190719_LP/"
# backgroundShots = "20190719_bg/"

laserShots = "20190719_postOpt/"

pixelSize = calibration_with_grating()

# =============================================================================
# Loading in the files
# =============================================================================
laserFiles = func.FilesInFolder(folderPath + laserShots, ".tiff")

if False:
    backgroundFiles = func.FilesInFolder(folderPath + backgroundShots, ".tiff")

    # Background summation
    background_sums = []
    for i, b in enumerate(backgroundFiles[:]):
        print ("Background File: ", b)
        image = io.imread(folderPath + backgroundShots + b)
        fs = focal_spot(image, plot_raw_input = False)    
        background_sums.append(fs.return_energy())
        if i == 0:
            background_aveImage = fs.im
        else:
            background_aveImage += fs.im
    background_aveImage = background_aveImage / (i+1)
else:
    background_aveImage = np.loadtxt(folderPath + "Background_average_20190719_bg.txt")


# Laser focal spot calculations
data_accumulative = []
for l in laserFiles[:]:
    print ("\nLaser file: ", l)

    image = io.imread(folderPath + laserShots + l)
    fs = focal_spot(image, plot_raw_input = False, 
                    backgroundImage = background_aveImage) 
    
    d = fs.fit_2DGaus(umPerPixel = pixelSize * 1e6, plotting = False)
    fwhm = fs.fwhm_gaussian(d)
    peakIntensity = (fs.im/sum(fs.im)/(pixelSize*1e-2)**2 * 4.92).max()
    data_accumulative.append([d[0], d[1], 0.5*(d[0]+ d[1]),  fs.return_energy(), peakIntensity])

# Plot the 
data_accumulative = np.array(data_accumulative)
f, ax = plt.subplots(nrows = 3, sharex = True, gridspec_kw={'height_ratios': [3, 1, 1]})
ax[0].plot(data_accumulative[:,0], label = "x size")
ax[0].plot(data_accumulative[:,1], label = "y size")
ax[0].plot(data_accumulative[:,2], label = "ave size")

ax[0].legend()
ax[1].plot(data_accumulative[:,3], 'r', label = "energy")
ax[0].set_ylabel(r"Spot size ($\mu m$)")
ax[1].set_ylabel(r"Energy (arb units)")


ax[2].set_ylabel("Peak Intensity\n" + r"($W cm^{-2}$)")
ax[2].plot(data_accumulative[:,4], label = "peakIntensity")
ax[2].set_xlabel("Shot Number")
func.saveFigure(folderPath + laserShots[:-1] + "_focalspots.png")
plt.show()

ave_size = np.average(data_accumulative[:,2])
std_size = np.std(data_accumulative[:,2])
ave_energy = np.average(data_accumulative[:,3])
ave_intensity = np.average(data_accumulative[:,4])

outDict = {}
outDict['ave_size'] = ave_size
outDict['std_size'] = std_size
outDict['ave_energy'] = ave_energy
outDict['ave_intensity'] = ave_intensity
outDict['data'] = data_accumulative


func.saveDictionary(folderPath + laserShots[:-1] + "_focalspots.json", outDict)

aveBackgroundCount = np.average(background_sums)
print ("\n\nAverage FWHM (um)\t {:.2e} +/-{:.2e} ".format( ave_size, std_size),
       "\nAve Energy in shot (pixel counts)\t", ave_energy - aveBackgroundCount,
       "\nAve Intensity \t {:.2e}".format( ave_intensity))


# # I = I / np.sum(I) / (dx * dx) * 4.92
# plt.imshow(fs.im/sum(fs.im)/(pixelSize*1e-2)**2 * 4.92)
# plt.colorbar()
