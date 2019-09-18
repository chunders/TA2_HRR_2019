#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Wed Jun 19 09:38:20 2019

@author: chrisunderwood

    Phase shift to Density
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [8.0, 4.0]
import matplotlib.pyplot as plt

# Load my module of functions
import CUnderwood_Functions3 as func
from scipy.optimize import curve_fit
import abel


class deltaPhaseToDensity():
    
    def __init__(self, laserwavelength_m, mPerPix):
        self.constants(laserwavelength_m, mPerPix)
        
    def constants(self, laserwavelength_m, mPerPix):
        self.lambda_l = laserwavelength_m
        self.mu_0 =  12.57e-7
        self.m_e = 9.11e-31
        self.q_e = 1.6e-19
        self.e_0 = 8.85e-12
        self.c = 3e8
        self.sizePerPixel = mPerPix
        
        self.w_0 = (2 * np.pi * self.c ) / self.lambda_l
        self.n_c = (self.w_0**2 * self.m_e * self.e_0) / self.q_e**2   

        self.r_e = 2.81794e-15        
    
    def load_data(self, filePath, plotting = False):
        import loadDataToNumpy_class
        ld = loadDataToNumpy_class.loadInDataToNumpy(filePath)
        self.phase = ld.loadData()
        self.phaseShape = np.shape(self.phase)
        if plotting:
            plt.figure(figsize = (8,4))
            plt.title("Input phase")
            plt.imshow(self.phase, cmap = plt.cm.seismic,
                       norm = func.MidpointNormalize(midpoint = 0))
            plt.colorbar()
            plt.show()
            
    def load_arrIntoClass(self, arr):
        self.phase = arr
        self.phaseShape = np.shape(self.phase)             
            
    
    def inverse_abel_transform(self, plotting = False, method = "hansenlaw"):
        # Flip the image, so the abel transform work
        image = self.phase.T
        
        # Using the inbuilt gaussian method to find the axis
        self.inverse_abel = abel.Transform(image,
                                      #center = (50, 200),
                                      method =  method, 
                                      center = "gaussian",
                                      center_options = {'axes' : 1, "verbose":True},
                                      direction='inverse', verbose = True).transform.T        
            
        if plotting:
            plt.figure(figsize = (8,4))            
            plt.imshow(self.inverse_abel, cmap = plt.cm.seismic, 
                   norm = func.MidpointNormalize(midpoint = 0),
                   ) 
            cbar = plt.colorbar()
            cbar.set_label("Raw Abel Transform")
            plt.title("Inverse Abel Transform")
            plt.xlabel("Pixels")
            plt.ylabel("Pixels")        
            plt.show()

            lineout_ave = np.average(self.inverse_abel[ 10:-10, :], axis = 0)
            plt.figure(figsize = (8,4))
            
            plt.plot(np.arange(len(lineout_ave)) * self.sizePerPixel * 1e3, lineout_ave)
            plt.xlabel("Distance (mm)")
            plt.show()

    def convert_Inverse_Abel_to_Ne(self, plotting = True, pixelsAroundPlasmaChannel = 10):
        perCm3 = False
        if perCm3: self.inverse_abel *= 1e-6
        
        # Using N_e = 1/(r_e * l) * pyabelInverse, from chat with Rob Shaloo
    
        self.n_e = self.inverse_abel / (self.r_e *  self.lambda_l * self.sizePerPixel)
        
        # Take an average cropping to the center, this could be cleverer
        
        print ("Taking average lineout in region of size {}mm around axis".format(2 * pixelsAroundPlasmaChannel * self.sizePerPixel *1e3))
        lineout_ave = np.average(self.n_e[ self.phaseShape[0]//2 - pixelsAroundPlasmaChannel: self.phaseShape[0]//2 + pixelsAroundPlasmaChannel, :], 
                                 axis = 0)
        if plotting:
            f, ax = plt.subplots(nrows=2, sharex = True, figsize = (8,6))
            ax[0].set_title("Number Density")
            im1 = ax[0].pcolormesh( 
                    np.arange(self.phaseShape[1]) * self.sizePerPixel *1e3,
                    np.arange(self.phaseShape[0]) * self.sizePerPixel *1e3 - self.phaseShape[0] * 0.5 * self.sizePerPixel *1e3,
                    self.n_e, cmap = plt.cm.seismic,
                    norm = func.MidpointNormalize(midpoint = 0) )
            for height in [self.phaseShape[0]//2 - pixelsAroundPlasmaChannel - self.phaseShape[0] * 0.5, self.phaseShape[0]//2 + pixelsAroundPlasmaChannel - self.phaseShape[0] * 0.5]:
                print (self.phaseShape, height)
                ax[0].hlines(height * self.sizePerPixel *1e3, 0, self.phaseShape[1]* self.sizePerPixel *1e3)
            
            cax = f.add_axes([0.95, 0.25, 0.05, 0.5])
            plt.colorbar(im1, cax = cax)
            ax[1].plot(np.arange(len(lineout_ave)) * self.sizePerPixel *1e3, lineout_ave)
            ax[1].set_xlabel("Distance (mm)")
            ax[1].set_ylim([0, None])
            plt.show()        
            
        return self.n_e, np.c_[np.arange(len(lineout_ave)) * self.sizePerPixel *1e3, lineout_ave]

if __name__ == "__main__":
    filePath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/"
    fileName = "phaseShift_unwrapped_2.txt"
    rhoCalc = deltaPhaseToDensity(laserwavelength_m = 800e-9, mPerPix = 42.2e-6)
    rhoCalc.load_data(filePath + fileName, plotting=True)
    # rhoCalc.rotatePlasmaChanel_to_axis(plotting=True, plottingLineouts = True)
    # rhoCalc.rotateByGivenAngle(-2.5, plotting=True)
    widen = 15
    midpoint = 48

    # rhoCalc.cropPhaseToCenter([midpoint - widen, midpoint + widen, 0, 300])
    # rhoCalc.shift_to_centre()
    rhoCalc.inverse_abel_transform(plotting = False)
    
    _ = rhoCalc.convert_Inverse_Abel_to_Ne()
    
    