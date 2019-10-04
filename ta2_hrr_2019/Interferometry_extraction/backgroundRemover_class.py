#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Wed Jun 19 09:36:27 2019

@author: chrisunderwood
    
    Background Remover and crop to ROI
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [6.0,4.0]
import matplotlib.pyplot as plt

# Load my module of functions
import CUnderwood_Functions3 as func



class backgroundRemoverAndCrop():
    
    def __init__(self, plotting = False):
        self.plotting = plotting


    def LoadData(self, fileName):
        """ Load the data using the class
        """
        import loadDataToNumpy_class
        ld = loadDataToNumpy_class.loadInDataToNumpy(fileName)
        
        self.im = ld.loadData()
        self.shape = np.shape(self.im)
        if self.plotting:
            plt.imshow(self.im, vmin = np.average(self.im) - (self.im.max() - np.average(self.im)))
            plt.colorbar()
            plt.show()
        
    def load_arrIntoClass(self, arr):
        """ Directly take a numpy array
        """
        self.im = arr
        self.shape = np.shape(self.im)
        if self.plotting:
            plt.imshow(self.im, vmin = np.average(self.im) - (self.im.max() - np.average(self.im)))
            plt.colorbar()
            plt.show()
     
        
    def blur(self, image, sigSize):
        """ Blur the image with a gaussain filter to create a bg to remove
        This removes large order deffects in the beam, but should leave the fringes.
        Make sure that sigSize is large enough
        """
    # This cannot be run after the mask has been applied (???What mask)
        import scipy.ndimage as ndimage

        self.im_gblur = ndimage.gaussian_filter(image,
                                sigma=(sigSize, sigSize), order=0)
        if self.plotting:
            plt.imshow(self.im_gblur, vmin = 33000)
            plt.show()
        
    def sub_blurredIm(self, bot = None, top= None, left= None, right= None):
        """ Subtract the blurred image from the raw image.
        If bot, top, left, right are given, also crops to a larger region of interest
        """
        self.im_bgrm = self.im - self.im_gblur
        
        imageRange = 1000
        
        if self.plotting:
            if not None in [bot, top, left, right]:
                plt.vlines(left, bot, top)
                plt.vlines(right, bot, top)
                plt.hlines(top, left, right)
                plt.hlines(bot, left, right)
            
            plt.imshow(self.im_bgrm, vmin=-imageRange,  vmax=imageRange)
            plt.colorbar()
            plt.show()
        return self.im_bgrm[bot:top, left:right]

        

if __name__ == "__main__":
    #Load original data
    LoadPath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/verticalFrings.txt"
    savePath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/"
        # Call the class
    bc = backgroundRemoverAndCrop(plotting = True)
    bc.LoadData(LoadPath)
    bc.blur(bc.im, 11)
    bot = 40
    top = 300
    left = 100
    right = 500
    out = bc.sub_blurredIm(bot, top, left, right)    
    # out = bc.sub_blurredIm()    

    
    np.savetxt(savePath + "im_bgRemoved_croppedROI.txt", out)
    plt.imshow(out)
    plt.colorbar()