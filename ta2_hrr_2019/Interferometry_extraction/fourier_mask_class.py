#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Thu Oct  3 13:38:31 2019

@author: chrisunderwood
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [5.0,4.0]
import matplotlib.pyplot as plt

# Load my module of functions
import sys
sys.path.insert(0, '/Users/chrisunderwood/Documents/Python/')
import CUnderwood_Functions3 as func

class fourier_filter_for_mask():
    def __init__(self, image):
        self.image = image
        self.shape = np.shape(image)
    
    def wrapImage(self, image):
       wrappedImage = np.zeros_like(image)
       shape = np.shape(image)
       for i in range(shape[0]):
           for j in range(shape[1]):
               # Wrapping image
               i_w = (i + shape[0]//2)%shape[0]
               j_w = (j + shape[1]//2)%shape[1]
               # print 
               wrappedImage[i_w][j_w] = image[i][j]
       return wrappedImage
   
    def create_mask(self, x = 0.9, y = 0.9, mask_percentage = 0.5):
        
        self.F_image= np.fft.fft(self.image)
        self.crop_in_fourier_space(x, y)
        self.filtered_image = np.fft.ifft(self.F_image)
        self.norm_filter = func.normaliseArr(np.real(self.filtered_image))
        self.mask = self.norm_filter < mask_percentage
        self.antimask = self.norm_filter > mask_percentage
        
        
        
    def crop_in_fourier_space(self, x = 0.5, y = 0.5):
        y = int(self.shape[0]//2 * y)
        x = int(self.shape[1]//2 * x)
        print (x, y)
    
        for i in range(self.shape[0]//2-y, self.shape[0]//2+y):
            for j in range(self.shape[1]//2-x, self.shape[1]//2+x):
                self.F_image[i][j] = 0
                
    def plot(self, image, cbarOn = True, title = ""):
        plt.pcolormesh(image)
        if cbarOn:
            plt.colorbar()
        plt.title(title)
        plt.show()        
    
if __name__ == "__main__":
    loadPath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/"    
    phi = np.loadtxt(loadPath + "phaseShift_unwrapped.txt")    
        
    
    mask = fourier_filter_for_mask(phi)
    mask.create_mask(mask_percentage = 0.55)
    
    mask.plot(mask.image)
    mask.plot( abs(mask.filtered_image))
    mask.plot( abs(mask.norm_filter))
    mask.plot(mask.mask, cbarOn = False)
    mask.plot(mask.mask * mask.image, title = "Bg")
