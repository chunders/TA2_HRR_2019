#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Tue Aug 27 10:29:13 2019

@author: chrisunderwood

    E Profile analysis
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [8.0,4.0]
import matplotlib.pyplot as plt
import os
from skimage import io
from scipy.signal import medfilt2d
from scipy.optimize import curve_fit

# =============================================================================
# Functions for loading in the file paths
# =============================================================================
def FilesInFolder(DirectoryPath, fileType, starts = ""):
    """ Finds all the files in a directory that are of the filetype given.
    Has the option of putting a prefix too.    
    """
    files = os.listdir(DirectoryPath)
    shots = []
    for i in files:
        if not i.startswith('.') and i.endswith(fileType) and i.startswith(starts):
            shots.append(i)
    return shots

def SpliceArr(arr, start, end):
    """ Slice the strings in an array
    The index to slice to are start and end. Use None if you do not want one of these.
    """
    
    splicedArr = []
    for a in arr:
        item = a[start:end]
        if item.isdigit():          
            splicedArr.append(float(item))
        else:
            splicedArr.append(item)
    return splicedArr

def SplitArr(arr, stringSplit, indexToReturn=0):
    """ Split the strings in an array by the stringSplit char/string
    """
    splitArr = []
    for a in arr:
        item = a.split(stringSplit)[indexToReturn]
        if item.isdigit():
            splitArr.append(float(item))
        else:
            splitArr.append(item)
    return splitArr

def sortArrAbyB(A, B):
    B, A = zip(*sorted(zip(B, A)))
    return np.array(A), np.array(B)


# =============================================================================
# Class for analysis of the eprofile images
# =============================================================================
class eprofile():
    def __init__(self):
        pass
    
    def load_file(self, filePath):
        """ Load the image file
        Convert image to double, as one function required it to not be unit16
        """
        self.image = io.imread(filePath)
        self.image = self.image.astype('float64') 

        
    def crop_to_lanex(self, botLeft,topRight):
        """ Crop the image to the square lanex region
        The botLeft should be the smaller of the two coordinates
        """
        bot = botLeft[0]
        left = botLeft[1]
        top = topRight[0]
        right = topRight[1]
        self.image = self.image[left:right, bot:top]
        
    def simplePointing(self, kernelsize = 7, plotting = True):
        """ Filter the 2D array by a kernel, and then return the pixel coordinates
        of the max value. Average the position if multiple are found
        """
        f_image = medfilt2d(self.image, kernelsize)
        max_coor = np.where(f_image == f_image.max())     
        # These are back to front for some reason
        x_arr = max_coor[1]
        y_arr = max_coor[0]
        if len(x_arr) > 1:
            # print( "Multiple pixels found having the maximum value\nAveraging there positions together")
            x_arr = [np.average(x_arr)]
            y_arr = [np.average(y_arr)]
        if plotting:
            plt.imshow(f_image, vmin = self.image.min(), origin='lower')
            plt.plot(x_arr, y_arr, 'rx-')        
            plt.show()
        return x_arr[0], y_arr[0]
    
    def gausPointing(self, plotting = True):
        """ Returns the coordinates of a gaussian beam center
        Filter the 2D array by a kernel to remove hard hits
        Fit a guassian
        """
        out = self.fit_2DGaus(plotting = plotting)
        if out[0] is not None:
            return out[0], out[1]
        else:
            return None, None
        
    
    def twoD_Gaussian(self,xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
        """ A 2D gaussian function that works with curvefit
        """
        (x, y) = xdata_tuple                                                        
        xo = float(xo)                                                              
        yo = float(yo)                                                              
        a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)   
        b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)    
        c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)   
        g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)         
                            + c*((y-yo)**2)))                                   
        return g.ravel()
        
    def fit_2DGaus(self, plotting = True, view_n_std = 6, kernelsize = 7):
        """ Fit a 2D gaussian to the beam profile
        Returns the xo, yo, sigma_x, sigma_y
        """
        # Flip the image here so the end result is rotated correctly
        f_image = medfilt2d(self.image.T, kernelsize)
        imShape = np.shape(f_image)        
        
        maxCoors = self.simplePointing(plotting = False)
        
        x = np.linspace(0, imShape[0], imShape[0])
        y = np.linspace(0, imShape[1], imShape[1])
        x, y = np.meshgrid(x, y)
        
        #               amplitude, xo, yo 
        initial_guess = [f_image.max(), maxCoors[0], maxCoors[1],
                         # sigma_x, sigma_y, theta, offset
                         imShape[0]/3, imShape[1]/3, 0, 0]
        try:
            popt, pcov = curve_fit(self.twoD_Gaussian, (x, y), f_image.T.ravel(), 
                                   p0=initial_guess)
            
            data_fitted = self.twoD_Gaussian((x, y), *popt)
    
            if plotting:
                plt.pcolormesh(x, y, f_image.T, vmin = self.image.min())
                plt.contour(x, y, data_fitted.reshape(imShape[1], imShape[0]),
                            5, cmap = plt.cm.Spectral)
                xLims = [popt[1] - view_n_std * popt[3], popt[1] + view_n_std * popt[3]]
                yLims = [popt[2] - view_n_std * popt[4], popt[2] + view_n_std * popt[4]]            
                    
                for lim, shape  in zip([xLims, yLims], imShape):
                    if lim[0] < 0:
                        lim[0] = 0
                    if lim[1] > shape -1:
                        lim[1] = shape - 1
                        
                # plt.xlim(xLims)
                # plt.ylim(yLims)
                plt.axes().set_aspect('equal')
                plt.show()
        except RuntimeError:
            print ("The fit has not converged")
            return [None, None, None, None]
        # Return xo, yo, sigma_x, sigma_y        
        return popt[1:5]
    
    def gaussian_divergence(self):
        """ Use fit_2DGaus and average the x and y widths of the gaussian to get 
        an averaged divergence
        """
        out = self.fit_2DGaus(plotting = False)
        if out[2] is not None:
            # The width of a gaussian should be positive, so abs to remove
            # any negatives that may have occurred
            return np.average(abs(out[2:]))
        else:
            return None
        


if __name__ == "__main__":
    folderPath = "/Users/chrisunderwood/Downloads/run005/"

    # Loading the files    
    fileList = FilesInFolder(folderPath, ".TIFF")
    shotNumber = SpliceArr(SplitArr(fileList, ".", indexToReturn=0), 4, None)
    fileList, shotNumber = sortArrAbyB(fileList, shotNumber)
    
    ep = eprofile()
    ep.load_file(folderPath + fileList[52])
    ep.crop_to_lanex([0, 80], [410, 340])
    
    plt.title("Raw image cropped to lanex")
    plt.imshow(ep.image, vmax = 33500, origin='lower')
    plt.colorbar()
    plt.axes().set_aspect('equal')
    plt.show()
    
    x_arr, y_arr = ep.simplePointing(kernelsize=13)
    print ("Simple pointing", x_arr, y_arr)
    
    g_pointing = ep.gausPointing()
    print ("Gaus pointing", g_pointing)

    out = ep.fit_2DGaus()

    div = []
    for file in fileList:
        ep = eprofile()
        ep.load_file(folderPath + file)
        ep.crop_to_lanex([0, 80], [410, 340])
        div.append(ep.gaussian_divergence())
    plt.plot(shotNumber, div)        
    plt.xlabel("Shot Number")
    plt.ylabel("Divergence (Nos Pixels)")
    plt.show()
