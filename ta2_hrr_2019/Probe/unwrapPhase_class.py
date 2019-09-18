#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Wed Jun 19 09:37:31 2019

@author: chrisunderwood

    Unwrap Phase
    
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [6.0,4.0]
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.signal import find_peaks 
import rotateArray

# Load my module of functions
import CUnderwood_Functions3 as func

class unwrapPhase():
    def __init__(self):
        pass

    def load_arrIntoClass(self, phase, plotting = False):
        
        self.phase = phase
        self.phaseShape = np.shape(self.phase)

        if plotting:
            plt.pcolormesh(self.phase)       
            plt.colorbar()
            plt.show()
    
    def load_data(self, loadPath, phaseFile, plotting = False):
        import loadDataToNumpy_class

        ld_image = loadDataToNumpy_class.loadInDataToNumpy(loadPath + phaseFile)
        self.phase = ld_image.loadData()
        self.phaseShape = np.shape(self.phase)
        if plotting:
            plt.pcolormesh(self.phase)       
            plt.colorbar()
            plt.show()

    def unwrap_skimage(self, plotting = False):
        from skimage.restoration import unwrap_phase
        self.unwrappedPhase = unwrap_phase(self.phase, wrap_around=False)
        
        if plotting:
            plt.pcolormesh(self.unwrappedPhase, cmap = plt.cm.seismic,
                               norm = func.MidpointNormalize(midpoint = 0))   
            plt.title("Unwrapped phase using skimage")
            plt.colorbar()
            plt.show()      
            
    def fitting_peak_of_phi_lineout(self, l, plotting = True):
        peaks_pos = find_peaks(l, distance=30, height=0)
        if len(peaks_pos[0]) > 0:
            maxPeak = peaks_pos[1]['peak_heights'].max()
            maxInd = func.nearposn(peaks_pos[1]['peak_heights'], maxPeak)
            guess = [maxPeak, peaks_pos[0][maxInd], 0.1, np.average(l)]
            xl = range(len(l))
            try:
                popt, pcov = curve_fit(gaus, xl,  l, p0 = guess)
                r2 = func.rSquared_from_curveFit(gaus, xl, l, popt)
                if plotting:
                    plt.figure(figsize=(6,3))
                    plt.plot(xl, l, label = "data")
                    plt.plot(xl, gaus(xl, *popt), label = "fit")        
                    plt.title("Detecting Sign of Phase peak")                    
                    plt.legend()
                    plt.show()                
            except RuntimeError:
                print ("Fit failed of gaussian to phase peak")
                r2 = -10
            return r2
        else:
            return -10

            
    def correct_phase_sign(self, plotting = False):
        l = np.sum(self.unwrappedPhase, axis= 1) - np.average(self.unwrappedPhase, axis= 1)
        pos = self.fitting_peak_of_phi_lineout(l)
        neg = self.fitting_peak_of_phi_lineout(-l)
        if neg > pos:
            print ("inverting Phase")
            self.unwrappedPhase *= -1
                
    def mask_pc(self, peakCenterFraction = 0.35, plotMask = True, mask_threshold_level = 0.01, plotLineoutFits = False):
        """ Make a mask around the plasma channel
        For each col in the image, fit a gaussian and exclude the area around the peak
        """
        peakRelHeight = 0.9
        
        # for row in range(self.phaseShape[0]):
        xaxis = np.arange(self.phaseShape[0])
        maskArr = []
        self.pc_centres = []
        maxPhase = self.unwrappedPhase.max()
        for row in range(self.phaseShape[1]):            
            line = self.unwrappedPhase[:,row]            
            xpeaks = find_peaks(line, height = line.max() * peakRelHeight, prominence = 0.6)  
            if len(xpeaks[0]) == 1:
                # Fit a gaussian to the peak to locate it
                guess = [xpeaks[1]['peak_heights'][0], xpeaks[0][0], 0.010, 0]
                try:
                    popt, pcov = curve_fit(gaus, xaxis, line, p0 = guess, 
                                        bounds=([xpeaks[1]['peak_heights'][0] * 0.1, -1000, 0, -1000], 
                                                [xpeaks[1]['peak_heights'][0] * 10, 1000, 1000, 1000])
                                       )
                    # Make the region around pc large
                    popt[2] /= 1.2
                    popt[3] = 0
                    mask = gaus(xaxis, *popt) > mask_threshold_level 
                    maskArr.append(mask)
                    # Append the position of the peak if it is prominent enough
                    if xpeaks[1]['peak_heights'][0] > maxPhase * peakCenterFraction:
                        self.pc_centres.append([row, xpeaks[0][0]])
                    if plotLineoutFits and row % 10 == 0:
                        plt.plot(xaxis, line, "-")
                        plt.plot(xaxis, gaus(xaxis, *popt), "--")
                        
                except RuntimeError:
                    # if the gaussian fit fails
                    maskArr.append(np.ones(self.phaseShape[0], dtype=bool))                    
            else:
                maskArr.append(np.ones(self.phaseShape[0], dtype=bool))
        if plotLineoutFits:
            plt.show()
        self.maskArr = np.array(maskArr).T
        self.pc_centres = np.array(self.pc_centres)
        if plotMask:
            plt.imshow(self.maskArr)
            plt.title("Mask")
            plt.show()        
        return self.maskArr, self.pc_centres

    def fit_background(self):
        self.bg = np.ma.array(self.phase, mask = self.maskArr)
        m = self.phaseShape[0]
        n = self.phaseShape[1]
        X1, X2 = np.mgrid[:m, :n]    
        
        #Regression
        X = np.hstack(   ( np.reshape(X1, (m*n, 1)) , np.reshape(X2, (m*n, 1)) ) )
        X = np.hstack(   ( np.ones((m*n, 1)) , X ))
        YY = np.reshape(self.bg, (m*n, 1))
        
        theta = np.dot(np.dot( np.linalg.pinv(np.dot(X.transpose(), X)), X.transpose()), YY)
        
        self.bg_plane = np.reshape(np.dot(X, theta), (m, n));
        plt.title("Removing background plane")
        plt.imshow(self.bg_plane, vmin = self.unwrappedPhase.min(), vmax = self.unwrappedPhase.max(), 
                   norm = func.MidpointNormalize(midpoint = 0),
                   cmap = plt.cm.seismic)    
        plt.colorbar()
        plt.show()
        
        plt.imshow(self.unwrappedPhase - self.bg_plane, norm = func.MidpointNormalize(midpoint = 0),
                   cmap = plt.cm.seismic)    
        plt.plot(self.pc_centres[:,0], self.pc_centres[:,1], "o-g", lw = 2)
        plt.colorbar()
        plt.show()
        self.unwrappedPhase_with_bg = self.unwrappedPhase
        self.unwrappedPhase = self.unwrappedPhase - self.bg_plane
        
    def rotate_plasmaC_to_Horz(self, start=None, end = None, plot_linfit = True, angle = False):
        if angle:
            rot = rotateArray.rotateFringesToVertical(self.unwrappedPhase)            
            self.unwrappedPhase = rot.applyRotation(angle)            
        else:
            startInd = None
            endInd = None
            if start is not None:
                startInd = func.nearposn(self.pc_centres[:,0], start)            
            if end is not None:
                endInd = func.nearposn(self.pc_centres[:,0], end)      
            print(start, end, startInd, endInd)
            guess = [0, np.average(self.pc_centres[:,1][startInd:endInd])]
            print ("Straight line guess", guess)
            popt, pcov = curve_fit(lin, self.pc_centres[:,0][startInd:endInd], self.pc_centres[:,1][startInd:endInd], p0 = guess)
            if plot_linfit:
                plt.title("Fitting peak centres for rot angle")
                plt.ylim([0, self.phaseShape[0]])
                plt.plot(self.pc_centres[:,0], self.pc_centres[:,1], '.-', label = "Peak Pos")
                plt.plot(self.pc_centres[:,0], lin(self.pc_centres[:,0], *popt), '-', label = "fit")
                plt.legend()
                plt.show()
    
            rot = rotateArray.rotateFringesToVertical(self.unwrappedPhase)
            angle = -rot.rotationAngleFromGradient(popt[0])        
            self.unwrappedPhase = rot.applyRotation(angle)
        return angle
            
       

def lin(x, *params):
    # linear fit function
    m = params[0]
    c = params[1]
    return x *m + c    
    
def gaus(x, *p):
    # Gaus terms
    a = p[0]
    xo = p[1]
    w = p[2]
    c = p[3]
    return a*np.exp( -( (x-xo)**2/2*w**2) ) + c



    
if __name__ == "__main__":
    loadPath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/"
    file = "phaseShift_needsUnwrapping_2.txt"

    unwrap = unwrapPhase()
    unwrap.load_data(loadPath, file, False)
    unwrap.unwrap_skimage(plotting = True)
    
    np.savetxt(loadPath + "phaseShift_unwrapped.txt", unwrap.unwrappedPhase)
    unwrap.correct_phase_sign()
    maskArr, centers = unwrap.mask_pc()
    c = np.array(centers)
    # plt.plot(c[:,0], c[:,1])    
    
    unwrap.fit_background()
    unwrap.rotate_plasmaC_to_Horz(start = 60)

     # np.savetxt('/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/phaseShift_unwrapped_2.txt',
     #            unwrap.unwrappedPhase)
    
    
    