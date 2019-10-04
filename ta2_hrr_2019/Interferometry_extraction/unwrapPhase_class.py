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
from fourier_mask_class import fourier_filter_for_mask

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
        pos = self.fitting_peak_of_phi_lineout(l, plotting)
        neg = self.fitting_peak_of_phi_lineout(-l, plotting)
        if neg > pos:
            print ("inverting Phase")
            self.unwrappedPhase *= -1
                
    def mask_pc(self, plotMask = True, plotPlasmaChannel = False,  mask_threshold_level = 0.55):
        """ Make a mask around the plasma channel
        Create a fourier mask
        """


        mask = fourier_filter_for_mask(self.unwrappedPhase)
        mask.create_mask(mask_percentage = mask_threshold_level)
        self.maskArr = mask.mask
        if plotMask:
            plt.imshow(self.maskArr)
            plt.title("Mask")
            plt.show()    
            

        plasma_channel = mask.antimask * self.unwrappedPhase
        if plotPlasmaChannel:
            plt.imshow(plasma_channel)
            plt.title("masked plasma channel")
            plt.show()  
        
        self.pc_centres = []
        for row in range(self.phaseShape[1]):            
            line = plasma_channel[:,row]       
            if line.max() > 0:
                xpeak = func.nearposn(line, line.max())
                self.pc_centres.append([row, xpeak])
        self.pc_centres = np.array(self.pc_centres)
            
        return self.maskArr, self.pc_centres

    def fit_background(self, plotting = False):
        """ Fit the background.
        Fit a 2D plane and subtract
        """
        self.bg = np.ma.array(self.phase, mask = self.maskArr)
        m = self.phaseShape[0]
        n = self.phaseShape[1]
        X1, X2 = np.mgrid[:m, :n]    
        
        #Regression
        X = np.hstack(   ( np.reshape(X1, (m*n, 1)) , np.reshape(X2, (m*n, 1)) ) )
        X = np.hstack(   ( np.ones((m*n, 1)) , X ))
        YY = np.reshape(self.bg, (m*n, 1))
        
        theta = np.dot(np.dot( np.linalg.pinv(np.dot(X.transpose(), X)), X.transpose()), YY)
        
        self.bg_plane = np.reshape(np.dot(X, theta), (m, n))
        if plotting:
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

        if False:        
            # Plot in 3d to check the background is subtracted properly
            from mpl_toolkits.mplot3d import Axes3D        
            jet = plt.get_cmap('jet')
            fig = plt.figure(figsize = (8,8))
            ax = fig.add_subplot(1,1,1, projection='3d')        
            ax.plot_surface(X1,X2,self.bg_plane)
            ax.plot_surface(X1,X2,self.bg, rstride = 60, cstride = 10, cmap = jet, linewidth = 0)
            ax.view_init(elev=15., azim=-35)
            ax.set_zlim(-1,1)
            plt.show()

        
        
    def rotate_plasmaC_to_Horz(self, start=None, end = None, plot_linfit = True, plot_rotation = False,  angle = False):
        if angle:
            rot = rotateArray.rotateFringesToVertical(self.unwrappedPhase)    
            print (plot_rotation)        
            self.unwrappedPhase = rot.applyRotation(angle, plotting = plot_rotation)            
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
            self.unwrappedPhase = rot.applyRotation(angle, plotting = plot_rotation)

        print ("Angle Rotated: ", angle)
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
    
    unwrap.fit_background(plotting=True)
    unwrap.rotate_plasmaC_to_Horz(start = 60, end = 300, plot_rotation = True)

     # np.savetxt('/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/phaseShift_unwrapped_2.txt',
     #            unwrap.unwrappedPhase)
    
    
    plt.pcolormesh(unwrap.unwrappedPhase_with_bg * unwrap.maskArr)
    plt.colorbar()
    
    print ( np.average(unwrap.unwrappedPhase_with_bg * unwrap.maskArr) )
    print ( np.average(unwrap.unwrappedPhase_with_bg ) )