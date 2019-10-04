#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Wed Jun 19 09:33:32 2019

@author: chrisunderwood

    Rotate array to vertical
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [6.0,4.0]
import matplotlib.pyplot as plt

# Load my module of functions
import CUnderwood_Functions3 as func

def mode(arr):
    if arr==[]:
        return None
    else:
        return max(set(arr), key=arr.count)


class rotateFringesToVertical():
    def __init__(self, image):
       self.im_raw = image 
       # im_raw should be an np.array, otherwise this can crash the code. 
       self.im_raw = np.array(self.im_raw)       
       self.im_shape = np.shape(image)
        
    def applyRotation(self, angleDeg, plotting = True, verbose = False):    
    # Rotate the image, by the amount given by the angle of the fringes
    # Set the rotated image to be the new image.
        import cv2
        
        M = cv2.getRotationMatrix2D((self.im_shape[0]/2,self.im_shape[1]/2), 
                                    angleDeg, # Angle of rotation in Degrees
                                    1   #Scale
                                    )

        self.im_vert = cv2.warpAffine(self.im_raw ,M,(self.im_shape[1] ,self.im_shape[0]))
        if verbose:
            print ('\tRotation of image')
            print ('tRotation matrix, \n' ,  M)
            print ('Angle Rads' , np.deg2rad(angleDeg))
            print ('Angle Degrees', angleDeg)

        if plotting:
            # Show that the rotation of the image is correctly done - straightline on rotated one should be with fringes
            f, ax = plt.subplots(ncols = 2, figsize=(10, 4))
            im1 = ax[0].imshow(self.im_vert, vmin = np.average(self.im_vert) - (self.im_vert.max() - np.average(self.im_vert)))
            ax[0].vlines(self.im_shape[1]/2 , 0, self.im_shape[0]-1)
            ax[0].set_title('Rotated Image')
            plt.colorbar(im1, ax = ax[0])
            
            im2 = ax[1].imshow(self.im_raw, vmin = np.average(self.im_vert) - (self.im_vert.max() - np.average(self.im_vert)))
            plt.colorbar(im2, ax = ax[1])            
            ax[1].vlines(self.im_shape[1]/2, 0, self.im_shape[0]-1)
            ax[1].set_title('Raw Image')
            plt.show()
            if False:
                #Check that the rotated image is actually different to the orgiinal one
                plt.title('Rotated subtracted from standard')
                plt.imshow(self.im_raw - self.im_vert)
                plt.show()
        #Return the image, which now is rotated such that the fringes are vertical
        return self.im_vert
    

    
    def findRotation(self, startIndex = 0, endIndex = 1000, horzStart = 200, horzEnd = 450, 
                     plotting = False):
        # The angle of the image fringes is returned
        # Create lineout region 
        self.lineoutPixelHeight = 5
        numberOfLineouts = (endIndex - startIndex) // self.lineoutPixelHeight # Find the number of rows of pixel width 5 we can fit in
        print ("Number of lineouts in region ", startIndex, endIndex, numberOfLineouts)
        heights_2DBoxes, width_AreaOfFringes = self.createLineoutRegion_unperturbed(startIndex,
                                                            numberOfLineouts, horzStart, horzEnd)
        # print (heights_2DBoxes, width_AreaOfFringes)
        # Make linesouts
        self.Create_lineouts_of_2Dimage_at_given_heights_Over_width_given(heights_2DBoxes,
                                    width_AreaOfFringes, plotting = plotting)
        
        peaks_positionOfFringes, LineOutHeights = self.locating_Peaks_in_lineouts()
        m_ave, c_ave = self.finding_equation_ofStraightline_for_each_Fringe(peaks_positionOfFringes, LineOutHeights,
                                                                            plotting = plotting)
        # Returns the angles in degrees
        return self.rotationAngleFromGradient(m_ave)
        
    
    def createLineoutRegion_unperturbed(self, startIndex, numberOfLineouts, horzStart, horzEnd):
    # Create a region where the fringes are unperturbed.
    # This region is what we will used to find the rotation required.
        heights_2DBoxes = []
        for i in range(numberOfLineouts):
            heights_2DBoxes.append([startIndex + i * self.lineoutPixelHeight, 
                                    startIndex + i * self.lineoutPixelHeight + self.lineoutPixelHeight])
        width_AreaOfFringes = [horzStart, horzEnd]
        self.horzStart = horzStart
        
        return heights_2DBoxes, width_AreaOfFringes
    
    def Create_lineouts_of_2Dimage_at_given_heights_Over_width_given(self, 
                                            height_boxes, horrizontal_startEnd, plotting=False):
    # Creates lineouts of a 2Dimage over a set of height boxes, all over the same
    # widths.
        relHeightsPeaks = 0.3
        distance = 5
        self.relHeights = relHeightsPeaks
        self.findDistance = distance
        self.lineouts = {}
        if plotting: plt.pcolormesh(self.im_raw, cmap = plt.cm.Blues_r)
        
        for heights in height_boxes:
            croppedIm = self.im_raw[heights[0]:heights[1], horrizontal_startEnd[0]:horrizontal_startEnd[1]]
            index = (heights[0]+heights[1] )* 0.5
            self.lineouts[index] = np.sum(croppedIm, axis = 0)
            
            if plotting:
                # print (heights)
                plt.hlines(heights[0], horrizontal_startEnd[0], horrizontal_startEnd[1])
                plt.pcolormesh(range(horrizontal_startEnd[0], horrizontal_startEnd[1]),
                               range(heights[0], heights[1]),
                               croppedIm)
                plt.hlines(heights[1], horrizontal_startEnd[0], horrizontal_startEnd[1])    
                
                # l = np.sum(croppedIm, axis = 0)
                # peaks = find_peaks(l, rel_height=self.relHeights, distance = self.findDistance)
                # peaksLocations = np.array(peaks[0])
                # plt.plot(horrizontal_startEnd[0] + peaksLocations, np.ones_like(peaksLocations) * (heights[1] + heights[0])*0.5,
                #          'rx--')
            
        if plotting:
            plt.title("The lineouts across the fringes")
            plt.axes().set_aspect('equal', 'datalim')
            plt.show()  
            
    def locating_Peaks_in_lineouts(self, plotting = False):
    # From the lineouts at each height, find the peaks of each sine wave.
    # Return the peak locations at that height
        from scipy.signal import find_peaks

        # The central heights of the lineouts
        centralHeights_lineouts = list(self.lineouts)
        centralHeights_lineouts.sort()
        
        # Find the peaks of the centralHeights_lineouts
        LineOutHeights = []
        peaks_positionOfFringes = []
        for line in centralHeights_lineouts:
            l = self.lineouts[line] # The lineout data
            peaks = find_peaks(l, rel_height=self.relHeights, distance = self.findDistance)
            peaksLocations = np.array(peaks[0]) + self.horzStart # The peak locations appended to a list

            # The heights and peaks are appended to lists.
            LineOutHeights.append(line) 
            peaks_positionOfFringes.append(peaksLocations)
            if plotting: plt.plot(peaksLocations, np.ones_like(peaksLocations) * line, 'o-')
        if plotting: plt.show()

        # For the next stage of the analysis I need a rectangular array of data
        # Therefore if an incorrect number of peaks has been found in a lineout
        # it needs to be removed.
        lenOfPeaks = []        
        for row in peaks_positionOfFringes:
            lenOfPeaks.append(np.shape(row))
        numberOfPeaks = mode(lenOfPeaks)
        
        # Remove the lineouts without enough peaks
        peaks_positionOfFringes_cropped = []
        LineOutHeights_cropped = []
        for i, row in enumerate(peaks_positionOfFringes):
            if len(row) == numberOfPeaks[0]:
                peaks_positionOfFringes_cropped.append(row)
                LineOutHeights_cropped.append(LineOutHeights[i])        
        peaks_positionOfFringes_cropped = np.array(peaks_positionOfFringes_cropped)    
        
        # Return the cropped arrays.
        return peaks_positionOfFringes_cropped, LineOutHeights_cropped
    
    
    def finding_equation_ofStraightline_for_each_Fringe(self, peaks_positionOfFringes, heights,
                                                        plotting = True):
    # Finding the equation of the straight line that goes through each peak/fringe
    # from this angle the rotation of the image can be calculated
        from scipy.optimize import curve_fit
        
        if plotting:    
            plt.figure(figsize=(10,7))
            color=iter(plt.cm.rainbow(np.linspace(0,1, np.shape(peaks_positionOfFringes)[1] )))
        Average_Fit = []
        
        assert len(np.shape(peaks_positionOfFringes)) == 2, "The length of each row in peaks_positionOfFringes is not the same\nTherefore the fringes are not sharp enough for each lineout to have a peak found"
        
        for i in range(np.shape(peaks_positionOfFringes)[1]):
            p0 = [1e5, -1e6]
            popt, pcov = curve_fit(lin, heights, peaks_positionOfFringes[:,i], p0=p0)  
            y_for_fitting = np.linspace(0,self.im_shape[0])
            Average_Fit.append(popt)
            if plotting: 
                # print (popt)
                c=next(color)
                plt.plot(peaks_positionOfFringes[:,i], heights, 'x-', color=c)
                plt.plot(lin(y_for_fitting, *popt), y_for_fitting, '-', color=c)
        Average_Fit = np.array(Average_Fit)
        m_ave = np.average(Average_Fit[:,0])
        c_ave = np.average(Average_Fit[:,1])             
        print ('Ave gradient and intercept', m_ave, c_ave)
        
        if plotting:
            plt.pcolormesh(range(self.im_shape[1]),
                          range(self.im_shape[0]), 
                          self.im_raw)
            plt.xlim([0, self.im_shape[1]])
            plt.ylim([0, self.im_shape[0]])     
            plt.show()
        return m_ave, c_ave
    
    def rotationAngleFromGradient(self, m_ave):
        if m_ave != 0.0:
            # print ('Rotating Image, as angle of fringes is not zero')
            angle = -np.rad2deg(np.arctan(m_ave))
        else:
            angle = 0.
        return angle
    
    
    
    
def plotImage(im, name):
    plt.title(name)
    plt.imshow(im)
    plt.show()    
    

def lin(x, *params):
    # linear fit function
    m = params[0]
    c = params[1]
    return x *m + c
    
    
if __name__ == "__main__":
    # Load the data
    import loadDataToNumpy_class
    rootExperimentFolder = "/Volumes/GoogleDrive/Shared drives/Murphy Group/GeminiRR_2019/"
    filePath = rootExperimentFolder + "20190208/20190208r011/20190208r011s001_Probe_Interfero.tiff"
    ld = loadDataToNumpy_class.loadInDataToNumpy(filePath)
    im_raw = ld.loadData()
    plotImage(im_raw, "Raw Data")
    
    # Run the rotateArray class
    rotAnal = rotateFringesToVertical(im_raw)
    angle = rotAnal.findRotation(startIndex = 5, endIndex = 200,
                                 plotting = True)
    # print (angle)
    rotatedImage = rotAnal.applyRotation(angle)
    
    savePath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/"
    np.savetxt(savePath + "verticalFrings.txt", rotatedImage)
    
