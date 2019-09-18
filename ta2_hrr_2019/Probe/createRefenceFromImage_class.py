#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Wed Jun 19 16:36:41 2019

@author: chrisunderwood

    Create reference image
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [6.0,4.0]
import matplotlib.pyplot as plt

# Load my module of functions
import CUnderwood_Functions3 as func


class createReference():
    
    def __init__(self):
        pass
    
    def loadData(self, filePath, line_of_pchannel, width_pchannel):
        import loadDataToNumpy_class
        ld = loadDataToNumpy_class.loadInDataToNumpy(filePath)
        self.im = ld.loadData()
        self.imShape = np.shape(self.im)
        self.line_of_pchannel = line_of_pchannel
        self.width_pchannel = width_pchannel     
        
    def load_arrIntoClass(self, arr, line_of_pchannel, width_pchannel):
        self.im = arr
        self.imShape = np.shape(self.im)
        self.line_of_pchannel = line_of_pchannel
        self.width_pchannel = width_pchannel             
        
    def createAveragedLineout(self):
    # Creates an average lineout of the image, not including the plasma channel.
    # This assumes that the region up to the plasma channel is big enough to create a
    # reference from.
        indexToCropTo = self.line_of_pchannel - self.width_pchannel
        try:
            self.lineout = self.im[:indexToCropTo].mean(axis=0)
        except:
            self.lineout = np.average(self.im[:indexToCropTo], axis = 0)        
            
    def create_ref_From_lineouts(self, plotting=False):
    # Work out where the phase shifts are from the fringes
    # This is calcuated from the lineouts, so assumes
    # that the fringes are vertical
        self.reference_CreatedFromLineouts = []
        for i in range(self.imShape[0]):
            self.reference_CreatedFromLineouts.append(self.lineout)
        self.reference_CreatedFromLineouts = np.array(self.reference_CreatedFromLineouts)
        
        if plotting:
            plt.title('The refernce image')
            plt.imshow(self.reference_CreatedFromLineouts)
            plt.show()
            
    def display_marking_PlasmaChannel(self, centeredOnZero=True):

    # Plot image, with indicators of where I have said the plasma channel is.
        plt.hlines(self.line_of_pchannel, 0, self.imShape[1])
        plt.hlines(self.line_of_pchannel + self.width_pchannel, 0, self.imShape[1],
                   linestyle = 'dashed')
        plt.hlines(self.line_of_pchannel - self.width_pchannel, 0, self.imShape[1],
                   linestyle = 'dashed')
        plt.title('Displaying image with the region of the plasma\nchannel marked')
        if centeredOnZero:
            plt.imshow(self.im, norm = func.MidpointNormalize(midpoint = 0), cmap = plt.cm.seismic)
            plt.colorbar()
            plt.show()
        else:   
            plt.imshow(self.im)
            plt.show()
        


if __name__ == "__main__":
    filePath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/im_bgRemoved_croppedROI.txt"
    savePath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/"

    ref = createReference()
    
    line_of_pchannel = 160
    width_pchannel = 18

    ref.loadData(filePath, line_of_pchannel, width_pchannel)
    ref.display_marking_PlasmaChannel()

    ref.createAveragedLineout()
    
    ref.create_ref_From_lineouts(plotting=True)
    np.savetxt(savePath + "reference.txt",  ref.reference_CreatedFromLineouts)
