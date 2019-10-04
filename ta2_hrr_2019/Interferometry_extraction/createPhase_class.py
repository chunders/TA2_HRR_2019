#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Thu Aug 29 10:56:50 2019

@author: chrisunderwood
    Recreate Phase
    
    Fourier space variables will start with F_
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [8.0,6.0]
import matplotlib.pyplot as plt

# Load my module of functions
import sys
sys.path.insert(0, '/Users/chrisunderwood/Documents/Python/')
import CUnderwood_Functions3 as func

class phaseShift():
    
    def __init__(self):
        pass
    
    def plot_raw_input(self):
        """ Plot the raw image with a colorbar
        """
        if self.refExists:
            f, ax = plt.subplots(ncols = 2)
            im1 = ax[0].pcolormesh(self.im , cmap = plt.cm.seismic)
            im2 = ax[1].pcolormesh(self.ref, cmap = plt.cm.seismic)
            plt.colorbar(im1, ax = ax[0])
            plt.colorbar(im2, ax = ax[1])                
            ax[0].set_title("Im")
            ax[1].set_title("Ref")                
            plt.show()
        else:
            plt.imshow(self.im , cmap = plt.cm.seismic)
            plt.colorbar()
            plt.show()
    
    def load_arrIntoClass(self, im, ref = None, plotting = False):
        """ Load data into the class from arrays
        """
        self.im = im
        self.imShape = np.shape(self.im)
        if ref is not None:
            self.ref = ref
            self.refExists = True
        else:
            self.refExists = False

        if plotting:
            self.plot_raw_input()
    
    def load_data(self, loadPath, imFile, refFile = None, plotting = False):
        """ Load data into the class from file
        """        
        import loadDataToNumpy_class
        ld_image = loadDataToNumpy_class.loadInDataToNumpy(loadPath + imFile)
        self.im = ld_image.loadData()
        self.imShape = np.shape(self.im)
        if refFile is not None:
            ld_image = loadDataToNumpy_class.loadInDataToNumpy(loadPath + refFile)
            self.ref = ld_image.loadData()
            self.refExists = True
        else:
            self.refExists = False
        if plotting:
            self.plot_raw_input()        

    def superGaussian(self, x, x0, w, power):
        """ Super Gaussian function
        """
        return np.exp(- (2 * (x- x0)/(w))**power)
    
    def createGaussianCroppingWindow(self, image, bot, top, left, right, power = 8):
        """ Create a window to crop in fourier space with
        """
        print ("Creating Gaussian To Crop. Image shape", np.shape(image))
        cropGaussian = np.zeros_like(image, dtype = float)
        s = np.shape(image)
        sgY = []
        sgX = []
        for y in range(s[1]):
            sgY.append( self.superGaussian(y, (left + right)*0.5, abs(right-left), power))
        for x in range(s[0]):
            sgX.append( self.superGaussian(x, (top + bot)*0.5, abs(top - bot), power))         

        for i, x in enumerate(sgX):
            for j, y in enumerate(sgY):
                cropGaussian[i][j] = x * y
    
        cropGaussian = np.real(cropGaussian)
        return cropGaussian            
            
    def zeroPadImages(self, padSize = 100):
        """ Pads an image with zeros
        """        
        self.padSize = padSize
        self.im_PlasmaChannel = np.pad(self.im_PlasmaChannel, padSize, 'constant')    
        if self.refExists:
            self.ref_PlasmaChannel = np.pad(self.ref_PlasmaChannel, padSize, 'constant')    
        
    def check_btlf_coors_in_image(self,bot, top, left, right, shape):
        """ Check the four coordinates that they are within the image
        Corrects them to be within if too large/small
        """
        bot, top, left, right = [int(bot), int(top), int(left), int(right)]
        # print (bot, top, left, right, shape)
        if bot < 0:
            bot = 0
        if left < 0:
            left = 0
        if top > shape[0]-1:
            top = shape[0]-1
        if right > shape[1]-1:
            right = shape[1]-1   
        return int(bot), int(top), int(left), int(right)
        
    
    def cropToPlasmaChannel(self, pc_crop_coors, plotting = False, 
                            paddingX = 30, paddingY = 10, 
                            padSize = 100,
                            verbose = False):
        """ Crop to the plasma channel, using the four corner coors
        The crop is done with a window function to remove the fourier efffect of
        sharp edges
        """     # Would like to add choice of windows
        bot, top, left, right = pc_crop_coors  
        bot, top, left, right = self.check_btlf_coors_in_image(bot, top, left, right, self.imShape)
        self.paddingX = paddingX
        self.paddingY = paddingY    
        self.padSize = padSize           
        if plotting:
            for y in [top, bot]:
                plt.hlines(y, left, right, color = "black")
            for x in [left, right]:
                plt.vlines(x, top, bot, color = "black" )  
            plt.imshow(self.im)
            plt.show()
        power = 2*6
        gaus_cropping = self.createGaussianCroppingWindow(self.im, bot, top , left, right , power)
        if False:
            plt.imshow(gaus_cropping)
            plt.colorbar()
            plt.title("Gaussian cropping window")
            plt.show()
        self.im_PlasmaChannel = self.im * gaus_cropping
        if self.refExists:
            self.ref_PlasmaChannel = self.ref * gaus_cropping
        
        
        b_pc, t_pc, l_pc, r_pc = bot - paddingY, top + paddingY, left - paddingX, right + paddingX
        b_pc, t_pc, l_pc, r_pc = self.check_btlf_coors_in_image(b_pc, t_pc, l_pc, r_pc, self.imShape)
        
        # Crop to the new region of the image
        self.im_PlasmaChannel = self.im_PlasmaChannel[b_pc:t_pc, l_pc:r_pc]
        if self.refExists:
            self.ref_PlasmaChannel = self.ref_PlasmaChannel[b_pc:t_pc, l_pc:r_pc]

        # Pad the image with zeros
        self.zeroPadImages(padSize)
        
        if plotting:
            plt.clf()
            print ("Plotting the raw images cropped with a window")
            plt.title("Raw image with Gauss window applied and padding")
            plt.pcolormesh(self.im_PlasmaChannel , cmap = plt.cm.seismic, norm = func.MidpointNormalize(midpoint=0))
            plt.axes().set_aspect('equal')
            plt.colorbar()            
            plt.show()       
            
    def fft_of_plasma(self, plotting = True):
    # Do 2D fft 
        self.F_im_PC = np.fft.fft2(self.im_PlasmaChannel)
        self.F_imShape = np.shape(self.F_im_PC)
        self.F_im_PC = self.wrapImage(self.F_im_PC)
        if self.refExists:
            self.F_ref_PC = np.fft.fft2(self.ref_PlasmaChannel)
            self.F_ref_PC = self.wrapImage(self.F_ref_PC)
                    
        if plotting:    
            ## # It is not helpful to see the output of the reference cropping
            # if self.refExists:
            #     f, ax = plt.subplots(ncols = 2)
            #     im1 = ax[0].pcolormesh( abs(self.F_im_PC), cmap = plt.cm.seismic, norm = mpl.colors.LogNorm())
            #     im2 = ax[1].pcolormesh( abs(self.F_ref_PC), cmap = plt.cm.seismic, norm = mpl.colors.LogNorm())
            #     plt.colorbar(im1, ax = ax[0])
            #     plt.colorbar(im2, ax = ax[1])                
            #     ax[0].set_title("Im PC FFT")
            #     ax[1].set_title("Ref PC FFT")                
            #     plt.show()
            # else:
            plt.pcolormesh( abs(self.F_im_PC), cmap = plt.cm.seismic, norm = mpl.colors.LogNorm())
            plt.colorbar()
            plt.title("Im PC FFT")
            plt.show()
            
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
            
    def moveCenterOfImage(self, image, moveX, moveY):
        wrappedImage = np.zeros_like(image)
        shape = np.shape(image)
        for i in range(shape[0]):
            for j in range(shape[1]):
                # Wrapping image
                i_w = (i + moveX)%shape[0]
                j_w = (j + moveY)%shape[1]
                # print 
                wrappedImage[i_w][j_w] = image[i][j]
        return wrappedImage
    
    def plot_FT_space(self, bot, top, left, right, peakRelHeight):
        f, ax = plt.subplots(nrows = 2, figsize = (8, 6), sharex = True, 
                             gridspec_kw={'height_ratios': [3,1]})
        ax[0].pcolormesh(abs(self.F_im_PC), cmap = plt.cm.seismic
                    , norm=mpl.colors.LogNorm()
                    )
        lenData = len(np.sum(abs(self.F_im_PC), axis = 0))
        xData = range(lenData)
        yData = np.sum(abs(self.F_im_PC), axis = 0)
        # xData = xData[lenData//3: 2 * lenData//3]
        # yData = yData[lenData//3: 2 * lenData//3]
        ax[1].plot(xData, yData)
        for y in [top, bot]:
            ax[0].hlines(y, left, right, color = "black")
        for x in [left, right]:
            ax[0].vlines(x, top, bot, color = "black" )  
        ax[1].set_xlim([xData[0], xData[-1]])

        from scipy.signal import find_peaks            
        peaks = find_peaks(np.sum(abs(self.F_im_PC), axis = 0) , 
                      # width = 20
                      height = np.sum(abs(self.F_im_PC), axis = 0).max() * peakRelHeight
                      )
        print ("Peaks in f space", peaks)
        ax[1].plot(peaks[0], peaks[1]['peak_heights'], "x")
        plt.suptitle("The crop region in F space")
        func.tightLayout_with_suptitle()
        plt.show()
        
    def plot_cropWindow_and_cropped_Fpeak(self, bot, top, left, right, gaussCrop):
        f, ax = plt.subplots(nrows = 2, sharex = True, gridspec_kw={'height_ratios': [1, 3]} )
        ax[0].set_title("Cropping window")            
        im_c = ax[0].pcolormesh(gaussCrop)
        for y in [top, bot]:
            ax[0].hlines(y, left, right, color = "white")
        for x in [left, right]:
            ax[0].vlines(x, top, bot, color = "white" )        

        ax[1].set_title("Cropped peak in Fourier space")
        im = ax[1].pcolormesh( abs(self.F_cropped_im)  , cmap = plt.cm.viridis
                , norm=mpl.colors.LogNorm(), vmin  = 1
                ) 
        plt.colorbar(im_c, ax = ax[0])                    
        plt.colorbar(im, ax = ax[1])        
        for y in [top, bot]:
            ax[1].hlines(y, left, right, color = "white")
        for x in [left, right]:
            ax[1].vlines(x, top, bot, color = "white" )
        plt.show()
            
    def crop_to_FFT_peak(self, crop_coors, GausPower = 6,  plot_crop_window_and_peak = True,
                         plot_fft_space = True, peakRelHeight = 0.3):
        
        bot, top, left, right = crop_coors
        F_shape = np.shape(self.F_im_PC)
        bot, top, left, right = self.check_btlf_coors_in_image(bot, top, left, right, F_shape)
        self.F_cropped_im = np.zeros_like(self.F_im_PC)
        
        gaussCrop = self.createGaussianCroppingWindow(self.F_im_PC, bot, top, left, right, GausPower)
        self.F_cropped_im = gaussCrop * self.F_im_PC
        if self.refExists:
            self.F_cropped_ref = gaussCrop * self.F_ref_PC        

        if plot_fft_space:
            self.plot_FT_space(bot, top, left, right, peakRelHeight)
            
        if plot_crop_window_and_peak:
            self.plot_cropWindow_and_cropped_Fpeak(bot, top, left, right, gaussCrop)
            
            
    def find_peaks_in_lineouts(self, peakRelHeight):
        """ Create a lineout of the image and find the peaks
        """
        from scipy.signal import find_peaks           
        x_range = range(len(np.sum(abs(self.F_im_PC), axis = 0)))
        xData = np.sum(abs(self.F_im_PC), axis = 0)
        y_range = range(len(np.sum(abs(self.F_im_PC), axis = 1)))
        yData = np.sum(abs(self.F_im_PC), axis = 1)
        
        # Locate the peaks
        xpeaks = find_peaks(xData , height = xData.max() * peakRelHeight)         
        ypeaks = find_peaks(yData , height = yData.max() * peakRelHeight)  
        # Take the first and last peak to be the peaks of interest. Tune the variable peakRelHeight
        # so this is the case
        if len(xpeaks[0]) > 2:
            xpeaks = [ [xpeaks[0][0], xpeaks[0][-1]], 
                      {'peak_heights': [xpeaks[1]['peak_heights'][0], xpeaks[1]['peak_heights'][-1]]} ]

        return x_range, xData, y_range, yData, xpeaks, ypeaks
        
            
    def auto_select_FT_peak(self, GausPower = 6, peakRelHeight = 0.6, yPercRange = 0.3, xPercRange = 0.3, 
                            plot_fft_space = True, plotting_cropping = True, plot_found_peaks = True):
        """ Auto select the FT crop region from the peak locations
        """
        F_shape = np.shape(self.F_im_PC)
        self.F_cropped_im = np.zeros_like(self.F_im_PC)

        # Find the peaks in the lineouts
        x_range, xData, y_range, yData, xpeaks, ypeaks = self.find_peaks_in_lineouts(peakRelHeight)        
        if plot_found_peaks:        
            plt.figure(figsize = (6, 4))
            plt.title("Finding the peaks in the FFT of the image")
            plt.plot(x_range, xData, label = "x")
            plt.plot(y_range, yData, label = "y")
            plt.plot(xpeaks[0], xpeaks[1]['peak_heights'], "x")
            plt.plot(ypeaks[0], ypeaks[1]['peak_heights'], "x")
            plt.legend()
            plt.show()        

        xrange = xpeaks[0][1] - xpeaks[0][0]
        yrange = F_shape[0]
        left = xpeaks[0][1] - xPercRange * xrange
        right = xpeaks[0][1] + xPercRange * xrange        

        # print (F_shape)
        bot = ypeaks[0][0] - yrange * yPercRange
        top = ypeaks[0][0] + yrange * yPercRange

        bot, top, left, right = self.check_btlf_coors_in_image(bot, top, left, right, F_shape)
        print ("Found crop coors", bot, top, left, right)
        
        gaussCrop = self.createGaussianCroppingWindow(self.F_im_PC, bot, top, left, right, GausPower)
        self.F_cropped_im = gaussCrop * self.F_im_PC
        if self.refExists:
            self.F_cropped_ref = gaussCrop * self.F_ref_PC

        if plot_fft_space:
            self.plot_FT_space(bot, top, left, right, peakRelHeight)
            
        if plotting_cropping:
            self.plot_cropWindow_and_cropped_Fpeak(bot, top, left, right, gaussCrop)
        
        # return the found crop coors
        return bot, top, left, right 
            
    
    def createPhase_inverseFT(self, plotting = False):
        self.result_im = np.fft.ifft2(self.F_cropped_im)
        
        if self.refExists:
            self.result_ref = np.fft.ifft2(self.F_cropped_ref)
            self.phaseShift = np.angle(self.result_im / self.result_ref)
        else:
            self.phaseShift = np.angle(self.result_im)
        
        self.phaseShift = self.phaseShift[self.padSize:-self.padSize 
                                          ,self.padSize:-self.padSize]
        self.phaseShift = self.phaseShift[self.paddingY:-self.paddingY,
                                          self.paddingX:-self.paddingX]
        if plotting:
            plt.title("Phase Shift Due to Plasma")
            plt.pcolormesh(self.phaseShift , cmap = plt.cm.seismic, vmin = -np.pi, vmax = np.pi)
            cbar = plt.colorbar(ticks=[-np.pi, 0, np.pi])
            cbar.ax.set_yticklabels([r'$-\pi$', '0', r'$\pi$'])
            plt.show()

    def move_CroppedRegion_To_Centre(self, plotting = True):
        F_lineout = np.average(abs(self.F_cropped_im), axis = 0)
        maxInd = func.nearposn(F_lineout, F_lineout.max())
        
        if plotting:
            f, ax = plt.subplots(nrows = 2, sharex = True, figsize = (6, 3), gridspec_kw={'height_ratios': [1, 3]})
            ax[0].set_title("Peak in fourier space before moving")
            ax[0].plot(F_lineout)
            ax[1].pcolormesh(abs(self.F_cropped_im), norm = mpl.colors.LogNorm(), vmin = 1)
            plt.show()        
            
        self.F_cropped_im = self.moveCenterOfImage(self.F_cropped_im, 0, self.F_imShape[1]//2 - maxInd)
        if self.refExists:
            self.F_cropped_ref = self.moveCenterOfImage(self.F_cropped_ref, 0, self.F_imShape[1]//2 - maxInd)            

        if plotting:
            F_lineout = np.average(abs(self.F_cropped_im), axis = 0)
            f, ax = plt.subplots(nrows = 2, sharex = True, figsize = (6, 3), gridspec_kw={'height_ratios': [1, 3]})
            ax[0].set_title("Peak in fourier space after moving")            
            ax[0].plot(F_lineout)
            ax[1].pcolormesh(abs(self.F_cropped_im), norm = mpl.colors.LogNorm() , vmin = 1)
            plt.show()       
    
        self.F_cropped_im = self.wrapImage(self.F_cropped_im)
        if self.refExists:
            self.F_cropped_ref = self.wrapImage(self.F_cropped_ref)


if __name__ == "__main__":
    loadPath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/"
    imFile = "im_bgRemoved_croppedROI.txt"
    refFile = "reference.txt"

    
    ps = phaseShift()
    ps.load_data(loadPath, imFile, refFile,  plotting = False)

    # Crop to the plasma channel
    pc_crop_coors = [120, 200, 30, 360] # bot, top, left, right
    ps.cropToPlasmaChannel( pc_crop_coors, plotting=False,
                            paddingX=20, paddingY=20, padSize = 50
                           )    
    ps.fft_of_plasma(plotting = True)
    
    FT_crop_coors = [50, 140, 250, 330] # bot, top, left, right

    # ps.crop_to_FFT_peak(FT_crop_coors)
    
    ps.auto_select_FT_peak(yPercRange = 0.3, xPercRange = 0.25, plot_fft_space = True, plotting_cropping = False, plot_found_peaks = False)
    ps.move_CroppedRegion_To_Centre(plotting = False) 
    
    ps.createPhase_inverseFT(plotting = True)
    
    # np.savetxt('/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/phaseShift_needsUnwrapping_2.txt',
    #            ps.phaseShift)
    
