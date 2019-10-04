#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Tue Jul  9 19:32:35 2019

@author: chrisunderwood
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [6.0,4.0]
import matplotlib.pyplot as plt


import sys
sys.path.insert(0, '/Users/chrisunderwood/Documents/Experimental_Tools/Interferometry_extraction/')
# Load my module of functions
# import CUnderwood_Functions3 as func
from Interferometry_Extraction import Interferometry
# import recreatePhase
import graph_click_class
import time
from loadDataToNumpy import loadInDataToNumpy
from matplotlib.colors import LogNorm


Run_Selecting_Regions = True

class interferometery_interface(Interferometry):

    def Select_Region_With_Fringes(self):
        print ("Stage one: Raw image rotation for vertical fringes\n\tClick around region with no plasma channel to locate angle of fringes")

        self.rawMin = self.im_raw.min()
        self.rawMax = self.im_raw.max()        
        if Run_Selecting_Regions: 
            time.sleep(0.1)

            # Have a loop around the input, incase it is done incorrectly.
            inputFlag = True
            while inputFlag:
                fig, ax = plt.subplots(ncols = 1, figsize = (6,4))                
                self.plotRawImage(ax)
                g = graph_click_class.graphClick(fig, ax)
                #Draw box on image around clicks
                fig, ax = plt.subplots(ncols = 1, figsize = (6,4))
                self.plotRawImage(ax)
                g.plot_square_from_clicks(ax)
                plt.show() 

                userInput = int(input("Choose Again?\n\tYes = 1\n\tNo = 0\n"))
                if userInput == 0:
                    inputFlag = False
                else:
                    print ("Try selecting region again")
                    inputFlag = True
            print ("Coors: ", g.coors)
            # Format is: bot, top, left, right
            self.fringeCrop = g.correctly_return_coors()
        else:
            # Test set of coors
            self.fringeCrop = [29, 181, 270, 320]
        print (self.fringeCrop)      

    def plot_fringeCrop(self):
        plt.figure(figsize=(6,4))
        plt.imshow(self.im_raw[self.fringeCrop[0]:self.fringeCrop[1],
                   self.fringeCrop[2]: self.fringeCrop[3]])
        plt.show()


    def Select_cropROI(self, givenRegion = False, plotComparison = True):
        print ("Stage two: Crop rotated image to region of interest\n\tClick around region of interest")
        if Run_Selecting_Regions and not givenRegion:
            time.sleep(0.1)
            
            if hasattr(self, "vertFringesIm"):
                imageToSelectROI = self.vertFringesIm
            else:
                imageToSelectROI = self.im_raw
                self.rawMin = imageToSelectROI.min()
                self.rawMax = imageToSelectROI.max()
                        
            # Have a loop around the input, incase it is done incorrectly.
            inputFlag = True
            while inputFlag:
                fig, ax = plt.subplots(ncols = 1, figsize = (6,4)) 
                im = ax.imshow(imageToSelectROI, vmin = self.rawMin, vmax = self.rawMax )
                plt.colorbar(im, ax = ax)
                g = graph_click_class.graphClick(fig, ax)
                #Draw box on image around clicks
                fig, ax = plt.subplots(ncols = 1, figsize = (6,4))
                im = ax.imshow(imageToSelectROI, vmin = self.rawMin, vmax = self.rawMax )
                g.plot_square_from_clicks(ax)
                plt.show() 

                userInput = int(input("Choose Again?\n\tYes = 1\n\tNo = 0\n"))
                if userInput == 0:
                    inputFlag = False

                    # make sure the results are in the correct order.
                    t, b, l, r = [g.coors[0][1], g.coors[1][1], g.coors[0][0], g.coors[1][0]]
                    if t > b:
                        t1 = t
                        t = b
                        b = t1
                    if l > r:
                        l1 = l
                        l = r
                        r = l1
                    self.crop_ROI = [t, b, l, r]
                else:
                    print ("Try selecting region again")
                    inputFlag = True
            print ("Coors: ", self.crop_ROI )
        else:
            self.crop_ROI = givenRegion
        
        self.removeBackground_CropROI(ROI_btlr = self.crop_ROI, plotting = False)

        if plotComparison and hasattr(self, "ref"):
            f, ax = plt.subplots(ncols = 2)
            ax[0].imshow(self.im_raw)
            ax[1].imshow(self.ref)
            plt.show()

    def create_referenceImage(self):
        print ("Stage three: \n\tClick above the plasma channel and below the channel")
        if Run_Selecting_Regions: 
            time.sleep(0.1)   
            # Have a loop around the input, incase it is done incorrectly.
            inputFlag = True
            while inputFlag:
                fig, ax = plt.subplots(ncols = 1, figsize = (6,4))
                ax.imshow(self.im_bg_removed)
                g = graph_click_class.graphClick(fig, ax, 2)
                fig, ax = plt.subplots(ncols = 1, figsize = (6,4))
                ax.imshow(self.im_bg_removed)        
                center, width = g.plot_plasmaChannel_from_clicks(ax = ax, image = self.im_bg_removed)
                userInput = int(input("Choose Again?\n\tYes = 1\n\tNo = 0\n"))
                if userInput == 0:
                    inputFlag = False
                    self.center_height_PC = center
                    self.height_width_around_PC = width
                else:
                    print ("Try selecting region again")
                    inputFlag = True
        else:
            self.center_height_PC = 78
            self.height_width_around_PC = 20
        self.createReference(line_of_pchannel = self.center_height_PC,
                             width_pchannel = self.height_width_around_PC)

        # Produced: self.ref
    def load_reference_image(self, fileName):
        ld = loadInDataToNumpy(fileName)
        self.ref_file_path = fileName
        self.ref = ld.loadData()


    def createPhase_with_fft(self, plot_inputs = True):
        print ("Stage three: \n\tSelect the plasma channel with two clicks")

        if Run_Selecting_Regions: 
            time.sleep(0.1)   
            # Have a loop around the input, incase it is done incorrectly.
            inputFlag = True
            while inputFlag:
                fig, ax = plt.subplots(ncols = 1, figsize = (6,4))
                plt.title("Choose the plasma channel")
                ax.imshow(self.im_bg_removed)
                g = graph_click_class.graphClick(fig, ax, 2)    

                fig, ax = plt.subplots(ncols = 1, figsize = (6,4))
                im = ax.imshow(self.im_bg_removed)
                g.plot_square_from_clicks(ax)
                plt.show() 

                userInput = int(input("Choose Again?\n\tYes = 1\n\tNo = 0\n"))
                if userInput == 0:
                    inputFlag = False
                    self.crop_To_PlasmaChannel = g.correctly_return_coors()
                else:
                    print ("Try selecting region again")
                    inputFlag = True
        else:
            self.crop_To_PlasmaChannel = [8, 123, 4, 311] 

        print ("Cropping to plasma channel")
        print (self.crop_To_PlasmaChannel)
        paddingX = 400
        paddingY = 100
        padSize = 100
        self.crop_plasma_channel(self.crop_To_PlasmaChannel, plotting = plot_inputs,
                    paddingX = paddingX, paddingY = paddingY, padSize = padSize)

        print (self.ps.paddingX, self.ps.paddingY)


        if Run_Selecting_Regions:
            # self.ps = recreatePhase.phaseShift()
            self.ps.load_arrIntoClass(self.im_PlasmaChannel, self.ref_PlasmaChannel)
            
            self.ps.im_PlasmaChannel = self.im_PlasmaChannel
            self.ps.ref_PlasmaChannel = self.ref_PlasmaChannel

            self.F_im_PC = np.fft.fft2(self.im_PlasmaChannel)
            self.F_imShape = np.shape(self.F_im_PC)
            self.F_im_PC = self.ps.wrapImage(self.F_im_PC)

            # Have a loop around the input, incase it is done incorrectly.
            inputFlag = True
            while inputFlag:
                fig, ax = plt.subplots(ncols = 1, figsize = (12,6))
                plt.title("Select the FFT region")
                ax.pcolormesh( abs(self.F_im_PC), norm = LogNorm(), cmap = plt.cm.jet)
                g = graph_click_class.graphClick(fig, ax)  

                # Replot the figures with the selected box
                fig, ax = plt.subplots(ncols = 1, figsize = (12,6))
                im = ax.pcolormesh( abs(self.F_im_PC), norm = LogNorm(), cmap = plt.cm.jet)
                g.plot_square_from_clicks(ax)
                plt.show() 

                userInput = int(input("Choose Again?\n\tYes = 1\n\tNo = 0\n"))
                if userInput == 0:
                    inputFlag = False
                    self.fftCropRegion = [g.coors[0][1], g.coors[1][1], g.coors[0][0], g.coors[1][0]]
                else:
                    print ("Try selecting region again")
                    inputFlag = True    
        else:
            self.fftCropRegion = [78, 36, 173, 236]   
        self.createPhase(self.fftCropRegion, plotCropFFTPeak = False,
                    plot_fft_plasma = False,
                    plot_moveFFTCentre = False, zeroPad = None, 
                    plot_fft_space = False, plot_final_result = False)

    def unwrap_phase_image(self):
        print ("Stage four: Unwrap the phase.")
        f, ax = plt.subplots(ncols = 2)
        ax[0].pcolormesh(self.phase)
        phase_lineout = np.sum(self.phase, axis = 1)
        ax[1].plot(phase_lineout, range(len(phase_lineout)))
        plt.show()        
        
        pc_bot =  int(input("Bottom of plasma Channel (height)\n"))
        pc_top = int(input("Top of plasma Channel (height)\n"))
        background = int(input("Top where I can get background from (height)\n"))
        self.unwrapPhase_both_directions(pc_bot, pc_top, background)


    def abel_transform_phase(self, mPerPix,  plot_n_e_result = True, laserwavelength_m = 800e-9, rotationAngle = 1.9):
        self.mPerPix = mPerPix
        self.laserwavelength_m = laserwavelength_m
        self.recreateElectronDensity(plot_n_e_result = plot_n_e_result, laserwavelength_m  = self.laserwavelength_m,
                            mPerPix = self.mPerPix, 
                            rotationAngle = rotationAngle)





if __name__ == "__main__":

    from path_to_data import *
   
    gui = interferometery_interface(filePath)

    badInput = True
    while badInput:
        userInput_refernce = int(input("Is there a reference image?\n\tYes = 1\n\tNo = 0\n"))
        if userInput_refernce == 0:
            build_reference = True
            badInput = False
        elif userInput_refernce == 1:
            build_reference = False
            badInput = False            


    if build_reference:
        gui.Select_Region_With_Fringes()
        # gui.plot_fringeCrop()
        gui.rotateImage_VertFringes(cropRegion = gui.fringeCrop, plotting = False)    
        if not Run_Selecting_Regions:
            gui.Select_cropROI(givenRegion = [153, 233, 138, 461], plotComparison = False)
        else:
            gui.Select_cropROI(plotComparison = False)
        gui.create_referenceImage()
    else:
        gui.load_reference_image(userInput_ref_filePath)
        gui.Select_cropROI(plotComparison = False)

    gui.createPhase_with_fft(plot_inputs = False)
    gui.unwrap_phase_image()

    mPerPix = input("Enter the size that each of the pixels corresponds to in [m]:\n")
    rotationAngle = float(input("Enter the angle in degrees:\n"))
    gui.abel_transform_phase(mPerPix= float(mPerPix), rotationAngle=rotationAngle)


    saveString = ""
    print ("\n\nThings to save:")
    if build_reference:
        saveString += ","   # This would be the file Name
        # Things to do with creating the reference
        print("\tBUILDING OWN REFERENCE:\n\tFringe Crop: {}\n\tFringe Angle {}\n\tcrop_ROI: {}".format(
                        gui.fringeCrop, 
                        gui.fringe_rotationAngle,
                        gui.crop_ROI)
                )
        for num in gui.fringeCrop:
            saveString += "{},".format(num) 
        saveString += "{},".format(gui.fringe_rotationAngle) 
        for num in gui.crop_ROI:
            saveString += "{},".format(num)         
    else:
        print ("\tLoading_reference {}".format(gui.ref_file_path))

        saveString += "{},".format(gui.ref_file_path)
        saveString += ",,,,,,,,," # these would be the refence creation steps

    if build_reference:
        # Things to recreate the phase.
        print ("\n\tPlasma channel {} and width {}\n\tPlasma Channel Crop {} \n\tFFT Crop {}".format(gui.center_height_PC, gui.height_width_around_PC,
                gui.crop_To_PlasmaChannel, gui.fftCropRegion)
            # Abel transform
            + "\n\tLaser wavelength {} \n\tSize of Pixels (m) {}".format(
                        gui.laserwavelength_m, gui.mPerPix)

            )
        saveString += "{}, {},".format(gui.center_height_PC, gui.height_width_around_PC)
        for num in gui.crop_To_PlasmaChannel:
            saveString += "{},".format(num)
        for num in gui.fftCropRegion:
            saveString += "{},".format(num)

    saveString += "{}, {}".format(gui.laserwavelength_m, gui.mPerPix)

    print (saveString)
    for i, l in enumerate(saveString.split(",")):
        print (i, l)



    

    

    
    # analysis.rotateImage_VertFringes() # If angle is given it just rotates this angle



    
