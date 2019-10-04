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

    Testing recreate Phase
    
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Load my module of functions
import CUnderwood_Functions3 as func

import unittest
import createPhase_class

class test_recreatePhase(unittest.TestCase):

	def test_inverseFFT(self):
		pass

	def test_crop_to_FFT_peak(self):

		loadPath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/"
		imFile = "im_bgRemoved_croppedROI.txt"
		refFile = "reference.txt"
		ps = createPhase_class.phaseShift()
		ps.load_data(loadPath, imFile, refFile, plotting=False)
		bot = 140
		top = 180
		left = 30
		right = 360
		ps.cropToPlasmaChannel( bot, top, left, right, False)

		ps.fft_of_plasma(False)		
		ps.crop_to_FFT_peak(10, 30, 250, 300, False)

		self.assertEqual(np.shape(ps.F_cropped_im), np.shape(ps.F_cropped_ref))

	def test_wrapImage(self):
	# Test that the sum of the image remains the same.
		ps = createPhase_class.phaseShift()
		image = np.ones((100, 30))

		for i in range( 10, 20):
			for j in range(10, 20):
				image[i][j] = i * j * 0.2

		wrapped = ps.wrapImage(image)
		self.assertAlmostEqual(np.sum(wrapped), np.sum(image))


	def test_createCropSuperGaussian(self):
		ps = createPhase_class.phaseShift()
		image = np.zeros((1000, 1000))
		bot, top, left, right  = [10, 50, 200, 500]
		cropGaussian = ps.createCropSuperGaussian(image, bot, top, left, right)
		self.assertAlmostEqual(cropGaussian.max(), 1)
		self.assertAlmostEqual(cropGaussian.min(), 0)

if __name__ == "__main__":
	unittest.main()