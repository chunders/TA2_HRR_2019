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
    
    Background Remover
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Load my module of functions
import CUnderwood_Functions3 as func

import unittest
import backgroundRemover


class test_backgroundRemover(unittest.TestCase):
    
    def test_RemovingBackgroundReducesAverage(self):
        LoadPath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/verticalFrings.txt"
            # Call the class
        bc = backgroundRemover.backgroundRemoverAndCrop()
        bc.LoadData(LoadPath)
        bc.blur(bc.im, 11)
        bot = 40
        top = 300
        left = 100
        right = 500
        out = bc.sub_blurredIm(bot, top, left, right)    
        orig = bc.im[bot:top, left:right]
        
        self.assertTrue(abs(np.average(out)) < abs(np.average(orig)))

    def test_blur(self):
      # I can't think how to test this as it is just a scipy function....
      pass
      

        

if __name__ == "__main__":
	unittest.main()
