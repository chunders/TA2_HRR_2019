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

    Testing rotate array
    
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Load my module of functions
import CUnderwood_Functions3 as func

import rotateArray
import loadDataToNumpy_class


import unittest

class test_rotateArray(unittest.TestCase):
    
    def test_findRotation(self):
        filePath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/Test_Data/"
        fileName = filePath + "SyntheticReference.tiff"
        testData = loadDataToNumpy_class.loadInDataToNumpy(fileName)
        d = testData.loadData()

        rot = rotateArray.rotateFringesToVertical(d)
        rotation = rot.findRotation(startIndex = 0, numberOfLineouts = 15, 
                                    horzStart = 100, horzEnd = 500)
        print ("Angle in Degrees", rotation)
        print ("Expected Angle ", 0)
        self.assertAlmostEqual(rotation, 0, places = 5)
        
        
    

if __name__ == "__main__":
	unittest.main()