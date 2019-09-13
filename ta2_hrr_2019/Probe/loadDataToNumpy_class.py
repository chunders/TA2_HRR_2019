#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Wed Jun 19 09:32:01 2019

@author: chrisunderwood

    Load data class
    turns the data into a numpy array
"""
import numpy as np

class loadInDataToNumpy():
    def __init__(self, filePath):
        """ Initialise with the filepath
        """
        self.filePath = filePath

    def loadData(self):
        """ Loads the data
        This looks at the file extension and loads it appropiately
        """
        if self.filePath.endswith(".png"):
            print ("png files not yet supported")
            return None
        
        if self.filePath.endswith(".txt"):
            	data = np.loadtxt(self.filePath)
        if self.filePath.endswith(".tiff"):
            	from skimage import io
            	data = io.imread(self.filePath)
            	data = data.astype(float)
            	# print (type(data))
        return data

