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
        elif self.filePath.endswith(".txt"):
            print ("Loading .txt")           
            data = np.loadtxt(self.filePath)
        elif self.filePath.endswith(".tiff") or self.filePath.endswith(".TIFF"):
            print ("Loading .tiff/TIFF")
            from skimage import io
            data = io.imread(self.filePath)
            data = data.astype(float)
            # print (type(data)) 
        elif self.filePath.endswith(".tif"):
            print ("Loading .tif")            
            from skimage import io
            data = io.imread(self.filePath)
            data = data.astype(float)
            # print (type(data))            
        else:
            print ("The type is not an expected file")
            print ("File: ", self.filePath)
            print ("Please Edit: {}\n\tTo accept this file type".format("loadDataToNumpy_class.py"))
        return data

