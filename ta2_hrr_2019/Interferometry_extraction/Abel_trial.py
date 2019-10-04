#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Fri Jul  5 16:30:35 2019

@author: chrisunderwood
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [6.0,4.0]
import matplotlib.pyplot as plt
import abel
# Load my module of functions
import CUnderwood_Functions3 as func

xsize = 20
ysize = 30
testArray = np.zeros((ysize, xsize))

x = np.linspace(-5, 5, xsize)
y = np.linspace(-15, 15, ysize)
plt.pcolormesh(x,y, testArray)
plt.colorbar()
plt.show()


for i, yp in enumerate(y):
    for j, xp in enumerate(x):
        value = np.exp( - ((yp)/3.)**2) * abs(xp + 5)
        testArray[i][j] = value

testArray = (testArray / testArray.max()) * np.pi
plt.pcolormesh(x,y, testArray)
plt.colorbar()
plt.show()

testArray = testArray.T
inverse_abel = abel.Transform(testArray, #center = (50, 200),
                        center = "gaussian",
                        center_options = {'axes' : 1, "verbose":True},
                      direction='inverse',
                      verbose = True).transform.T

plt.pcolormesh(x,y,inverse_abel)
plt.colorbar()
plt.show()
