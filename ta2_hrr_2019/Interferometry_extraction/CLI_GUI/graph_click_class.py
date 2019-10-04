#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Thu Jul 11 16:32:30 2019

@author: chrisunderwood
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

class graphClick():
    def __init__(self, fig,  ax, number_of_coors = 2):
        self.fig = fig
        self.ax = ax
        self.coors = []
        self.number_of_coors = number_of_coors
        cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        plt.show()
    
    def onclick(self, event):
        global ix, iy
        ix, iy = event.xdata, event.ydata
        print ('x = %d, y = %d'%(ix, iy) )
    
        global coords
        coords = [int(ix), int(iy)] #Turn into pixels, as these are images
    
        self.coors.append(coords)
        
        if len(self.coors) > self.number_of_coors - 1:         
            plt.close()
            return self.correctly_return_coors()

    def plot_square_from_clicks(self, ax):
        ax.vlines(self.coors[0][0], self.coors[0][1], self.coors[1][1])
        ax.vlines(self.coors[1][0], self.coors[0][1], self.coors[1][1])  
        ax.hlines(self.coors[0][1], self.coors[0][0], self.coors[1][0])
        ax.hlines(self.coors[1][1], self.coors[0][0], self.coors[1][0]) 
        plt.show() 

    def plot_plasmaChannel_from_clicks(self, ax, image):
    	shape = np.shape(image)
    	ax.hlines(self.coors[0][1], 0, shape[1])
    	ax.hlines(self.coors[1][1], 0, shape[1])
    	plt.show()

    	return [abs(self.coors[0][1] + self.coors[1][1]) // 2, abs(self.coors[0][1] - self.coors[1][1])]

    def correctly_return_coors(self):
        t, b, l, r = [self.coors[0][1], self.coors[1][1], self.coors[0][0], self.coors[1][0]]
        if t > b:
            t1 = t
            t = b
            b = t1
        if l > r:
            l1 = l
            l = r
            r = l1   
        return t, b, l, r  





