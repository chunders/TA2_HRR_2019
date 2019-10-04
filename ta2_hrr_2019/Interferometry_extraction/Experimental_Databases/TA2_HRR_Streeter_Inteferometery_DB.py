#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Wed Oct  2 15:53:59 2019

@author: chrisunderwood
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [8.0,8.0]
import matplotlib.pyplot as plt
import os

# Load my module of functions
import sys
sys.path.insert(0, '/Users/chrisunderwood/Documents/Python/')
import CUnderwood_Functions3 as func

Header = "day,run,burst,shot,ref,plasmaChannel,,,,paddingX,paddingY,padSize,fftCropRegion,,,,rotationAngle,mPerPix,mask_threshold_level\n"

def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs



class probeDatabase():
    def __init__(self, fileName):
        
        self.fileName = fileName
        self.appendingToDB = False
        if os.path.exists(self.fileName):
            self.f=open(self.fileName, "r")
            self.appendingToDB = True
            self.lines = self.f.readlines()
            self.lines = self.lines[1:] # Remove the header
            self.f.close()
            
    def create_string_to_write(self, day, run, burst, shot, ref, plasmaChannel,
            paddingX, paddingY, padSize, fftCropRegion,rotationAngle, mPerPix, mask_threshold_level):
        """ Create a list of all the variables that need to be inputted for a 
        repeatable extraction.
        """
        output_csv = ""
        for item in [day, run, burst, shot, ref, plasmaChannel, paddingX, 
                     paddingY, padSize, fftCropRegion, rotationAngle, mPerPix, mask_threshold_level]:
            if type(item) in [int, str, float, np.float64]:
                output_csv += str(item) + ","
            elif item == None:
                output_csv += ","
            elif len(item) > 1:
                for i in item:
                    output_csv += str(i) + ","
        self.output_csv = output_csv[:-1] + "\n"    
        print (self.output_csv)
        self.day, self.run, self.burst, self.shot = [int(day), int(run.split("n")[1]),
                  int(burst.split("t")[1]), int(shot.split("t")[1].split(".")[0]) ]
        
        
    def save_line_to_file(self):
        """ Add the extraction data for this shot to the file
        """
        self.f = open(self.fileName, "w")
        if self.appendingToDB:
            if len(self.lines) < 1:
                self.appendingToDB = False
                
        if self.appendingToDB:
            print ("Appending")
            self.check_lines_for_duplicates()
            self.lines.append(self.output_csv)
            self.sortLines()
            self.f.write(Header)
            for l in self.lines:
                self.f.write(l)
            self.f.close()
        else:
            print ("New")            
            self.lines = [self.output_csv]
            self.f.write(Header)
            for l in self.lines:
                self.f.write(l)
            self.f.close()        

    def check_lines_for_duplicates(self):
        indexToKeep = []
        for i, l in enumerate(self.lines):
            day = int(l.split(",")[0])
            run = int(l.split(",")[1].split("n")[1])
            burst = int(l.split(",")[2].split("t")[1])
            shot = int(l.split(",")[3].split("t")[1].split(".")[0])                

            if self.day == day and self.run == run and self.burst == burst and self.shot == shot:
                print ("Duplicate", i)
            else:
                indexToKeep.append(i)
        self.lines = np.array(self.lines)[indexToKeep].tolist()
                
            


    def sortLines(self):
        """ Sort the lines so they are in the order of
        - Day
            - Run 
                - Burst
                    - Shot
        Take into account the header
        """
        print ("Sorting the lines")
        print (self.lines)
        unsortedList = self.lines
        sortingArr = [] # Dont sort the header
        for l in unsortedList: # Dont sort the header
            print (l)
            print ("Options", l.split(",")[0], l.split(",")[1].split("n")[1], l.split(",")[2].split("t")[1], l.split(",")[3].split("t")[1].split(".")[0])

            day = int(l.split(",")[0])
            run = int(l.split(",")[1].split("n")[1])
            burst = int(l.split(",")[2].split("t")[1])
            shot = int(l.split(",")[3].split("t")[1].split(".")[0])    
            print (day, run, burst, shot)            
            sortingArr.append(day * 1e8 + run *1e6 + burst*1e3 + shot)
    
        l = self.output_csv 
        day = int(l.split(",")[0])
        run = int(l.split(",")[1].split("n")[1])
        burst = int(l.split(",")[2].split("t")[1])
        shot = int(l.split(",")[3].split("t")[1].split(".")[0])    
        new_item = day * 1e8 + run *1e6 + burst*1e3 + shot
            
        print (sortingArr)
        print (list_duplicates_of(sortingArr, new_item))
        
        sortedList, sortingArr = func.sortArrAbyB(unsortedList, sortingArr)
        self.lines = []
        for l in sortedList:
            self.lines.append(l)
            
            

if __name__ == "__main__":            
    savePath = "/Volumes/GoogleDrive/My Drive/2019_Streeter_TA2/Probe_Interferometry/"
    fileName = "TA2_HRR_Probe_Database.csv"
    
    day = "20190113"
    run = "run016"
    burst = "Burst7"
    shot = "Shot35.TIFF"
    ref = "/Shot50.TIFF"
    plasmaChannel = [350,500,1100,1200]
    paddingX = 200
    paddingY = 10
    padSize = 100
    fftCropRegion = [120,260,360,420]
    rotationAngle = 2.657329950002826
    mPerPix = 4.98107e-06
    
    pdb = probeDatabase(savePath + fileName)
    pdb.create_string_to_write(day, run, burst, shot, ref, plasmaChannel, paddingX, paddingY, padSize, fftCropRegion,
           rotationAngle, mPerPix)
    
    pdb.save_line_to_file()
    
