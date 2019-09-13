#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Thu Jul 18 15:33:19 2019

@author: chrisunderwood

    Updating the shot data base with the parameters used.
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [8.0,8.0]
import matplotlib.pyplot as plt

# Load my module of functions
import CUnderwood_Functions3 as func

def readDatabase_from_csv(file):
    f = open(file, "r")
    data = f.readlines()
    f.close()
    d = []
    for l in data:
        d.append(l.split(','))
    return d

def search_db_shot(db, day, run, s):
    for row, shot in enumerate(db):
        if shot[0] == day:
            if shot[1] == run:
                if shot[2] == s:
                    print (row, shot)
                    return row
    
def update_db(db, index, string):
    updates = string.split(",")
    print (len(updates))
    print (len(db[index]))
    
    for i in range(len(updates)):
        db[index][i+6] = updates[i]

def save_db_csv(db, name):
    f = open(name, "w")
    for row in db:
        outString = ''
        for r in row:
            outString += "{},".format(r)
        f.write(outString[:-1])
    f.close()
    

databasePath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/2019RadiationReactionDatabase.csv"
db = readDatabase_from_csv(databasePath)
new_databasePath = "/Volumes/GoogleDrive/My Drive/Experimental_Codes/Interferometry_Data_Extraction/2019RadiationReactionDatabase_new.csv"

if __name__ == "__main__":
    
    StringToUpdate = ",31,133,224,440,0.9081420138919725,122,249,129,447,78, 20,8,123,4,311,78,36,173,236,8e-07, 4.22e-05\n"
    print ("String to add", StringToUpdate, "\n")
    
    pathToData = "/Volumes/GoogleDrive/Shared drives/Murphy Group/GeminiRR_2019/20190208/20190208r011/20190208r011s002_Probe_Interfero.tiff"
    data_ID = pathToData.split("/")[-1].split("_")[0]
    print (data_ID)
    
    shotIndex = search_db_shot(db, '20190212', '004', '011')
    
    print (shotIndex)
    update_db(db, shotIndex, StringToUpdate)
    
    save_db_csv(db, new_databasePath)

