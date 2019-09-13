import getpass
from ftplib import FTP_TLS
import os  
import numpy as np
import time
import sys
import re


class MirageBoxSync():
    def __init__(self,user='m.streeter09@imperial.ac.uk',dataFolder = r'E:\Streeter2019\MIRAGE',boxFolder=r'/TA2_HRR_2019/MIRAGE'):
        self.user = user
        self.dataFolder = dataFolder
        self.boxFolder = boxFolder
        self.connect2box()
    def connect2box(self):
        """
        Connects to the box sync server 
        """

        ftp = FTP_TLS('ftp.box.com')  

        print('logging in')
        ftp.login(user=self.user, passwd = getpass.getpass())
        # move to destination directory
        ftp.cwd(self.boxFolder)  
        print('FTP: moving to folder ' + self.boxFolder +' \n')
        self.ftp  = ftp

    def disconnectBox(self):
        self.ftp.close()
        print('FTP: logged out')

    def copyFile(self,filePath):
        filename = os.path.split(filePath)[1]
        if filename not in self.ftp.nlst():
            with open(filePath, 'rb') as fb:
                try:
                    self.ftp.storbinary("STOR " + filename,fb)
                except error_perm:
                    connect2box()
                    self.ftp.storbinary("STOR " + filename,fb)
            print('   Done  \n')
        else:
            print(filename + ' already exists, skipping')

    def ftpMakeDir(self,dirName):
        try:
            self.ftp.mkd(dirName)
        except:
            print('FTP: could not create ' + dirName + '. Maybe it exists?')

    def findDiags(self):
        self.diagList = []
        diagList = next(os.walk(self.dataFolder))[1]
        for d in diagList:
            if d[0] is not '.':                
                self.diagList.append(d)

    def syncDate2box(self,dateStr = None):
        ftp = self.ftp
        if dateStr is None:
            print('Enter a date string i.e.  20190904')
            return None
        print('Finding diagnostics...')
        self.findDiags()
        print(self.diagList)
        for diag in self.diagList:
            print('FTP: moving to folder ' + self.boxFolder +' \n')
            ftp.cwd(self.boxFolder)  
            print('FTP: making diag folder ' + diag + ' \n')
            self.ftpMakeDir(diag)
            diagFolder = self.boxFolder + '/' + diag
            print('FTP: moving to diag folder ' + diagFolder + ' \n')
            ftp.cwd(diagFolder)
            print('FTP: making date folder ' + dateStr + ' \n')
            self.ftpMakeDir(dateStr)
            destPath = diagFolder + '/' + dateStr
            print('FTP: moving to date folder ' + destPath + ' \n')
            ftp.cwd(destPath)
            newPath = destPath
            sourcePath = os.path.join(self.dataFolder,diag,dateStr)
            
            print('FTP: walking through date folder \n')
            # walk through directories copying files
            for root, dirs, files in os.walk(sourcePath, topdown=True):
                relPath = os.path.relpath(root, sourcePath).replace('\\','/')
                if relPath is not '.':
                    newPath = destPath + '/' + relPath
                    print('relPath: ' + relPath)
                    print('FTP: moving to folder ' + newPath + ' \n')
                    ftp.cwd(newPath)
                    

                existingFiles=[]
                existingDirs = []
                for name, facts in ftp.mlsd(newPath,["type"]):
                    if facts["type"] == "file":
                        existingFiles.append(name)
                    elif facts["type"] == "dir":
                        existingDirs.append(name)
                for fileName in files:
                    if fileName in existingFiles:
                        print('FTP: skipping existing file ' + fileName + ' \n')
                    else:
                        print('FTP: writing file ' + fileName)
                        self.copyFile(os.path.join(root,fileName))
                        

                for dirName in dirs:
                    if dirName in existingDirs:
                        print('FTP: remote directory exists' + dirName + ' \n')
                    else:   
                        print('FTP: making remote directory ' + dirName + ' \n')
                        self.ftpMakeDir(dirName)

            print('All files copied' + ' \n')




            





