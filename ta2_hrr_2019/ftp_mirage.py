import getpass
from ftplib import FTP_TLS
import os  
import numpy as np
import time
import sys


class MirageBoxSync():
    def __init__(self,user='m.streeter09@imperial.ac.uk',dataFolder = r'E:\Streeter2019\MIRAGE',boxFolder=r'TA2_HRR_2019/MIRAGE'):
        self.connect2box()
        self.dataFolder = dataFolder
        self.boxFolder = boxFolder
    def connect2box(self):
        """
        Connects to the box sync server 
        """

        ftp = FTP_TLS('ftp.box.com')  

        print('logging in')
        ftp.login(user=self.user, passwd = getpass.getpass())
        # move to destination directory
        ftp.cwd(self.boxFolder)  
        print('FTP: moving to folder ' + boxFolder +' \n')
        self.ftp  = ftp

    def disconnectBox(self):
        self.ftp.close()
        print('FTP: logged out')

    def copyFile(self,filePath):
        filename = os.path.split(filePath)[1]
        with open(filePath, 'rb') as fb:
            self.ftp.storbinary("STOR " + filename,fb)

    def ftpMakeDir(self,dirName):
        try:
            self.ftp.mkd(dirName)
        except:
            print('FTP: could not create ' + dirName + '. Maybe it exists?')

    def findDiags(self):
        self.diagList = []
        for root, dirs, files in os.walk(self.dataFolder, topdown=True):
            self.diagList = dirs

    def syncDate2box(self,dateStr = None):
        if dateStr is None:
            print('Enter a date string i.e.  20190904')
            return None

        self.findDiags()

        for diag in self.diagList:
            ftp.cwd(self.boxFolder)  
            print('FTP: making diag folder ' + diag + ' \n')
            ftpMakeDir(diag)
            print('FTP: moving to diag folder ' + diag + ' \n')
            ftp.cwd(os.path.join(self.boxFolder,diag))
            print('FTP: making date folder ' + dateStr + ' \n')
            ftpMakeDir(dateStr)
            print('FTP: moving to date folder ' + dateStr + ' \n')
            destPath = os.path.join(self.boxFolder,diag,dateStr)
            ftp.cwd(destPath)
            newPath = destPath
            sourcePath = os.path.join(self.dataFolder,diag,dateStr)
            
            print('FTP: walking through date folder \n')
            # walk through directories copying files
            for root, dirs, files in os.walk(sourcePath, topdown=True):
                relPath = os.path.relpath(root, sourcePath)
                if relPath is not '.':
                    newPath = os.path.join(destPath,relPath)
                    ftp.cwd(newPath)
                    print('FTP: moving to folder ' + newPath + ' \n')

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
                        print('   Done  \n')

                for dirName in dirs:
                    if dirName in existingDirs:
                        print('FTP: remote directory exists' + dirName + ' \n')
                    else:   
                        print('FTP: making remote directory ' + dirName + ' \n')
                        self.ftpMakeDir(dirName)

                print('All files copied' + ' \n')




            





