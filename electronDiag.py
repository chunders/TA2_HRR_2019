import numpy as np
import cv2
import scipy.signal as sig
from scipy.interpolate import interp2d,interp1d
import scipy.io as sio

## Should make an object to store the calibration info and interpolation functions

class screenCalibration(calFilePath):
    __init__

def specCalFromAngle(x,bAng,calFilePath):
    calData = sio.loadmat(calFilePath)
    specCalFun = interp2d(calData['xScreen'], calData['exitAngle'], calData['electronEnergy'], kind='linear', copy=True, bounds_error=False, fill_value=0)
    return specCalFun(x,bAng)

def img2spec(eSpecScreen,bAng,eAxis):
    """ Converts electron spectromter screen data to the energy axis
    arguments: 
        eSpecScreen: a dictionary containing 'x_mm' and 'img' 
        bAng: the angle of the beam relative to the center in the plane of dispersion
        eAxis: the energy axis which the output will be interpolated onto
    """
    x = (eSpecScreen['x_mm'])*1e-3
    img = eSpecScreen['img']
    eEnergy = specCalFromAngle(x,bAng)
    imgSize = np.shape(img)
    spec = img/np.abs(np.tile(np.gradient(eEnergy),(imgSize[0],1)))
    specFun = interp1d(eEnergy,spec, bounds_error=False, fill_value=None)
    specLin = specFun(eAxis)
    return specLin

def img2screen(img,tForm):
    """ De-warps the raw image to give the data on linear spatial axes
    Arguments:
         img: the raw image data
         tForm: dictionary containing transform information of the form
         tForm = {  
            'description': description
            'H': result of cv2.findHomography
            'newImgSize':number of pixels in resultant image (Nx_new,Ny_new)
            'x_mm': output spatial axes x in mm,
            'y_mm': output spatial axes y in mm,
            'imgArea0': pixel area in original image,
            'imgArea1': pixel area in new image
            }
    """

    imgArea0 = tForm['imgArea0']
    H = tForm['H']
    newImgSize = tForm['newImgSize']
    imgArea1 = tForm['imgArea1']
    img = img - np.median(img)
    img = sig.medfilt2d(img,kernel_size=(3,3))
    imgCountsPerArea = img/imgArea0
    imgCountsPerArea[imgArea0==0] =0
    imgCountsPerArea[np.isinf(imgCountsPerArea)] = 0
    imgCountsPerArea[np.isnan(imgCountsPerArea)] = 0

    imgScreen = cv2.warpPerspective(imgCountsPerArea, H, newImgSize)*imgArea1
    return imgScreen


def calcTransform(imgP_pixels,imgP_real,imgSize,newSize,newLims):
    """ Calculate transoform properties using cv2.findHomography
    x is horizontal (1st axis) and y is vertical (0th axis)
    arguments:
        imgP_pixels np array of [y,x] pixel coordinates of control points for transform
        imgP_real np array of [y,x] spatial coordinates of control points for transform
        imgSize (Ny,Nx) size of orginal image in pixels
        newSize (Ny,Nx) size of resultant image
        newLims [yMin,xMin;yMax,xMax] screen limits in real units
    """

    # physical range covered by output image
    yRange = newLims[1,0]-newLims[0,0]
    xRange = newLims[1,1]-newLims[0,1]

    # physical resolution of output image
    dy = yRange/newSize[0]
    dx = xRange/newSize[1]
    y_mm = newLims[0,0] + np.linspace(0,newSize[0],num=newSize[0])*dy
    x_mm = newLims[0,1] + np.linspace(0,newSize[1],num=newSize[1])*dx

    # Pixel values in new image
    imgP_trans = (imgP_real-[newLims[0,0],newLims[0,1]]) / [dx,dy]

    H, status = cv2.findHomography(imgP_pixels,imgP_trans)

    (Ny,Nx) = imgSize

    # calculate pixel areas in original image
    retval,H_inv = cv2.invert(H)
    (X,Y) = np.meshgrid(x_mm,y_mm)
    X_raw = cv2.warpPerspective(X, H_inv, (Nx,Ny))
    Y_raw = cv2.warpPerspective(Y, H_inv, (Nx,Ny))
    imgArea0 = np.gradient(X_raw,axis=1)*np.gradient(Y_raw,axis=0)

    tForm = {  
        'H': H,
        'newImgSize':newSize,
        'x_mm': x_mm,
        'y_mm': y_mm,
        'imgArea0': imgArea0,
        'imgArea1': dx*dy
    }

    return tForm