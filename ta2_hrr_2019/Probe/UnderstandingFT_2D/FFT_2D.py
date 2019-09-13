#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""    _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Thu Jun 20 10:10:36 2019

@author: chrisunderwood
"""
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = [8.0,5.0]
import matplotlib.pyplot as plt

# Load my module of functions
import sys
sys.path.insert(0, '/Users/chrisunderwood/Documents/Python/')
import CUnderwood_Functions3 as func

def wrapImage(image):
    wrappedImage = np.zeros_like(image)
    shape = np.shape(image)
    for i in range(shape[0]):
        for j in range(shape[1]):
            # Wrapping image
            i_w = (i + shape[0]//2)%shape[0]
            j_w = (j + shape[1]//2)%shape[1]
            # print 
            wrappedImage[i_w][j_w] = image[i][j]
    return wrappedImage

def returnHalfFFTimage(fimage):
    newImage = np.zeros_like(fimage)
    shape = np.shape(fimage)
    for i in range(shape[0]):
        for j in range(shape[1]//2):
            newImage[i][j] = fimage[i][j]
    return newImage
    
def superGaussian(x, x0, w, power):
    return np.exp(- (2 * (x- x0)/(w))**power)

def cropSuperGaussian(image, bot, top, left, right, power = 8):
    cropGaussian = np.zeros_like(image)
    s = np.shape(image)
    for x in range(s[0]):
        for y in range(s[1]):
            sgX = superGaussian(x, (left + right)*0.5, abs(right-left), power)
            sgY = superGaussian(y, (top + bot)*0.5, abs(top - bot), power)
            cropGaussian[x][y] = sgX * sgY
    cropGaussian = np.real(cropGaussian)
    return image * cropGaussian, cropGaussian
    
def createSinWaves_atAngle():
    image = [] 
    for x in xArr:
        output = []
        for y in yArr:
            value = np.sin(x + 5 * y) 
            output.append( value)
        image.append(output)
    return image

def gaussAndSin():
    image = [] 
    for x in xArr:
        output = []
        for y in yArr:
        # The plasma channel will not make the signal brighter, just shift it
        # if anything it will be dimmer, but thats too much to worry about at the 
        # moment.
            plasmaStart = 1
            plasmaEffect = np.heaviside(x+plasmaStart,0) * 1.
            plasmaEffect *= (np.sin(4*(x+3.5 +  x*0.05)) -  backgroundSinWaves(x)) * np.exp(- (y / 2.)**6 )
            value = plasmaEffect + backgroundSinWaves(x) 
            output.append( value)
        image.append(output)
    image = np.array(image)
    return image.T

def backgroundSinWaves(x):
    return np.sin(4*(x-(1.5)))

def referenceToPlasma():
    image = [] 
    for x in xArr:
        output = []
        for y in yArr:
            value = backgroundSinWaves(x)
            output.append( value)
        image.append(output)
    image = np.array(image)
    return image.T


def plotAbsWithCbar(image, title = ""):
    plt.title(title)
    plt.imshow( abs(image), cmap = plt.cm.gnuplot)
    plt.colorbar()
    plt.show()       

def wrapped_FFT(im):
    fim = np.fft.fft2(im)
    F_wrapped = wrapImage(fim)
    return F_wrapped

def calculatePhase(im):
    phaseImage = np.arctan(np.imag(im)/ np.real(im))
    return phaseImage

xArr = np.linspace(-5, 5, 100)
yArr = np.linspace(-5, 5, 100)

image = gaussAndSin()
refIm = referenceToPlasma()
# image = createSinWaves_atAngle()

print ("Two images to look at, reference and mock plasma channel")
plotAbsWithCbar(image, "Image to play with")
plotAbsWithCbar(refIm, "refIm to play with")

print ("Images in fourier space")
F_image = wrapped_FFT(image)
plotAbsWithCbar(F_image, "FFT")

F_refIm = wrapped_FFT(refIm)
plotAbsWithCbar(F_refIm, "Ref FFT")

print ("Recontruction of one peak, hard cut off")
# Hard Cut off half the box
croppedInFSpace = returnHalfFFTimage(F_image)
plotAbsWithCbar(croppedInFSpace, "Cropped F Space - Hard cut")
imageCropped = np.fft.ifft2(croppedInFSpace)
plotAbsWithCbar(imageCropped, "Reconstructed Image Cropped in FT space - Phase with wrapping" )

# # Looking at the reference image.
# crop_FSpace_refIm = returnHalfFFTimage(F_refIm)
# plotAbsWithCbar(crop_FSpace_refIm, "Ref Cropped F Space")
# refCropped = np.fft.ifft2(crop_FSpace_refIm)
# plotAbsWithCbar(refCropped, "Ref Reconstruction Cropped in FT space" )

# =============================================================================
# This is the method that we want to be using 
# =============================================================================
print ("2D Gauss filter --  Hamming window")
# 2D Gauss filter
croppedInFSpace, gauss = cropSuperGaussian(F_image, 20, 50, 40, 60, power= 6)
croppedInFSpace_ref, gauss = cropSuperGaussian(F_refIm, 20, 50, 40, 60, power= 6)

plotAbsWithCbar(gauss, "Cropped gauss")
plotAbsWithCbar(croppedInFSpace, "Cropped F Space with gauss")
imageCropped = np.fft.ifft2(croppedInFSpace)
plotAbsWithCbar(imageCropped, "Reconstructed Image Cropped in FT space" )

plotAbsWithCbar(croppedInFSpace_ref, "Cropped F Space with gauss")
imageCropped_ref = np.fft.ifft2(croppedInFSpace_ref)
plotAbsWithCbar(imageCropped_ref, "Reconstructed Image Cropped in FT space" )

print ("The re calculated image and reference, \nand then calculated the angle between them for the phase")
showFourierReconstructions = False
phaseImage = calculatePhase(imageCropped)
if showFourierReconstructions: plotAbsWithCbar(phaseImage, "Phase Image")

phaseRef = calculatePhase(imageCropped_ref)
if showFourierReconstructions: plotAbsWithCbar(phaseRef, "Phase Ref")

plotAbsWithCbar(phaseImage - phaseRef, "Phase difference")

phaseShift_2 = np.angle(imageCropped / imageCropped_ref)
plotAbsWithCbar(phaseShift_2, "Phase by angle between two images")



