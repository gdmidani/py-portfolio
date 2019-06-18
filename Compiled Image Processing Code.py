# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 17:28:40 2019

@author: gaith
"""

"""COMPILED IMAGE PROCESSING CODE"""

# -*- coding: utf-8 -*-

# By Gaith Midani

"""Astrophotographic Filter Response Project"""

import numpy as np
from scipy import *
import matplotlib.pyplot as plt
import pylab as plb


lambd = range(0,119)
LambdPeak = range(0,119)
#LambdLow = range(0,119)
#LambdHigh = range(0,119) 

fLambd = 60	# Arbitrary constant selected for source flux

##def R():
##        DataArray = np.genfromtxt('IRAS60mu.dat')
##        return DataArray

class R:
        def __call__(self):
                DataArray = np.genfromtxt(['IRAS60mu.dat'])
                return DataArray

R = np.genfromtxt('IRAS60mu.dat')
Rx = np.genfromtxt('IRAS60mu.dat', usecols = 0)
Ry = np.genfromtxt('IRAS60mu.dat', usecols = 1)

Rx_max = np.nanmax(Rx)
Ry_max = np.nanmax(Ry)
Rmax = np.nanmax(R)

# Peak wavelength λpeak's role in response function:
Ry_max = R([LambdPeak])

# Low wavelength λlow and high wavelength λhigh's roles in response function
# Double-check the math from other sources... obviously LambdLow and LambdHigh are different:
(1/2)*Ry_max = LambdLow
(1/2)*Ry_max = LambdHigh

FWHM = LambdHigh - LambdLow

# Using Simpson's Rule to integrate:
def integrate(Ry, h):
    i=1
    total=y_vals[0]+Ry[-1]
    for y in Ry[1:-1]:
        if i%2 == 0:
            total+=2*y
        else:
            total+=4*y
        i+=1
    return total*(h/3.0)

# Ry as the stand-in for the general response function due to the requirements of integral definition above
Bandwidth = (1/Rmax)*integrate(Ry) 

LambdCenter = 0.5*(Bandwidth)

LambdMean = (Lambd*integrate(Ry))/(fLambd*integrate(Ry))

LambdIso = (1/Rmax)*(fLambd*integrate(Ry))

# Color correction will adjust the effective wavelength to its proper place. Even without knowing the source
# shape, we can use a single spectral shape for all detected sources and report measurements at nominal
# survey wavelengths (in this case, 60.0 microns). For IRAS, the approximation is:
def S(Lambd):
	Src = (1/Lambd)
	return Lambd

"""Astrophotographic Standardization"""

##An observer used B and V filters to obtain four exposures of a standard field at different air masses: two B exposures at air masses 1.05 and 2.13, and two V exposures at airmasses 1.10 and 2.48. In this field, four photometric standard stars were included. Their measured magnitudes are given in the following table.
##  	
##		(B-V)	V	b(1)	b(2)	v(1)	v(2)
##Airmass		        1.05    2.13    1.10    2.48
##
##Star A	-0.07  12.01   9.853   10.687   8.778   9.427
##
##Star B	0.36   12.44   10.693  11.479   9.160   9.739
##
##Star C	0.69   12.19   10.759  11.462   8.873   9.425
##
##Star D	1.15   12.89   11.898  12.547   9.522   10.001
##
##    Calculate extinction coefficients for this instrumental system at B and V bands.
##    Compute the standard transformation coefficients αV and αB-V (or αB)
##    Calculate standard magnitudes of Obj1 (i.e., V and B-V) whose instrumental magnitudes are v=9.850 and b=10.899 taken at airmass=1.50
##
# Plot (b(2)-b(1))/(X2-X1) and find slope for k1(B-V), as appears in this equation:
# mAtm = m + k(Lambd)*X ~ m + (k0+k1(B-V))*X

import numpy as np
import matplotlib.pyplot as plt

Ab21 = 10.687-9.853
Bb21 = 11.479-10.693
Cb21 = 11.462-10.759
Db21 = 12.547-11.898

Xb21 = 2.13-1.05

Av21 = 9.427-8.778
Bv21 = 9.739-9.160
Cv21 = 9.425-8.873
Dv21 = 10.001-9.552

Xv21 = 2.48-1.10

# Here's the setup to calculate values of k1(B-V) for each star in B-band:
AbSetup = Ab21/Xb21
BbSetup = Bb21/Xb21
CbSetup = Cb21/Xb21
DbSetup = Db21/Xb21

# ...And for the V-Band:
AvExt = Av21/Xv21
BvExt = Bv21/Xv21
CvExt = Cv21/Xv21
DvExt = Dv21/Xv21

# Class definition for future reference:
class Line(object):

    def __init__(self, data):
            self.first, self.second = data

    def slope(self):
            (x1, y1), (x2, y2) = self.first, self.second
            try:
                    return (float(y2)-y1)/(float(x2)-x1)
            except ZeroDivisionError:
                    # line is vertical
                    return None

    def yintercept(self, slope):
            if slope != None:
                    x, y = self.first
                    return y - self.slope * x
            else:
                    return None

# Let's consider the current unknowns, m and k0*X, in terms of knowns:
# m + k0*X = (b1 or b2 or v1 or v2) - Ext*X

# For Star A, where b1/b2/v1/v2 represent the atmospheric magnitudes:
AUnknownSumb1 = 9.853 - AbExt*1.05
AUnknownSumb2 = 10.687 - AbExt*2.13
AUnknownSumv1 = 8.778 - AvExt*1.10
AUnknownSumv2 = 9.427 - AvExt*2.48

# For Star B:
BUnknownSumb1 = 10.693 - BbExt*1.05
BUnknownSumb2 = 11.479 - BbExt*2.13
BUnknownSumv1 = 9.160 - BvExt*1.10
BUnknownSumv2 = 9.739 - BvExt*2.48

# For Star C:
CUnknownSumb1 = 10.759 - CbExt*1.05
CUnknownSumb2 = 11.462 - CbExt*2.13
CUnknownSumv1 = 8.873 - CvExt*1.10
CUnknownSumv2 = 9.425 - CvExt*2.48

# For Star D: 
DUnknownSumb1 = 11.898 - DbExt*1.05
DUnknownSumb2 = 12.547 - DbExt*2.13
DUnknownSumv1 = 9.522 - DvExt*1.10
DUnknownSumv2 = 10.001 - DvExt*2.48

# Which can be used to find the extinction coefficients

# Color indices B-V for the stars:
ABV = -0.07
BBV = 0.36
CBV = 0.69
DBV = 1.15

# Stand-in for mSTD - m as a function of color index:
funcNumer = mSTD - m
funcA = (funcNumer, ABV) 
funcB = (funcNumer, BBV)
funcC = (funcNumer, CBV)
funcD = (funcNumer, DBV)

Alpha1A = slope(funcA)
Alpha10A = yintercept(funcA)

Alpha1B = slope(funcB)
Alpha10B = slope(funcB)

...

# Standard transformation coefficients:
# For stars at b1
AmSTDb1 = 9.853 + Alpha1A*(ABV) + Alpha10A
BmSTDb1 = 10.693 + Alpha1B*(BBV) + Alpha10B
CmSTDb1 = 10.759 + Alpha1C*(CBV) + Alpha10C
DmSTDb1 = 11.898 + Alpha1C*CBV + Alpha10D

...

"""Astrophotographic Bad Pixel Correction"""

# Mini-Project 1:
# By Gaith Midani
# ASTR 3010

# Result: successfully processed images, but ds9 display has not been properly utilized

import numpy as np
import pyfits
from ds9 import *

def fixbads(data,mask,fwhm=2.5,boxsize=11,silent=True,fixneg=True):
  """
  ; Fix bad pixels (from a mask) using Gaussian weighted interpolation.
  ;
  ; Bad pixels are : 
  ;    - pixels with mask value of zero
  ;    - pixels with pixel value of zero or negative
  ; 
  ; Interpolation scheme :
  ;    using a Gaussian kernel window of size '2*boxsize+1' with FWHM
  ;
  ; Badpixel mask :
  ;    0 = good pixels
  ;    1 = bad pixels
  ;
  ; fixneg = False --> Do not attemp to fix negative flux values
  ;
  ; Created by Inseok Song, Sep. 2010
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  """
  from numpy import arange,exp,outer,zeros,real,where,ones
  
  if boxsize > 11.0:
     if not silent: print 'Gaussian kernel size is too large. Reduced to 11.0 pixels!'
     boxsize=11.0

  if fixneg:
     mask[ where(data <= 0) ] = 1
 
  # create a Gaussian kernel of size [-boxsize <--> boxsize] with FWHM
  x=arange(2*boxsize+1)
  y=exp(-0.5*(x-boxsize)**2 / (fwhm/2.354)**2)
  kernel = outer(y,y)

  # Preparing a bigger array with padding on all four sides of the data array.
  #   This step is necessary to correct bad pixels that are too close
  #   to edges. 
  # mask2 is a padded mask array where values are -1 in the padded area.
  BIG = zeros( tuple( real(data.shape) + 2*boxsize ) )
  mask2 = -ones( tuple( real(data.shape) + 2*boxsize ) )

  BIG[boxsize:-boxsize,boxsize:-boxsize] = data.copy()
  mask2[boxsize:-boxsize, boxsize:-boxsize] = mask.copy()
  
  badx, bady = where (mask2 > 0)    # x,y positions of badpixels

  # perform a bad pixel interpolation using the Gaussian weight kernel
  badpos = zip(badx,bady)
  for (ix,iy) in badpos:
      sub_data =  BIG[ix-boxsize:ix+boxsize+1, iy-boxsize:iy+boxsize+1]
      sub_mask = (mask2[ix-boxsize:ix+boxsize+1, iy-boxsize:iy+boxsize+1]).copy()
      bads = where( sub_mask != 0 )
      sub_mask = 1.0 + sub_mask*0
      sub_mask[bads] = 0.0   # so that bad pixels are not used in interpolation
      weight = kernel * sub_mask
      if weight.sum() != 0:
         weight /= weight.sum()
         new_value = (weight*sub_data).sum()
         BIG[ix,iy] = new_value

  return BIG[boxsize:-boxsize,boxsize:-boxsize].copy()

def rdfitsfiles(filelist,silent=False,**kwargs):
   """Read multiple FITS files from a filelist into a 3D array.
   """
   import pyfits,numpy
 
   if not silent: print 'Reading ',filelist[0]
   imdata,header = pyfits.getdata(filelist[0],header=True,**kwargs)
 
   try: nfiles = len(filelist)
   except: return imdata
 
   for file in filelist[1:]:
       if not silent: print 'Reading ',file
       im = pyfits.getdata(file,**kwargs)
       imdata = numpy.dstack((imdata,im))
 
   return imdata
  

# Collection of biases define a master bias:

bias_files = ['D46.fits',\
'D47.fits',\
'D48.fits',\
'D49.fits',\
'D50.fits',\
'D51.fits',\
'D52.fits']
biases = rdfitsfiles(bias_files)
master_bias = np.median(biases,axis=2) 


# Collection of flats for each band, to subtract from master bias:

# B-band
Bflat_files = ['D01.fits',\
'D02.fits',\
'D03.fits',\
'D04.fits',\
'D05.fits']
flatsB = rdfitsfiles(Bflat_files) 
Bflat = np.median(flatsB,axis=2) 
Bflat -= master_bias
# V-band
Vflat_files = ['D06.fits',\
'D07.fits',\
'D08.fits',\
'D09.fits',\
'D10.fits']
flatsV = rdfitsfiles(Vflat_files) 
Vflat = np.median(flatsV,axis=2) 
Vflat -= master_bias
# I-band
Iflat_files = ['D11.fits',\
'D12.fits',\
'D13.fits',\
'D14.fits',\
'D15.fits']
flatsI = rdfitsfiles(Iflat_files) 
Iflat = np.median(flatsI,axis=2) 
Iflat -= master_bias

# Creating master flats by division of median:
master_BFlat = Bflat/np.median(Bflat)
master_VFlat = Vflat/np.median(Vflat)
master_IFlat = Iflat/np.median(Iflat)

# Sky files: 
sky_files=['D22.fits',\
'D23.fits',\
'D24.fits',\
'D25.fits',\
'D26.fits',\
'D27.fits']
skys=rdfitsfiles(sky_files)

# Ratioed sky files: last image / first image
ratio=skys[:,:,-1] / skys[:,:,0]
r_avg=np.median(ratio)
r_std=np.std(ratio)
BPM=np.zeros( ratio.shape )
# Bad pixels at +/- 4 standard deviations:
bads=np.where( (ratio < r_avg - 4.0*r_std) | (ratio > r_avg + 4.0*r_std) )
BPM[bads]=1

# Science target frames:
tgtB = pyfits.getdata('D43.fits')
tgtV = pyfits.getdata('D44.fits')
tgtI = pyfits.getdata('D45.fits')

# Bias subtraction:
tgtB -= master_bias
tgtV -= master_bias
tgtI -= master_bias
# Flat division:
tgtB /= master_BFlat
tgtV /= master_VFlat
tgtI /= master_IFlat
# Bad pixel correction:
final_tgt_B = fixbads(tgtB,BPM)
final_tgt_V = fixbads(tgtV,BPM)
final_tgt_I = fixbads(tgtI,BPM)

# Attempt at displaying image via DS9:

hduPrimary = pyfits.PrimaryHDU()
hduB = pyfits.ImageHDU(data=final_tgt_B, header = None)
hduV = pyfits.ImageHDU(data=final_tgt_V, header = None)
hduI = pyfits.ImageHDU(data=final_tgt_I, header = None)
hduAll = pyfits.HDUList(hdus=[hduPrimary,hduB,hduV,hduI])
d=ds9()
d.set_pyfits(hduAll)


# FileList.txt:

##D01.fits[986,995][real][flat][B_639]:DOME
##D02.fits[986,995][real][flat][B_639]:DOME
##D03.fits[986,995][real][flat][B_639]:DOME
##D04.fits[986,995][real][flat][B_639]:DOME
##D05.fits[986,995][real][flat][B_639]:DOME
##D06.fits[986,995][real][flat][V_641]:DOME
##D07.fits[986,995][real][flat][V_641]:DOME
##D08.fits[986,995][real][flat][V_641]:DOME
##D09.fits[986,995][real][flat][V_641]:DOME
##D10.fits[986,995][real][flat][V_641]:DOME
##D11.fits[986,995][real][flat][i_705]:DOME
##D12.fits[986,995][real][flat][i_705]:DOME
##D13.fits[986,995][real][flat][i_705]:DOME
##D14.fits[986,995][real][flat][i_705]:DOME
##D15.fits[986,995][real][flat][i_705]:DOME
##D16.fits[986,995][real][unknown][B_639]:SKY,FLAT
##D17.fits[986,995][real][unknown][B_639]:SKY,FLAT
##D18.fits[986,995][real][unknown][B_639]:SKY,FLAT
##D19.fits[986,995][real][unknown][B_639]:SKY,FLAT
##D20.fits[986,995][real][unknown][B_639]:SKY,FLAT
##D21.fits[986,995][real][unknown][B_639]:SKY,FLAT
##D22.fits[986,995][real][unknown][V_641]:SKY,FLAT
##D23.fits[986,995][real][unknown][V_641]:SKY,FLAT
##D24.fits[986,995][real][unknown][V_641]:SKY,FLAT
##D25.fits[986,995][real][unknown][V_641]:SKY,FLAT
##D26.fits[986,995][real][unknown][V_641]:SKY,FLAT
##D27.fits[986,995][real][unknown][V_641]:SKY,FLAT
##D28.fits[986,995][real][unknown][i_705]:SKY,FLAT
##D29.fits[986,995][real][unknown][i_705]:SKY,FLAT
##D30.fits[986,995][real][unknown][i_705]:SKY,FLAT
##D31.fits[986,995][real][unknown][i_705]:SKY,FLAT
##D32.fits[986,995][real][unknown][i_705]:SKY,FLAT
##D33.fits[986,995][real][unknown][i_705]:SKY,FLAT
##D34.fits[986,995][real][object][B_639]:PG1323-086
##D35.fits[986,995][real][object][V_641]:PG1323-086
##D36.fits[986,995][real][object][i_705]:PG1323-086
##D37.fits[986,995][real][object][B_639]:PG1525-071
##D38.fits[986,995][real][object][V_641]:PG1525-071
##D39.fits[986,995][real][object][i_705]:PG1525-071
##D40.fits[986,995][real][object][B_639]:PG1633+099
##D41.fits[986,995][real][object][V_641]:PG1633+099
##D42.fits[986,995][real][object][i_705]:PG1633+099
##D43.fits[986,995][real][object][B_639]:Target
##D44.fits[986,995][real][object][V_641]:Target
##D45.fits[986,995][real][object][i_705]:Target
##D46.fits[986,995][real][zero][i_705]:BIAS
##D47.fits[986,995][real][zero][i_705]:BIAS
##D48.fits[986,995][real][zero][i_705]:BIAS
##D49.fits[986,995][real][zero][i_705]:BIAS
##D50.fits[986,995][real][zero][i_705]:BIAS
##D51.fits[986,995][real][zero][i_705]:BIAS
##D52.fits[986,995][real][zero][i_705]:BIAS
##

"""Centroiding"""

#ASTR 3010
#Centroiding homework
#By Gaith Midani

import numpy as np
from scipy import ndimage

# Center of Mass calculation:
# Iterate the CM/centroid with different S/N 
# S/N: (1.0:20.0:0.1). Also, show the result of centroiding
# accuracy [offset = true position - measured position] with a plot
# of centroiding accuracy as a function of S/N (x-axis).

y = 0

noise=np.random.random((32,32))

# Signal components
x0=np.random.random()*32
y0=np.random.random()*32

# Calculating S/N 
SigAvg = (x0*y0)/2
StanDevBG = np.std(noise)
SN = SigAvg/StanDevBG

#Values intended to affect offset (incomplete):
SNnew1 = 1.0
SNnew2 = 20.0
SNnew3 = 0.1

def gauss2d(size, fwhm, center):
    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

src=gauss2d(32,5.0,[x0,y0])

fake_image = noise + src*20
#plt.imshow(noise,vmin=0,vmax=20,interpolation='nearest')
#plt.imshow(fake_image,vmin=0,vmax=20,interpolation='nearest')
#plt.imshow(fake_image) #Try plt.ndimage() to bypass axis inversion from imshow

CM = ndimage.measurements.center_of_mass(fake_image)

#
# Offset/Centroiding intended to be defined here
#

plt.plot(SN, offset)


"""Near-IR Plots"""

#
# Purpose: To plot NIR-mid IR (range of 1-28 microns) atmospheric transmissions
# 	   for airmass AM = 1.0 and precipitable water vapor PWV = 1.0mm and 5.0mm
#
# Source: Lord, S. D., 1992, NASA Technical Memorandum 103957
# Credit to Gemini Observatory
# 

from pylab import *
from numpy import *

#NIR transmission at am = 1 and pwv = 1 mm
NIR11trans = NIR11wl = loadtxt('/home/gaith/Dropbox/ASTR 3010 (Astro Observ)/NIRmktrans_airmass1_pwv1.dat',unpack=True,usecols=[0])

#NIR wavelength at am = 1 and pwv = 1 mm
NIR11wl = loadtxt('/home/gaith/Dropbox/ASTR 3010 (Astro Observ)/NIRmktrans_airmass1_pwv1.dat',unpack=True,usecols=[1])

#NIR transmission at am = 1 and pwv = 5 mm
NIR15trans = loadtxt('/home/gaith/Dropbox/ASTR 3010 (Astro Observ)/NIRmktrans_airmass1_pwv5.dat',unpack=True,usecols=[0])

#NIR wavelength at am = 1 and pwv = 5 mm
NIR15wl = loadtxt('/home/gaith/Dropbox/ASTR 3010 (Astro Observ)/NIRmktrans_airmass1_pwv5.dat',unpack=True,usecols=[1])

#Mid-IR transmission at am = 1 and pwv = 1 mm
mIR11trans = loadtxt('/home/gaith/Dropbox/ASTR 3010 (Astro Observ)/mIRmktrans_airmass1_pwv1.dat',unpack=True,usecols=[0])

#Mid-IR wavelength at am = 1 and pwv = 1 mm
mIR11wl = loadtxt('/home/gaith/Dropbox/ASTR 3010 (Astro Observ)/mIRmktrans_airmass1_pwv1.dat',unpack=True,usecols=[1])

#Mid-IR transmission at am = 1 and pwv = 5 mm
mIR15trans = loadtxt('/home/gaith/Dropbox/ASTR 3010 (Astro Observ)/mIRmktrans_airmass1_pwv5.dat',unpack=True,usecols=[0])

#Mid-IR wavelength at am = 1 and pwv = 5 mm
mIR15wl = loadtxt('/home/gaith/Dropbox/ASTR 3010 (Astro Observ)/mIRmktrans_airmass1_pwv5.dat',unpack=True,usecols=[1])

plt.plot(NIR11trans, NIR11wl)
plt.plot(NIR15trans, NIR15wl)
plt.plot(mIR11trans, mIR11wl)
plt.plot(mIR15trans, mIR15wl)

plt.legend(['NIR at AM = PWV = 1', 'NIR at AM = 1 & PWV = 5', 'mid-IR at AM = PWV = 1', 'mid-IR at AM = 1 & PWV = 5'], loc='upper right')

xlabel('Wavelength (microns)') 
ylabel('Transmission')
show()


"""
Using the filter response data of the IRAS 60mu filter (IRAS60mu.dat in astr3010 directory), write a python script to calculate following parameters.
 

    Bandwidth (W0), λLOW, λHIGH
    λcenter, λmean,  λpeak, λeff, and λiso
    Calculate the color-correction value for the case of a source spectrum is F(λ)~constant. The nominal wavelength of the IRAS 60mu filter is 60.0micron and remember that IRAS assumed the source spectrum to be λF(λ)~constant.

In the calculations of λeff & λiso and #3, assume a constant spectrum (F(λ)~constant) for the source spectrum.

------
"""

# Set some ground rules:
R() = response function
# as indicated by IRAS60mu.dat 

# Peak wavelength's role in response function:
R(λpeak) = Rmax

# Low and high wavelengths' role in response function:
R(λlow) = R(λhigh) = Rmax/2

# Relation of low and high wavelengths to FWHM:
FWHM = λhigh - λlow

λcenter = (λl + λhigh)/2
# and/or
λcenter = 0.5*W0
# where
W0 = (1/Rmax)*(antiderivative of R(λ)dλ)

λmean = (Antiderivative of λ*R(λ)dλ)/(antiderivative of R(λ)dλ)

# This ONLY applies IF the function is completely SYMMETRICAL. If not, disregard this:
λpeak = λcenter = λmean

# The effective and isophotal wavelengths below depend on the source spectrum. We'll define it as an
# arbitrary constant: fλ ~ constant

# The effective wavelength is a weighted mean wavelength (weighted by the source flux fλ):
λeff = (antiderivative of λ*fλ*R(λ)dλ) / (antiderivative of fλ*R(λ)dλ)

# Isophotal wavelength is the one for which we have:
λiso = W0*fλiso = (1/Rmax)*(antiderivative of fλ*R(λ)dλ)
 

# Color correction will adjust the effective wavelength to its proper place. Even without knowing the source
# shape, we can use a single spectral shape for all detected sources and report measurements at nominal
# survey wavelengths (in this case, 60.0 microns). For IRAS:
S(λ) ~ (1/λ)

"""Astrophotographic Image Shifting"""

# -*- coding: utf-8 -*-
##In ../Project_Data/Shift_Combine_Data directory of the course server (chaos), there are five FITS files obtained from a single box5 dither observation. Write a python script to do the following steps.
##
##    create a median sky from five images
##    subtract the median sky from each five images
##    shift the sky-subtracted images and median combine them. In shifting images, you need to find the source centroid position using the cntrd routine and the routine for sub-pixel image shifting can be done by using mskpy.image.imshift
##
##* To keep all pixels from all dither images (i.e., larger median image as shown in the lecture note example), you need to "pad" the input image before you use the mskpy.image.imshift routine. Check the lecture note for padding example.
##
##You need to submit a python script and the name of the final FITS image stored in your server directory.

from numpy import *
from ds9 import *
import scipy as np
import pyfits
import cntrd
import mskpy

im1=pyfits.getdata('N01.fits',ignore_missing_end=True)
im2=pyfits.getdata('N02.fits',ignore_missing_end=True)
im3=pyfits.getdata('N03.fits',ignore_missing_end=True)
im4=pyfits.getdata('N04.fits',ignore_missing_end=True)
im5=pyfits.getdata('N05.fits',ignore_missing_end=True)

# Defining and combining median values of our FITS files
median1 = np.median(im1)
median2 = np.median(im2)
median3 = np.median(im3)
median4 = np.median(im4)
median5 = np.median(im5)
MedianSky = median1 + median2 + median3 + median4 + median5

# Median subtraction
Corr1 = im1 - MedianSky
Corr2 = im2 - MedianSky
Corr3 = im3 - MedianSky
Corr4 = im4 - MedianSky
Corr5 = im5 - MedianSky

# Peak values of x- and y-components of our corrected value for centroiding
cntrdx1 = np.amax(Corr1, axis=0)
cntrdy1 = np.amax(Corr1, axis=1)
cntrdx2 = np.amax(Corr2, axis=0)
cntrdy2 = np.amax(Corr2, axis=1)
cntrdx3 = np.amax(Corr3, axis=0)
cntrdy3 = np.amax(Corr3, axis=1)
cntrdx4 = np.amax(Corr4, axis=0)
cntrdy4 = np.amax(Corr4, axis=1)
cntrdx5 = np.amax(Corr5, axis=0)
cntrdy5 = np.amax(Corr5, axis=1)

# FWHM for centroiding
def FWHM(X,Y):
	half_max = max(Y) / 2
	d = sign(half_max - array(Y[0:-1])) - sign(half_max - array(Y[1:]))
	left_idx = find(d > 0)[0]
	right_idx = find(d < 0)[-1]
	return X[right_idx] - X[left_idx]

FWHM1 = FWHM(cntrdx1,cntrdy1)
FWHM2 = FWHM(cntrdx2,cntrdy2)
FWHM3 = FWHM(cntrdx3,cntrdy3)
FWHM4 = FWHM(cntrdx4,cntrdy4)
FWHM5 = FWHM(cntrdx5,cntrdy5)

# Putting it together...
ShiftLoc1 = cntrd.cntrd(Corr1,cntrdx1,cntrdy1,FWHM1)
ShiftLoc2 = cntrd.cntrd(Corr2,cntrdx2,cntrdy2,FWHM2)
ShiftLoc3 = cntrd.cntrd(Corr3,cntrdx3,cntrdy3,FWHM3)
ShiftLoc4 = cntrd.cntrd(Corr4,cntrdx4,cntrdy4,FWHM4)
ShiftLoc5 = cntrd.cntrd(Corr5,cntrdx5,cntrdy5,FWHM5)

# Drizzling
NewShiftLoc1 = mskpy.image.imshift(ShiftLoc1,cntrdx1,cntrdy1)
NewShiftLoc2 = mskpy.image.imshift(ShiftLoc2,cntrdx2,cntrdy2)
NewShiftLoc3 = mskpy.image.imshift(ShiftLoc3,cntrdx3,cntrdy3)
NewShiftLoc4 = mskpy.image.imshift(ShiftLoc4,cntrdx4,cntrdy4)
NewShiftLoc5 = mskpy.image.imshift(ShiftLoc5,cntrdx5,cntrdy5)

# Total results
ShiftTot = ShiftLoc1 + ShiftLoc2 + ShiftLoc3 + ShiftLoc4 + ShiftLoc5
medianShiftTot = np.median(ShiftTot)

pyfits.writeto('shift1.fits')
hdus = pyfits.open('shift1.fits')
d=ds9()
d.set_pyfits(hdus)
d.show()

"""Measuring Instrumental Magnitude"""

# -*- coding: utf-8 -*-

"""ASTR 3010 Mini-Project #2,
By Gaith Midani
November 2, 2014"""

# Science target frames:
##D34.fits[986,995][real][object][B_639]:PG1323-086
##D35.fits[986,995][real][object][V_641]:PG1323-086
##D36.fits[986,995][real][object][i_705]:PG1323-086
##D37.fits[986,995][real][object][B_639]:PG1525-071
##D38.fits[986,995][real][object][V_641]:PG1525-071
##D39.fits[986,995][real][object][i_705]:PG1525-071
##D40.fits[986,995][real][object][B_639]:PG1633+099
##D41.fits[986,995][real][object][V_641]:PG1633+099
##D42.fits[986,995][real][object][i_705]:PG1633+099

import pyfits as pf
import matplotlib.pyplot as plt
import numpy as np
import pylab
from cntrd import *
import photutils as ph
import matplotlib.patches as patches

### 1. PREPROCESSING
def fixbads(data,mask,fwhm=2.5,boxsize=11,silent=True,fixneg=True):
  """
  ; Fix bad pixels (from a mask) using Gaussian weighted interpolation.
  ;
  ; Bad pixels are : 
  ;    - pixels with mask value of zero
  ;    - pixels with pixel value of zero or negative
  ; 
  ; Interpolation scheme :
  ;    using a Gaussian kernel window of size '2*boxsize+1' with FWHM
  ;
  ; Badpixel mask :
  ;    0 = good pixels
  ;    1 = bad pixels
  ;
  ; fixneg = False --> Do not attemp to fix negative flux values
  ;
  ; Created by Inseok Song, Sep. 2010
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  """
  from numpy import arange,exp,outer,zeros,real,where,ones
  
  if boxsize > 11.0:
     if not silent: print 'Gaussian kernel size is too large. Reduced to 11.0 pixels!'
     boxsize=11.0

  if fixneg:
     mask[ where(data <= 0) ] = 1
 
  # create a Gaussian kernel of size [-boxsize <--> boxsize] with FWHM
  x=arange(2*boxsize+1)
  y=exp(-0.5*(x-boxsize)**2 / (fwhm/2.354)**2)
  kernel = outer(y,y)

  # Preparing a bigger array with padding on all four sides of the data array.
  #   This step is necessary to correct bad pixels that are too close
  #   to edges. 
  # mask2 is a padded mask array where values are -1 in the padded area.
  BIG = zeros( tuple( real(data.shape) + 2*boxsize ) )
  mask2 = -ones( tuple( real(data.shape) + 2*boxsize ) )

  BIG[boxsize:-boxsize,boxsize:-boxsize] = data.copy()
  mask2[boxsize:-boxsize, boxsize:-boxsize] = mask.copy()
  
  badx, bady = where (mask2 > 0)    # x,y positions of badpixels

  # perform a bad pixel interpolation using the Gaussian weight kernel
  badpos = zip(badx,bady)
  for (ix,iy) in badpos:
      sub_data =  BIG[ix-boxsize:ix+boxsize+1, iy-boxsize:iy+boxsize+1]
      sub_mask = (mask2[ix-boxsize:ix+boxsize+1, iy-boxsize:iy+boxsize+1]).copy()
      bads = where( sub_mask != 0 )
      sub_mask = 1.0 + sub_mask*0
      sub_mask[bads] = 0.0   # so that bad pixels are not used in interpolation
      weight = kernel * sub_mask
      if weight.sum() != 0:
         weight /= weight.sum()
         new_value = (weight*sub_data).sum()
         BIG[ix,iy] = new_value

  return BIG[boxsize:-boxsize,boxsize:-boxsize].copy()

def rdfitsfiles(filelist,silent=False,**kwargs):
   """Read multiple FITS files from a filelist into a 3D array.
   """
   import pyfits,numpy
 
   if not silent: print 'Reading ',filelist[0]
   imdata,header = pyfits.getdata(filelist[0],header=True,**kwargs)
 
   try: nfiles = len(filelist)
   except: return imdata
 
   for file in filelist[1:]:
       if not silent: print 'Reading ',file
       im = pyfits.getdata(file,**kwargs)
       imdata = numpy.dstack((imdata,im))
 
   return imdata

# Bands: [B_639][V_641][i_705]

# Collection of biases define a master bias:

bias_files = ['D46.fits',\
'D47.fits',\
'D48.fits',\
'D49.fits',\
'D50.fits',\
'D51.fits',\
'D52.fits']
biases = rdfitsfiles(bias_files)
master_bias = np.median(biases,axis=2) 


# Collection of flats for each band, to subtract from master bias:

# B-band
Bflat_files = ['D01.fits',\
'D02.fits',\
'D03.fits',\
'D04.fits',\
'D05.fits']
flatsB = rdfitsfiles(Bflat_files) 
Bflat = np.median(flatsB,axis=2) 
Bflat -= master_bias
# V-band
Vflat_files = ['D06.fits',\
'D07.fits',\
'D08.fits',\
'D09.fits',\
'D10.fits']
flatsV = rdfitsfiles(Vflat_files) 
Vflat = np.median(flatsV,axis=2) 
Vflat -= master_bias
# I-band
Iflat_files = ['D11.fits',\
'D12.fits',\
'D13.fits',\
'D14.fits',\
'D15.fits']
flatsI = rdfitsfiles(Iflat_files) 
Iflat = np.median(flatsI,axis=2) 
Iflat -= master_bias

# Creating master flats by division of median:
master_BFlat = Bflat/np.median(Bflat)
master_VFlat = Vflat/np.median(Vflat)
master_IFlat = Iflat/np.median(Iflat)

# Sky files: 
sky_files=['D22.fits',\
'D23.fits',\
'D24.fits',\
'D25.fits',\
'D26.fits',\
'D27.fits']
skys=rdfitsfiles(sky_files)

# Ratioed sky files: last image / first image
ratio=skys[:,:,-1] / skys[:,:,0]
r_avg=np.median(ratio)
r_std=np.std(ratio)
BPM=np.zeros( ratio.shape )
# Bad pixels at +/- 4 standard deviations:
bads=np.where( (ratio < r_avg - 4.0*r_std) | (ratio > r_avg + 4.0*r_std) )
BPM[bads]=1

tgtB = pf.getdata('D43.fits')
tgtV = pf.getdata('D44.fits')
tgtI = pf.getdata('D45.fits')
##PGB = pf.getdata('D34.fits',\
##'D37.fits',\
##'D40.fits')
##PGV = pf.getdata('D35.fits',\
##'D38.fits',\
##'D41.fits')
##PGI = pf.getdata('D36.fits',\
##'D39.fits',\
##'D42.fits')

# Where PGx_A represents PG1323-086 at x band, PGx_B represents PG1525-071
# at x band, and PGx_C represents PG1633+099 at x band:
PGB_A = pf.getdata('D34.fits')
PGB_B = pf.getdata('D37.fits')
PGB_C = pf.getdata('D40.fits')

PGV_A = pf.getdata('D35.fits')
PGV_B = pf.getdata('D38.fits')
PGV_C = pf.getdata('D41.fits')

PGI_A = pf.getdata('D36.fits')
PGI_B = pf.getdata('D39.fits')
PGI_C = pf.getdata('D42.fits')

# Bias subtraction:
tgtB -= master_bias
tgtV -= master_bias
tgtI -= master_bias

PGB_A -= master_bias
PGB_B -= master_bias
PGB_C -= master_bias

PGV_A -= master_bias
PGV_B -= master_bias
PGV_C -= master_bias

PGI_A -= master_bias
PGI_B -= master_bias
PGI_C -= master_bias

# Flat division:
tgtB /= master_BFlat
tgtV /= master_VFlat
tgtI /= master_IFlat

PGB_A /= master_BFlat
PGB_B /= master_BFlat
PGB_C /= master_BFlat

PGV_A /= master_VFlat
PGV_B /= master_VFlat
PGV_C /= master_VFlat

PGI_A /= master_IFlat
PGI_B /= master_IFlat
PGI_C /= master_IFlat

# Bad pixel correction:
final_tgt_B = fixbads(tgtB,BPM)
final_tgt_V = fixbads(tgtV,BPM)
final_tgt_I = fixbads(tgtI,BPM)

final_tgt_PGB_A = fixbads(PGB_A,BPM)
final_tgt_PGB_B = fixbads(PGB_B,BPM)
final_tgt_PGB_C = fixbads(PGB_C,BPM)

final_tgt_PGV_A = fixbads(PGV_A,BPM)
final_tgt_PGV_B = fixbads(PGV_B,BPM)
final_tgt_PGV_C = fixbads(PGV_C,BPM)

final_tgt_PGI_A = fixbads(PGI_A,BPM)
final_tgt_PGI_B = fixbads(PGI_B,BPM)
final_tgt_PGI_C = fixbads(PGI_C,BPM)

## 2. MEASURING INSTRUMENTAL MAGNITUDES, pt. 1

target_B = pf.getdata('D43.fits')
target_V = pf.getdata('D44.fits')
target_i = pf.getdata('D45.fits')
headerB = pf.getheader('D43.fits')
headerV = pf.getheader('D44.fits')
headeri = pf.getheader('D45.fits')

# Zooming in on source of interest
pylab.xlim(580,630)
pylab.ylim(625,665)

# Approximate positions
xapprx, yapprx = (607,645)

# Centroiding
xc,yc = cntrd(target_i,xapprx,yapprx,fwhm=6)
print xc,yc

# Setting aperture radius and sky annuli
ap_rad = 10.0
sky_in = 13.0
sky_out= 20.0

# Defining aperture radius and sky annulus region
aperture = ph.CircularAperture( (xc,yc), r=ap_rad)    
sky_ann  = ph.CircularAnnulus ( (xc,yc), r_in=sky_in, r_out=sky_out)

# Defining sky per pixel
n_pxls_src = np.pi*ap_rad**2
n_pxls_sky = np.pi*(sky_out**2 - sky_in**2)
sky_per_pixel = sky_raw['aperture_sum'] / n_pxls_sky
bg_flx_in_src = sky_per_pixel * n_pxls_src
print sky_per_pixel

# Defining source flux
src_flux = src_raw['aperture_sum'] - bg_flx_in_src
print src_flux

# Exposure time from FITS header
headerB=pf.getheader('D43.fits')
headerV=pf.getheader('D44.fits')
headeri=pf.getheader('D45.fits')
print headerB['EXPTIME']
print headerV['EXPTIME']
print headeri['EXPTIME']

# Measured flux rate of object
src_flux_rate = src_flux / header['EXPTIME']
print 'Source flux rate [cnts/sec] = %s' % src_flux_rate

"""To convert your measured count rates to instrumental magntiudes from the step #2,
you need to use a zero point. Generally, a count rate of 1 cnt/sec is regarded as the base (i.e., "zero") of the measurement.
For this "zero" flux, corresponding magnitude is set to 25 mag.
Then, you need to convert your flux rate by using the relation, m1 - m2 = 2.5*log10(F2/F1)
where m1=25 mag and F1=1 [cnts/sec]."""
# m2 is the instrumental magnitude inside the atmosphere and 2.5*log10(src_flux_rate) is the atmospheric absorption
# ==> m2 = 25 - 2.5*log10(src_flux_rate)

""" This part checks accuracy:"""
##imshow(target_i,interpolation='nearest')
##ylim(615,675)
##xlim(565,655)
##aperture.plot(color='red',lw=1.5,alpha=0.9)   # source aperture in red
##sky_ann.plot(color='green',lw=1.5,alpha=0.9)  # sky annuli in green

##plt.imshow(target_B, interpolation='nearest', origin='lower')
##plt.imshow(target_V, interpolation='nearest', origin='lower')
##plt.imshow(target_i, interpolation='nearest', origin='lower')
##plt.show()


## 3. OBTAINING EXTINCTION AND STANDARD TRANSFORMATION COEFFICIENTS
#Obtain extinction coefficients and standard transformation coefficients for three filters using the measurements from
##2 and catalog magntiudes of each standard stars appeared in StandardFields.pdf.

## Standard transformation magnitudes:

# Color index V for the stars, where A is PG1323-086, B is PG1525-071, C is PG1623+099 and V/BV/UB represent its corresponding bandpass:
AV = 13.481
BV = 15.053
CV = 14.397

# Color index B-V:
ABV = -0.140
BBV = -0.198
CBV = -0.192

# Color index U-B:
AUB = -0.681
BUB = -1.148
CUB = -0.974

# Color index V-R:
AVR = -0.048
BVR = -0.088
CVR = -0.093

# Color index R-I:
ARI = -0.078
BRI = -0.075
CRI = -0.116

# Color index V-I:
AVI = -0.127
BVI = -0.168
CVI = -0.212

# Stand-in for mSTD - m as a function of color index:
funcNumer = mSTD - m

funcAV = (funcNumer, AV) 
funcBV = (funcNumer, BV)
funcCV = (funcNumer, CV)

funcABV = (funcNumer, ABV)
funcBBV = (funcNumer, BBV)
funcCBV = (funcNumer, CBV)

funcAUB = (funcNumer, AUB)
funcBUB = (funcNumer, BUB)
funcCUB = (funcNumer, CUB)

funcAVR = (funcNumer, AVR)
funcBVR = (funcNumer, BVR)
funcCVR = (funcNumer, CVR)

funcARI = (funcNumer, ARI)
funcBRI = (funcNumer, BRI)
funcCRI = (funcNumer, CRI)

funcAVI = (funcNumer, AVI)
funcBVI = (funcNumer, BVI)
funcCVI = (funcNumer, CVI)

# Defining Alpha10 (color coefficient) and Alpha1 (zero-point constant),
# for each star at every color index
Alpha1AV = slope(funcAV)
Alpha10AV = yintercept(funcAV)

Alpha1ABV = slope(funcABV)
Alpha10ABV = yintercept(funcABV)

Alpha1AUB = slope(funcAUB)
Alpha10AUB = yintercept(funcAUB)

Alpha1AVR = slope(funcAVR)
Alpha10AVR = yintercept(funcAVR)

Alpha1ARI = slope(funcARI)
Alpha10ARI = yintercept(funcARI)

Alpha1AVI = slope(funcAVI)
Alpha10AVI = yintercept(funcAVI)

Alpha1BV = slope(funcBV)
Alpha10BV = yintercept(funcBV)

Alpha1BBV = slope(funcBBV)
Alpha10BBV = yintercept(funcBBV)

Alpha1BUB = slope(funcBUB)
Alpha10BUB = yintercept(funcBUB)

Alpha1BVR = slope(funcBVR)
Alpha10BVR = yintercept(funcBVR)

Alpha1BRI = slope(funcBRI)
Alpha10BRI = yintercept(funcBRI)

Alpha1BVI = slope(funcBVI)
Alpha10BVI = yintercept(funcBVI)

Alpha1CV = slope(funcCV)
Alpha10CV = yintercept(funcCV)

Alpha1CBV = slope(funcCBV)
Alpha10CBV = yintercept(funcCBV)

Alpha1CUB = slope(funcCUB)
Alpha10CUB = yintercept(funcCUB)

Alpha1CVR = slope(funcCVR)
Alpha10CVR = yintercept(funcCVR)

Alpha1CRI = slope(funcCRI)
Alpha10CRI = yintercept(funcCRI)

Alpha1CVI = slope(funcCVI)
Alpha10CVI = yintercept(funcCVI)

# Standard transformation coefficients for stars at each color index:
AstdV = 25 + Alpha1AV*(AV) + Alpha10AV
pf.writeto('AstdV.fits')
AstdBV = 25 + Alpha1ABV*(ABV) + Alpha10ABV
pf.writeto('AstdBV.fits')
AstdUB = 25 + Alpha1AUB*(AUB) + Alpha10AUB
pf.writeto('AstdUB.fits')
AstdVR = 25 + Alpha1AVR*(AVR) + Alpha10AVR
pf.writeto('AstdVR.fits')
AstdRI = 25 + Alpha1ARI*(ARI) + Alpha10ARI
pf.writeto('AstdRI.fits')
AstdVI = 25 + Alpha1AVI*(AVI) + Alpha10AVI
pf.writeto('AstdVI.fits')

BstdV = 25 + Alpha1BV*(BV) + Alpha10BV
pf.writeto('BstdV.fits')
BstdBV = 25 + Alpha1BBV*(BBV) + Alpha10BBV
pf.writeto('BstdBV.fits')
BstdUB = 25 + Alpha1BUB*(BUB) + Alpha10BUB
pf.writeto('BstdUB.fits')
BstdVR = 25 + Alpha1BVR*(BVR) + Alpha10BVR
pf.writeto('BstdVR.fits')
BstdRI = 25 + Alpha1BRI*(BRI) + Alpha10BRI
pf.writeto('BstdRI.fits')
BstdVI = 25 + Alpha1BVI*(BVI) + Alpha10BVI
pf.writeto('BstdVI.fits')
 
CstdV = 25 + Alpha1CV*(CV) + Alpha10CV
pf.writeto('CstdV.fits')
CstdBV = 25 + Alpha1CBV*(CBV) + Alpha10CBV
pf.writeto('CstdBV.fits')
CstdUB = 25 + Alpha1CUB*(CUB) + Alpha10CUB
pf.writeto('CstdUB.fits')
CstdVR = 25 + Alpha1CVR*(CVR) + Alpha10CVR
pf.writeto('CstdVR.fits')
CstdRI = 25 + Alpha1CRI*(CRI) + Alpha10CRI
pf.writeto('CstdRI.fits')
CstdVI = 25 + Alpha1CVI*(CVI) + Alpha10CVI
pf.writeto('CstdVI.fits')

""" Extinction coefficients:
 Monochromatic:
 mExt = (m2-25)/X

 Second-order:
 (m2 - 25)/X = k0 + k1*(color index)"""

##filter band 	k0 	k1
##B 	        0.283 	0.049 
##V 	        0.133 	0.032 
##I 	        0.047 	0.000 

# Second-order extinction coefficients, disregarding the constant X:

# Internal term for Star A, at B/V/I bands and different color indices
AmExtBV = 0.283 + 0.049*(AV)
AmExtBBV = 0.283 + 0.049*(ABV)
AmExtBUB = 0.283 + 0.049*(AUB)
AmExtBVR = 0.283 + 0.049*(AVR)
AmExtBRI = 0.283 + 0.049*(ARI)
AmExtBVI = 0.283 + 0.049*(AVI)

AmExtVV = 0.133 + 0.032*(AV)
AmExtVBV = 0.133 + 0.032*(ABV)
AmExtVUB = 0.133 + 0.032*(AUB)
AmExtVVR = 0.133 + 0.032*(AVR)
AmExtVRI = 0.133 + 0.032*(ARI)
AmExtVVI = 0.133 + 0.032*(AVI)

AmExtI = 0.047

# Internal term for Star B, at B/V/I bands and different color indices
BmExtBV = 0.283 + 0.049*(BV)
BmExtBBV = 0.283 + 0.049*(BBV)
BmExtBUB = 0.283 + 0.049*(BUB)
BmExtBVR = 0.283 + 0.049*(BVR)
BmExtBRI = 0.283 + 0.049*(BRI)
BmExtBVI = 0.283 + 0.049*(BVI)

BmExtVV = 0.133 + 0.032*(BV)
BmExtVBV = 0.133 + 0.032*(BBV)
BmExtVUB = 0.133 + 0.032*(BUB)
BmExtVVR = 0.133 + 0.032*(BVR)
BmExtVRI = 0.133 + 0.032*(BRI)
BmExtVVI = 0.133 + 0.032*(BVI)

BmExtI = 0.047

# Internal term for Star C, at B/V/I bands and different color indices
CmExtBV = 0.283 + 0.049*(CV)
CmExtBBV = 0.283 + 0.049*(CBV)
CmExtBUB = 0.283 + 0.049*(CUB)
CmExtBVR = 0.283 + 0.049*(CVR)
CmExtBRI = 0.283 + 0.049*(CRI)
CmExtBVI = 0.283 + 0.049*(CVI)

CmExtVV = 0.133 + 0.032*(CV)
CmExtVBV = 0.133 + 0.032*(CBV)
CmExtVUB = 0.133 + 0.032*(CUB)
CmExtVVR = 0.133 + 0.032*(CVR)
CmExtVRI = 0.133 + 0.032*(CRI)
CmExtVVI = 0.133 + 0.032*(CVI)

CmExtI = 0.047

### 4. OBTAINING INSTRUMENTAL MAGNITUDES
"""Obtain instrumental magnitudes of your science targets (4 stars) from three filter images (four stars are identified in the picture below).
Note that one of four targets was not detected so that you need to obtain 3sigma upper limits.

To convert your measured count rates to instrumental magntiudes from the step #2,
you need to use a zero point. Generally, a count rate of 1 cnt/sec is regarded as the base (i.e., "zero") of the measurement.
For this "zero" flux, corresponding magnitude is set to 25 mag.
Then, you need to convert your flux rate by using the relation, m1 - m2 = 2.5*log10(F2/F1) where m1=25 mag and F1=1 [cnts/sec]."""

# Instrumental magnitude of star A at B/V/I bands for each color index
AmInstBV = AmExtBV - 2.5*log10(src_flux_rate)
pf.writeto('AmInstBV.fits')
AmInstBBV = AmExtBBV - 2.5*log10(src_flux_rate)
pf.writeto('AmInstBBV.fits')
AmInstBUB = AmExtBUB - 2.5*log10(src_flux_rate)
pf.writeto('AmInstBUB.fits')
AmInstBVR = AmExtBVR - 2.5*log10(src_flux_rate)
pf.writeto('AmInstBVR.fits')
AmInstBRI = AmExtBRI - 2.5*log10(src_flux_rate)
pf.writeto('AmInstBVR.fits')
AmInstBVI = AmExtBVI - 2.5*log10(src_flux_rate)
pf.writeto('AmInstBVI.fits')

AmInstVV = AmExtVV - 2.5*log10(src_flux_rate)
pf.writeto('AmInstVV.fits')
AmInstVBV = AmExtVBV - 2.5*log10(src_flux_rate)
pf.writeto('AmInstVBV.fits')
AmInstVUB = AmExtVUB - 2.5*log10(src_flux_rate)
pf.writeto('AmInstVUB.fits')
AmInstVVR = AmExtVVR - 2.5*log10(src_flux_rate)
pf.writeto('AmInstVVR.fits')
AmInstVRI = AmExtVRI - 2.5*log10(src_flux_rate)
pf.writeto('AmInstVRI.fits')
AmInstVVI = AmExtVVI - 2.5*log10(src_flux_rate)
pf.writeto('AmInstVVI.fits')

AmInstI = AmExtI - 2.5*log10(src_flux_rate)
pf.writeto('AmInstI.fits')

# Instrumental magnitude of star B at B/V/I bands for each color index
BmInstBV = BmExtBV - 2.5*log10(src_flux_rate)
pf.writeto('BmInstBV.fits')
BmInstBBV = BmExtBBV - 2.5*log10(src_flux_rate)
pf.writeto('BmInstBBV.fits')
BmInstBUB = BmExtBUB - 2.5*log10(src_flux_rate)
pf.writeto('BmInstBUB.fits')
BmInstBVR = BmExtBVR - 2.5*log10(src_flux_rate)
pf.writeto('BmInstBVR.fits')
BmInstBRI = BmExtBRI - 2.5*log10(src_flux_rate)
pf.writeto('BmInstBRI.fits')
BmInstBVI = BmExtBVI - 2.5*log10(src_flux_rate)
pf.writeto('BmInstBVI.fits')

BmInstVV = BmExtVV - 2.5*log10(src_flux_rate)
pf.writeto('BmInstVV.fits')
BmInstVBV = BmExtVBV - 2.5*log10(src_flux_rate)
pf.writeto('BmInstVBV.fits')
BmInstVUB = BmExtVUB - 2.5*log10(src_flux_rate)
pf.writeto('BmInstVUB.fits')
BmInstBVR = BmExtBVR - 2.5*log10(src_flux_rate)
pf.writeto('BmInstBVR.fits')
BmInstBRI = BmExtBRI - 2.5*log10(src_flux_rate)
pf.writeto('BmInstBRI.fits')
BmInstBVI = BmExtBVI - 2.5*log10(src_flux_rate)
pf.writeto('BmInstBVI.fits')

BmInstI = BmExtI - 2.5*log10(src_flux_rate)
pf.writeto('BmInstI.fits')

# Instrumental magnitude of star C at B/V/I bands for each color index
CmInstBV = CmExtBV - 2.5*log10(src_flux_rate)
pf.writeto('CmInstBV.fits')
CmInstBBV = CmExtBBV - 2.5*log10(src_flux_rate)
pf.writeto('CmInstBBV.fits')
CmInstBUB = CmExtBUB - 2.5*log10(src_flux_rate)
pf.writeto('CmInstBUB.fits')
CmInstBVR = CmExtBVR - 2.5*log10(src_flux_rate)
pf.writeto('CmInstBVR.fits')
CmInstBRI = CmExtBRI - 2.5*log10(src_flux_rate)
pf.writeto('CmInstBRI.fits')
CmInstBVI = CmExtBVI - 2.5*log10(src_flux_rate)
pf.writeto('CmInstBVI.fits')

CmInstVV = CmExtVV - 2.5*log10(src_flux_rate)
pf.writeto('CmInstVV.fits')
CmInstVBV = CmExtVBV - 2.5*log10(src_flux_rate)
pf.writeto('CmInstVBV.fits')
CmInstVUB = CmExtVUB - 2.5*log10(src_flux_rate)
pf.writeto('CmInstVUB.fits')
CmInstVVR = CmExtVVR - 2.5*log10(src_flux_rate)
pf.writeto('CmInstVVR.fits')
CmInstVRI = CmExtVRI - 2.5*log10(src_flux_rate)
pf.writeto('CmInstVRI.fits')
CmInstVVI = CmExtVVI - 2.5*log10(src_flux_rate)
pf.writeto('CmInstVVI.fits')

CmInstI = CmExtI - 2.5*log10(src_flux_rate)
pf.writeto('CmInstI.fits')

# Taking the 3sigma upper limit of an array for star D's instrumental magnitude:
SigmaDmInstB = np.std(final_target_B)
SigmaDmInstV = np.std(final_target_V)
SigmaDmInstI = np.std(final_target_I)

DmInstB = 3*(SigmaDmInstB)
pf.writeto('DmInstB.fits')
DmInstV = 3*(SigmaDmInstV)
pf.writeto('DmInstV.fits')
DmInstI = 3*(SigmaDmInstI)
pf.writeto('DmInstI.fits')

# Value summary: 
print "AmInstBV: ", AmInstBV
print "AmInstBVV: ", AmInstBBV
print "AmInstBUB: ", AmInstBUB
print "AmInstBVR: ", AmInstBVR
print "AmInstBRI: ", AmInstBRI
print "AmInstBVI: ", AmInstBVI

print "AmInstVV: ", AmInstVV
print "AmInstVBV: ", AmInstVBV
print "AmInstVUB: ", AmInstVUB 
print "AmInstVVR: ", AmInstVVR
print "AmInstVRI: ", AmInstVRI
print "AmInstVVI: ", AmInstVVI

print "AmInstI: ", AmInstI

print "BmInstBV: ", BmInstBV 
print "BmInstBBV: ", BmInstBBV
print "BmInstBUB: ", BmInstBUB
print "BmInstBVR: ", BmInstBVR
print "BmInstBRI: ", BmInstBRI
print "BmInstBVI: ", BmInstBVI

print "BmInstVV: ", BmInstVV
print "BmInstVBV: ", BmInstVBV
print "BmInstVVB: ", BmInstVUB
print "BmInstVVR: ", BmInstVVR
print "BmInstVRI: ", BmInstVRI
print "BmInstVVI: ", BmInstVVI

print "BmInstI: ", BmInstI

print "CmInstBV: ", CmInstBV
print "CmInstBBV: ", CmInstBBV
print "CmInstBUB: ", CmInstBUB
print "CmInstBVR: ", CmInstBVR
print "CmInstBRI: ", CmInstBRI
print "CmInstBVI: ", CmInstBVI

print "CmInstVV: ", CmInstVV
print "CmInstVBV: ", CmInstVBV
print "CmInstVUB: ", CmInstVUB
print "CmInstVVR: ", CmInstVVR
print "CmInstVRI: ", CmInstVRI
print "CmInstVVI: ", CmInstVVI

print "CmInstI: ", CmInstI

print "DmInstB: ", DmInstB
print "DmInstV: ", DmInstV
print "DmInstI: ", DmInstI


"""Hnadling FITS data"""

import pyfits 

data = pyfits.getdata('file1.fits')	#Primary Header Data Unit (HDU)
data2 = pyfits.getdata('file1.fits',2)	#2nd extension
data2 = pyfits.getdata('file1.fits',ext=2)
header = pyfits.getdata('file1.fits')
# Get both data and header from a single command
data, hdr = pyfits.getdata('file1.fits', header=True)
# Access a specific keyword
print hdr['ra']
# Add your own keyword
hdr.update('TESTKEY', 'my first keyword')

# Accessing DS9 (go to http://hea-www.harvard.edu/RD/ds9/:

from ds9 import *

#data = pyfits.getdata('file1.fits')
d=ds9()
d.set_np2arr(data) #Will not work. Why?
d.set_n2parr( transpose(data) )
# numpy addresses 2D data with indexing [y,x] instead of [x,y]

# A better method is:
hdus = pyfits.open('file1.fits')
d.set_pyfits(hdus) #This even accesses the WCS (World Coordinate System)

# Pre-processing steps:
# Overscan correction and trimming
# Bias Subtraction
# Dark Subtraction
# Flatfielding

"""Square Gaussian kernel class"""

# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 19:04:34 2014

@author: Gaith
"""

import numpy as np

def makeGaussian(size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)
    
    