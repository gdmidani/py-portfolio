# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 18:26:46 2019

@author: gaith
"""

"""COMPILED EXPERIMENTAL OPTICS CODE"""

"""Tuning up"""

##from numpy import *
##
####import numpy as np
####
##x = array([0.0, 2.0, 4.0, 8.0,])
##y = array([1.1, 1.9, 3.2, 4.0, 5.9])
##yerr = array([0.1, 0.2, 0.1, 0.3, 0.3])
####
####print x
####print y
####print yerr
##
##DataIn = loadtxt('input.dat')
###DataIn1 = column_stack(DataIn)
###print DataIn1
##
##print DataIn

#x, y, yerr = loadtxt('input.dat', unpack=True)

#x, y = loadtxt('input.dat', unpack=True, usecols=[0,1]) 

#print x, y, yerr

"""---"""

##from numpy import *
##
##v = array([0.0, 1.0, 2.0, 3.0, 5.0, 12.0])
##t = array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
##p = 0.15 + v/10.0
##
##savetxt('output.dat', (t,p))
##
##DataOut = column_stack((t,p))
##savetxt('output.dat', DataOut)
##
##savetxt('output.dat', DataOut, fmt=('%3i', '%4.3f'))
##

"""---"""

"""General plotting principles (example: V vs. t)"""

from pylab import *

t = linspace(0.0, 2.0, 100) # makes a list
y = cos(2*pi*t) # makes a list with same length as t

figure() # opens a new figure
plot(t,y,ls='-', color='r', linewidth=1.0)

xlabel('time (s)')
ylabel('voltage (mV)')
title('A Simple Plot')
grid()

scatter(t,y,c='b',s=3,marker='o')
#errorbar(t,y,yerr,terr,marker='s',markersize=2)

semilogx(t,y,ls=':')
semilogy(t,y,marker='.')
loglog(t,y,ls='-')

grid(which='majorminor')
show()

"""Histogram principle"""

from pylab import *

t = array([3,2,2,3,4,5,6,4,3,23,5,54,4,7,90,8,9,9,9,20,4,7,74,37])

figure()
hist(t)
show()

events, edges, patches = hist(t)

"""Na washout, micrometer readings, linear fit for HeNe laser"""

import plotly.plotly as py
py.sign_in(username='limejuice66', api_key='i69hzcv1tr')
import plotly.graph_objs as go

from numpy import arange,array,ones,linspace
from scipy import stats

# For Washout
# Figure out y value of xi
#xi = linspace(0,0.300675,17)
#xi = linspace(0,0.002945,17)
xi = linspace(0,0.0050065,17)
A = array([xi,ones(17)])

# Verbatim, unadjusted
y = [11.87, 10.875, 10.37, 9.36, 8.875, 7.785, 7.405, 6.46, 5.97,
     4.98, 4.085, 3.535, 2.675, 2.23, 1.2, 0.79, 0.28]

slope, intercept, r_value, p_value, std_err = stats.linregress(xi, y)
line = slope*xi+intercept

trace1 = go.Scatter(
                  x=xi, 
                  y=y, 
                  mode='markers',
                  marker=go.Marker(color='rgb(255, 127, 14)'),
                  name='Data'
                  )

trace2 = go.Scatter(
                  x=xi, 
                  y=line, 
                  mode='lines',
                  marker=go.Marker(color='rgb(31, 119, 180)'),
                  name='Fit'
                  )

annotation = go.Annotation(
                  x=3.5,
                  y=23.5,
                  text='$R^2 = 0.9551,\\Y = 0.716X + 19.18$',
                  showarrow=False,
                  font=go.Font(size=16)
                  )
layout = go.Layout(
                title='HeNe: Mirror Displacement Vs. Micrometer Reading',
                plot_bgcolor='rgb(229, 229, 229)',
                  xaxis=go.XAxis(zerolinecolor='rgb(255,255,255)', gridcolor='rgb(255,255,255)'),
                  yaxis=go.YAxis(zerolinecolor='rgb(255,255,255)', gridcolor='rgb(255,255,255)'),
                  annotations=[annotation]
                )

layout = go.Layout(
    title='Na Washout: Mirror Displacement Vs. Micrometer Reading',
    xaxis=dict(
        title='Mirror Displacement [mm]',
        titlefont=dict(
            family='Courier New, monospace',
            size=18,
            color='#7f7f7f'
        )
    ),
    yaxis=dict(
        title='Micrometer Output [mm]',
        titlefont=dict(
            family='Courier New, monospace',
            size=18,
            color='#7f7f7f'
        )
    )
)

data = [trace1, trace2]
fig = go.Figure(data=data, layout=layout)


py.plot(fig, filename='Linear-Fit-IntHeNe')

"""Single slit experiment: intensity visualization"""

import matplotlib.pyplot as plt
from PIL import Image
import numpy

I = plt.imread('W04 S50.tiff')

print I.shape
print I.size

Inew = I.reshape((I.shape[0], -1), order='F')

print Inew.shape
print Inew.size

xlist = numpy.linspace(0,4.65,len(Inew))

InewList = (xlist, Inew)
DataOut = numpy.column_stack(InewList)

print InewList

plt.plot(Inew)
plt.title('4mm-width Single Slit Intensity Profile')
plt.show()



"""Analysis based on Fresnel equations"""

## Fresnel Equations: Plot your data with error bars.

import numpy as np
from operator import truediv
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit

# NOTE: All arrays skip the 90 deg data and goes to 92 degrees instead
xaxis = np.linspace(4,80,20)

"""Vertically polzarized light in mV
It must be established in the report whether this is s- or p- polarization"""
VertParray = np.array([0.12, 0.12, 0.125, 0.12, 0.12, 0.114, 0.106, 0.093,
           0.081, 0.06, 0.042, 0.036, 0.030, 0.028, 0.041, 0.078,
           0.153, 0.350, 0.605, 0.50])

"""Horizontally polarized light in mV"""
HorizParray = np.array([1.30, 1.64, 1.67, 1.68, 1.57, 1.43, 1.18,
                        1.00, 0.76, 0.85, 0.258, 0.242, 0.133, 0.068,
                        0.15, 0.56, 1.46, 3.20, 1.76, 3.60])

print VertParray, HorizParray

"""Evaluate uncertainty via ten measurements at 44 deg"""
VertParray44 = np.array([0.042, 0.042, 0.043, 0.046, 0.047,
                         0.047, 0.046, 0.047, 0.048, 0.046])
HorizParray44 = np.array([0.258, 0.259, 0.259, 0.259, 0.259,
                         0.260, 0.259, 0.258, 0.258, 0.258])

VertPuncertainty = np.std(VertParray44)
HorizPuncertainty = np.std(HorizParray44)

print VertPuncertainty, HorizPuncertainty

"""Vertically polarized light plot"""
plt.errorbar(xaxis,VertParray,yerr=VertPuncertainty,
         fmt='o',ecolor=None,elinewidth=None,
         capsize=3,barsabove=False,lolims=False,
         uplims=False,xlolims=False,xuplims=False)

plt.plot(xaxis,VertParray)
plt.xlabel('Degrees')
plt.ylabel('Power')
plt.title('Input Angle vs. Power (Vertical Polarization)')
plt.show()

"""Horizontally polarized light plot"""
plt.errorbar(xaxis,HorizParray,yerr=HorizPuncertainty,
         fmt='o',ecolor=None,elinewidth=None,
         capsize=3,barsabove=False,lolims=False,
         uplims=False,xlolims=False,xuplims=False)

plt.plot(xaxis,HorizParray)
plt.xlabel('Degrees')
plt.ylabel('Power')
plt.title('Input Angle vs. Power (Horizontal Polarization)')
plt.show()

### Two model functions needed: for s- and p-polarization
"""Model function"""
def Rs(thetaList = [], *args):
    for theta in thetaList:
        n1 = 1   # n value for air
        n2 = 1.5168    # n value for borosilicate glass
        theta_0 = 0
        C = 0
        return ((n1*np.math.cos(theta - (theta_0))-n2*np.math.sqrt(1-(n1/n2)**2*np.math.sin(theta-(theta_0))))/(n1*np.math.cos(theta - (theta_0))+n2*np.math.sqrt(1-(n1/n2)**2*np.math.sin(theta-(theta_0)))))**2

def Rp(thetaList = [], *args):
    for theta in thetaList:
        n1 = 1   # n value for air
        n2 = 1.5168    # n value for borosilicate glass
        theta_0 = 0
        C = 0
        return ((n2*np.math.cos(theta - (theta_0))+ n1*np.math.sqrt(1-(n1/n2)**2*np.math.sin(theta-(theta_0))))/(n2*np.math.cos(theta - (theta_0))+n1*np.math.sqrt(1-(n1/n2)**2*np.math.sin(theta-(theta_0)))))**2


"""Fit data to model function"""
def curve_fit_custom(f, xdata, ydata, p0=None, sigma=None, **kw):
    """
    Pass all arguments to curve_fit, which uses non-linear least squares
    to fit a function, f, to data.  Calculate the uncertaities in the
    fit parameters from the covariance matrix.
    """
    popt, pcov = curve_fit(f, xdata, ydata, p0, sigma, **kw)

    if sigma is None:
        chi2 = sum(((f(xdata,*popt)-ydata))**2)
    else:
        chi2 = sum(((f(xdata,*popt)-ydata)/sigma)**2)
    dof = len(ydata) - len(popt)
    rchi2 = chi2/dof
    print 'results of general_fit:'
    print '   chi squared = ', chi2
    print '   degrees of freedom = ', dof
    print '   reduced chi squared = ', rchi2

    # The uncertainties are the square roots of the diagonal elements
    punc = sp.zeros(len(popt))
    for i in sp.arange(0,len(popt)):
        punc[i] = sp.sqrt(pcov[i,i])
    return popt, punc, rchi2, dof

def noisegenerator(xsamples,control):
    stdev=1 #enter your assinged value for stdev here
    offset=0
    if control is 'usebias':
        usedoffset=offset
    else:#enter your assigned value for the offset here
        usedoffset=0
    if len(xsamples)>1:
        num=len(xsamples)
    else:
        num=1
    noiselist=sp.random.normal(usedoffset,stdev,size=num)
    return noiselist

ynoise=noisegenerator(xaxis,'usebias')

yunderlyingRs = Rs(xaxis)
yunderlyingRp = Rp(xaxis)

yobservedRs=yunderlyingRs+ynoise
yobservedRp=yunderlyingRp+ynoise

testmeasurements=noisegenerator(sp.zeros(10000),'usebias')

yobservedcorrectedRs=yobservedRs-sp.mean(testmeasurements)
yobservedcorrectedRp=yobservedRp-sp.mean(testmeasurements)

pguess = [12,5]
(poptRs, puncRs, rcRs, dRs) = curve_fit_custom(Rs, xaxis, yobservedcorrectedRs, pguess)
abestfitvalRs=poptRs[0]
abestfituncRs=puncRs[0]
bbestfitvalRs=poptRs[1]
bbestfituncRs=puncRs[1]
print "abestRs= %5.3f +/- %5.3f" % (abestfitvalRs,abestfituncRs)
print "bbestRs= %5.3f +/- %5.3f" % (bbestfitvalRs,bbestfituncRs)

(poptRp, puncRp, rcRp, dRp) = curve_fit_custom(Rp, xaxis, yobservedcorrectedRp, pguess)
abestfitvalRp=poptRp[0]
abestfituncRp=puncRp[0]
bbestfitvalRp=poptRp[1]
bbestfituncRp=puncRp[1]
print "abestRp= %5.3f +/- %5.3f" % (abestfitvalRp,abestfituncRp)
print "bbestRp= %5.3f +/- %5.3f" % (bbestfitvalRp,bbestfituncRp)

xdense=sp.linspace(0,92,1000)

ymodelpointsRs=Rs(xaxis,poptRs[0],poptRs[1])
ymodelsmoothRs=Rs(xdense,poptRs[0],poptRs[1])
ymodelpointsRp=Rp(xaxis,poptRp[0],poptRp[1])
ymodelsmoothRp=Rp(xdense,poptRp[0],poptRp[1])


residualsRs=ymodelpointsRs-yobservedRs
residualsRp=ymodelpointsRp-yobservedRp

smoothcurveRs=Rs(xdense)
smoothcurveRp=Rp(xdense)

#plt.plot(xsamples,residuals,'bo')
plt.plot(xaxis, yobservedcorrectedRs, 'go')
plt.plot(xaxis, yobservedcorrectedRp, 'yo')
#plt.plot(xdense,ymodelsmooth,'g-')
#plt.plot(xdense,smoothcurve, 'b-')
#plt.xlim([-.1,20.1])
#plt.ylim([-.1,20.1])
plt.show()

print 'residualsRs=', residualsRs
print 'residualsRp=', residualsRp
print 'meanRs=               ', sp.mean(residualsRs)
print 'meanRp=', sp.mean(residualsRp)
print 'standard deviation of Rs= ', sp.sqrt(sp.var(residualsRs))
print 'standard deviation of Rp= ', sp.sqrt(sp.var(residualsRp))

plt.figure()
(events, edges, patches)=plt.hist(testmeasurements,bins=100,normed=True)
plt.plot(edges, Rs(edges,0,1), linewidth=2, color='r')
plt.plot(edges, Rp(edges,0,1), linewidth=2, color='g')
plt.show()

from scipy.stats import chi2
print chi2.sf(114, 98)

