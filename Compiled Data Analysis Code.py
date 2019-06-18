# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 17:36:38 2019

@author: gaith
"""

"""COMPILED DATA ANALYSIS CODE"""

"""Curve fitting"""

import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# TASK 1: First, load in your data set
(xfromfile, yfromfile) = sp.loadtxt('group1.dat', unpack=True, usecols=[0,1]) 

### TASK 2: Plot the data
##plt.figure()
##plt.plot(xfromfile,yfromfile)
##plt.suptitle("Raw data set: group1")
###plt.show()

# TASK 3: Come up with a model function called peakfunc
"""We'll need to model the peak as a Gaussian function,
of form a*exp(((-x-b)**2)/(2*c**2)), where
sigma = (height of peak)
8.1 ~ (position of center of the peak, from previous plotting of raw data)
mu = 
... from our data plot, we see that height ~ 17.5, position ~ 8.1, and width ~ 3
so our approximate model function shall be...
"""
def peakfunc(x,mu,sigma):
    return (1/(sigma * sp.sqrt(2 * sp.pi)) * sp.exp( - (x - mu)**2 / (2 * sigma**2)))

# TASK 4: Define a yuncertaintyvector
"""Error of +/-0.7 (V)"""
yuncertainty=.7
plt.errorbar(xfromfile,yfromfile,yerr=yuncertainty,
         fmt='o',ecolor=None,elinewidth=None,
         capsize=3,barsabove=False,lolims=False,
         uplims=False,xlolims=False,xuplims=False)

xsamples = sp.linspace(0,20,100)
yuncertaintyvector=sp.ones(len(xsamples))*yuncertainty

print yuncertaintyvector

plt.plot(xfromfile,yfromfile)
#plt.plot(xfromfile,peakfunc(xfromfile, 0, -10, 3))
#plt.plot(xsamples,peakfunc(xsamples, 0, -10, 3))
#plt.xlim([-.1,20.1])
#plt.ylim([-.1,20.1])
#plt.show()


# TASK 5: Fit your data to your model function.
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

ynoise=noisegenerator(xsamples,'usebias')

yunderlying = peakfunc(xsamples,17.5,1)

yobserved=yunderlying+ynoise

testmeasurements=noisegenerator(sp.zeros(10000),'usebias')

yobservedcorrected=yobserved-sp.mean(testmeasurements)

(popt, punc, rc, d) = curve_fit_custom(peakfunc, xsamples, yobservedcorrected)
abestfitval=popt[0]
abestfitunc=punc[0]
bbestfitval=popt[1]
bbestfitunc=punc[1]
print "abest= %5.3f +/- %5.3f" % (abestfitval,abestfitunc)
print "bbest= %5.3f +/- %5.3f" % (bbestfitval,bbestfitunc)

xdense=sp.linspace(0.,20.,1000.)

ymodelpoints=peakfunc(xsamples,popt[0],popt[1])
ymodelsmooth=peakfunc(xdense,popt[0],popt[1])

residuals=ymodelpoints-yobserved

smoothcurve=peakfunc(xdense,17.5,1)

#plt.plot(xsamples,residuals,'bo')
plt.plot(xsamples, yobservedcorrected, 'go')
plt.plot(xdense,ymodelsmooth,'g-')
plt.plot(xdense,smoothcurve, 'b-')
plt.xlim([-.1,20.1])
plt.ylim([-.1,20.1])
plt.show()

print 'mean=               ', sp.mean(residuals)
print 'standard deviation= ', sp.sqrt(sp.var(residuals))

plt.figure()
(events, edges, patches)=plt.hist(testmeasurements,bins=100,normed=True)
plt.plot(edges, peakfunc(edges,0,1), linewidth=2, color='r')
plt.show()

from scipy.stats import chi2
print chi2.sf(114, 98)

"""Interactive sine series"""

# By Gaith Midani
# Sept. 26, 2014
# PHYS 2001

import numpy as np
import matplotlib.pyplot as plt

y1 = 10*np.pi
y2 = 20*np.pi
y3 = 50*np.pi
n = 40

# Minimum value of n required for >1% accuracy at t = 19pi/2: 40

def defsin1(t, order):
    a = t
    w = a
    for i in range(1, int(y1)):
        a *= -1 * t**2 / ((2 * i) * (2 * i + 1))
        w += a
    return w
    
def defsin2(t, order):
    a = t
    w = a
    for i in range(1, int(y2)):
        a *= -1 * t**2 / ((2 * i) * (2 * i + 1))
        w += a
    return w
    
def defsupersin(t, order):
    a = t
    w = a
    for i in range(1, n):
        a *= -1 * t**2 / ((2 * i) * (2 * i + 1))
        w += a
    return w
    
seriessin1 = np.vectorize(defsin1, excluded=['order'])
seriessin2 = np.vectorize(defsin2, excluded=['order'])
superseriessin = np.vectorize(defsupersin, excluded=['order'])

## First plot parameters
t1 = np.linspace(0, y1, 10000)
yn1 = seriessin1(t1, n)
y1 = np.sin(t1)

## Second plot parameters
t2 = np.linspace(0, y2, 10000)
yn2 = seriessin2(t2, n)
y2 = np.sin(t2)

# Third plot parameters 
t3 = np.linspace(-y3, y3, 10000)
yn3 = superseriessin(t3, n)
y3 = np.sin(t3)

#condExit1 = input('Would you like to continue? Press 1 for (Yes) or press 2 for (No)')
#condExit2 = int(condExit1)

while True:
    query = input('Which plot do you wish to view? Type "1" for seriessin(t,n) vs. sin(t), or type "2" for superseriessin(t,n) vs. sin(t)')
    selection = int(query)

    for query in range(selection):
        if selection == 1:
            query2 = input('Over which range? Type "1" for [0,10pi] or "2" for [0,20pi]')
            selection2 = int(query2)
        elif selection == 2:
            # Third plot
            plt.plot(t3, y3, label='sin(t)')
            plt.plot(t3, yn3, label='superseriessin(t,40)')                
        else:
            print('Please input either 1 or 2')
        break
#        condExit1 = input('Would you like to continue? Press 1 for (Yes) or press 2 for (No)')
#        condExit2 = int(condExit1)
#    if condExit2 == 2:
#        print('End of program')
#        break

    for query2 in range(selection2):
        if selection2 == 1:
            ## First plot
            plt.plot(t1, y1, label='sin(t)')
            plt.plot(t1, yn1, label='seriessin(t,40)')
        elif selection2 == 2:        
            ### Second plot
            plt.plot(t2, y2, label='sin(t)')
            plt.plot(t2, yn2, label='seriessin(t,40)')
        else:
            print('Please input either 1 or 2')
        break
    
#    condExit1 = input('Would you like to continue? Press 1 for (Yes) or press 2 for (No)')
#    condExit2 = int(condExit1)
#    if condExit2 == 2:
#        print('End of program')
#    break5

    plt.ylim([-2.5, 2.5])
    plt.legend()
    plt.show()
    

"""Chaotic system visualization"""
    
# PHYS 4102 Project
# Draft 1
# By Gaith Midani
# May 2, 2016

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Keeping m = 1, we need gamma = F_0/g to yield our desired value... this means
# for pt. 2), F_0 = 1.96
# for pt. 3), F_0 = 10.388
# for pt. 4), F_0 = 10.535 
# for pt. 5), F_0 = 10.5546 (with modification of initial condition, explained below)
b = 1
m = 1
F_0 = 1.96
g = 9.8
L = 0.11 # For part b, simple algebra tells us L (given g = 9.8) needs to be 0.11 in order to equal our desired value of 1.5*omega = 1.5*2pi

# Here, delta is the first-order time derivative of phi
def f(y,t,params):
    phi, delta = y
    twobeta, gamma, omega = params
    derivs = [delta, -delta*Q - ((omega_0)**2)*np.sin(phi) + gamma*((omega_0)**2)*np.cos(omega*t)]
    return derivs

twobeta = b/m
Q = 1/twobeta
gamma = F_0/(m*g)
omega_0 = np.sqrt(g/L) 
omega = 2*np.pi 

# Initial values (same for all parts except part 5, where phi0 = -np.pi/2) 
phi0 = 0.0
delta0 = 0.0

# Grouped parameters
params = [Q, gamma, omega]

# Grouped initial conditions
y0 = [phi0, delta0]

# Time array
tEnd = 50
tStep = 0.05
t = np.arange(0.,tEnd,tStep)

# ODE solution
psoln = odeint(f, y0, t, args = (params,))

fig = plt.figure(1, figsize=(8,8))

ax1 = fig.add_subplot(313)
ax1.plot(t,psoln[:,0])
ax1.set_xlabel('time')
ax1.set_ylabel('phi')
ax1.set_xlim(0.,6)

plt.tight_layout()
plt.show()