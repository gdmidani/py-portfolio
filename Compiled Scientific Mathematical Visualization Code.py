# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 18:37:18 2019

@author: gaith
"""

"""Compiled Scientific/Mathematical Visualization Code"""

"""Visualizing sine series"""

import numpy as np
#import matplotlib.pyplot as plt
from matplotlib import *

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

# First plot
plt.plot(t1, y1, label='sin(t)')
plt.plot(t1, yn1, label='seriessin(t, n)')
plt.ylim([-2.5, 2.5])
plt.xlabel('t')
plt.ylabel('f(t)')
plt.legend()
plt.show()

# Second plot
plt.plot(t2, y2, label='sin(t)')
plt.plot(t2, yn2, label='seriessin(t, n)')
plt.ylim([-2.5, 2.5])
plt.xlabel('t')
plt.ylabel('f(t)')
plt.legend()
plt.show()

# Third plot
plt.plot(t3, y3, label='sin(t)')
plt.plot(t3, yn3, label='superseriessin(t, n)') 
plt.ylim([-2.5, 2.5])
plt.xlabel('t')
plt.ylabel('f(t)')
plt.legend()
plt.show()

"""Visualizing slice sums"""

from pylab import *

v=zeros((100,100))
x = 20
y = 20
rsq = 1e-8
rinsq = 0
routsq = 0
bc = (rinsq, routsq)
bcs = ()
vinit = x*x + y*y

def __init__(vinit,x,y):
        vinit.x=x
        vinit.y=y

for i in range(100):
	for j in range(100):
		if ((rsq>=rinsq and rsq<=routsq)):
			v[i,j]=vinit
			bcs.append(i,j)

while True:
	v[1:99,1:99] = (v[1:99,0:98] + v[0:98,1:99] + v[1:99,0:98] + v[0:98,1:99])/4
	for bc in bcs:
		v[bc] = vinit

print i,vnow
imshow(v)
show()

"""Linear algebraic modeling"""

#To base it first on math and physics rather than code, here is another approach attempted.

from pylab import *
import numpy as np
from scipy.linalg import solve

###Variables:
#indicate array of size 100 for the uncharged border length/width and an array of area pi*50^2 for the charged (100V) center, OR
#subtract an array of area pi*50^2 for the charged (100V) center from an array of size 100 for the uncharged border length/width
#In defining an array, variable = [] (Now 'variable' refers to an empty list, or an array)

a = np.zeros((100,100)) #2D array with 100 rows, 100 columns, filled with zeroes
b = np.linspace(0,1,50) #Array with 50 points (or millimeters) in between 0 and 1 (arbitrary selection) 


x = 20
y = 20
h = 1 
i = 1
j = 1 
#rho = 1 #? 
#array name = array execution:
#Attempt 1, meant to represent the total structure:
print (a+b)

#Big question: how to relate the returned (a+b) to the V(i,j) function -- which in turn is prospectively to be
#related to the Jacobi algorithm?

#This is the key to measuring the value of a center as directed by its four adjacent components. To bracket it also implicates an array -- one which must be considered now in the variable.
#V(i,j) = (1/4)*(f(x+h,y)+f(x-h,y)+f(x,y+h)+f(x,y-h))

##What is 'f'? You may have to break down the relaxation method into components, define, then execute, such as:

#f(x,y) = (1/4)*(f1+f2+f3+f4)

#...This still doesn't answer what f is. Perhaps it is a function pertaining to potential.

##If the computer complains, then try:
#def f(x,y,h,n):
#	for i in range(n):
#		i += 1
#	return f(20,20,1)

#Attempt 2:
        
f1 = f(x+h,y)
f2 = f(x-h,y)
f3 = f(x,y+h)
f4 = f(x,y-h)

def f():
	if i <= n:
		return (1/4)*(f1 + f2 + f3 + f4)



#Another way to articulate this would be through the Jacobi method:
#f(i,j) = (1/4)*(f([i+1,j])+f([i-1,j])+f([i,j+1])+f([i,j-1])+4*pi*h**2*str(rho)(i,j))

###Jacobi practice:
#A = represents the matrix
#b = represents the solution to Ax=b
#x = what we are attempting to solve for (we first make an initial guess)
#n = number of iterations

#This may need to be applied in a loop for context
#'x' renamed to 'k'
def jacobi(A, b, k, n):
	D=np.diag(A)
	R = A - np.diagflat(D)

        for i in range(n):
                k = (b - np.dot(R,k))/ D
        return k

#A = np.array([[4.0, -2.0, 1.0], [1.0, -3.0, 2.0], [-1.0, 2.0, 6.0]]) #use a defined array instead
A = a + b
b = [1.0, 2.0, 3.0] #b will need to represent the solution to V[i,j]
k = [1.0, 1.0, 1.0] 
n = 25

print solve(A, b)
k = jacobi(A, b, k, n)

###Visual modeling:
#print v
print f()
imshow(V)
show()

"""Quantum mechanical solution with Hermite polynomial"""
#QM 10.1

import numpy as np
from matplotlib.pyplot import plt
import scipy as sp

#a

m = 1
w = 1
h = 1
thing = np.sqrt(m*w/h)

Pa = ((1/(16*np.sqrt(2)))*np.sqrt(thing/np.pi)*sp.special.hermite(2)**2

"""Quantum mechanical wavefunction probability solution"""

import numpy as np
import matplotlib.pyplot as plt

##a = 5.29e-11
##z = np.arange(0,100)
##fu = [z**2 for z in range(0,100)]
##P = np.array((1/(2*(a**2)))*[(1/2)*(np.exp(-2*z/a))*(a+2*z)
##                    + (2*z/(9*np.sqrt(2)*a))*(np.cos(0))*(np.exp(-3*z/a))*(a+3*z)
##                    + ((np.power(z,2))/(8*(a**2)))*(np.exp(-z/a))*(a+z)])

def P(z):
    fu = [z**2 for z in range(0,100)]
    z += 0
    P += (1/(2*(a**2)))*[(1/2)*(np.exp(-2*z/a))*(a+2*z)
                    + (2*z/(9*np.sqrt(2)*a))*(np.cos(0))*(np.exp(-3*z/a))*(a+3*z)
                    + ((z**2)/(8*(a**2)))*(np.exp(-z/a))*(a+z)]
    P.append()
    return P

plt.plot(P(100))
plt.show()

"""Classical system with circular motion"""

# P 6.18 Draft 4 rewrite v2

import numpy as np
import matplotlib.pyplot as plt

g = 9.8
# Alpha
a = np.radians(20)
# Initial position r_0
r = 0.1
subr = 0.999*r
supr = 1.001*r
# Initial velocity
r1 = 0
# Initial angle theta
th = 0
# Time bounds
t1 = 0 
t2 = 5
# Small interval epsilon
e = 0.001
# Omega for which motion is circular
w = np.sqrt(g/(r*(np.tan(a))))

T1 = np.arange(t1,t2,e)

rlist = []
for i in T1:
    r2 = r*(w**2)*(np.sin(a)**2)-g*(np.cos(a))*(np.sin(a))
    r = r + e*r1
    r1 = r1 + e*r2
    rlist.append(r)

T1 = np.arange(t1,t2,e)

plt.plot(T1,rlist,label='r(t) given r(0)=r_0')
plt.xlabel('t')
plt.ylabel('r(t)')
plt.legend()
plt.show()

subrlist = []
for i in T1:
    r2 = subr*(w**2)*(np.sin(a)**2)-g*(np.cos(a))*(np.sin(a))
    r = subr + e*r1
    r1 = r1 + e*r2
    subrlist.append(r)

print subrlist

plt.plot(T1,subrlist,label='r(t) given r(0)=0.999r_0')
plt.xlabel('t')
plt.ylabel('r(t)')
plt.legend()
plt.show()

suprlist = []
for i in T1:
    r2 = supr*(w**2)*(np.sin(a)**2)-g*(np.cos(a))*(np.sin(a))
    r = supr + e*r1
    r1 = r1 + e*r2
    suprlist.append(r)

plt.plot(T1,suprlist,label='r(t) given r(0)=1.001r_0')
plt.xlabel('t')
plt.ylabel('r(t)')
plt.legend()
plt.show()

print rlist, subrlist, suprlist

"""Classical mechanical system with circular motion, alt. version"""

# P 6.18 Draft 4 rewrite v2

import numpy as np
import matplotlib.pyplot as plt

g = 9.8
# Alpha
a = np.radians(20)
# Initial position
#r = 0.1
#subr = 0.999*r
#supr = 1.001*r
# Initial velocity
r1 = 0
# Initial angle theta
th = 0
# Time bounds
t1 = 0 
t2 = 3
# Small interval epsilon
e = 0.1

# Omega for which motion is circular
w = np.sqrt(g/(r*(np.tan(a))))


rlist = []
for i in np.arange(t1,t2,e):
        r2 = 0.1*(w**2)*(np.sin(a)**2)-g*(np.cos(a))*(np.sin(a))
        r = 0.1 + e*r1
        r1 = r1 + e*r2
        rlist.append(r)

T1 = np.arange(t1,t2,e)

plt.plot(T1,rlist)
plt.show()

subrlist = []
for i in np.arange(t1,t2,e):
    r2 = 0.999*0.1*(w**2)*(np.sin(a)**2)-g*(np.cos(a))*(np.sin(a))
    r = 0.999*0.1 + e*r1
    r1 = r1 + e*r2
    subrlist.append(r)

print subrlist

plt.plot(T1,subrlist)
plt.show()

suprlist = []
for i in np.arange(t1,t2,e):
        r2 = 1.001*0.1*(w**2)*(np.sin(a)**2)-g*(np.cos(a))*(np.sin(a))
        r = 1.001*0.1 + e*r1
        r1 = r1 + e*r2
        suprlist.append(r)

plt.plot(T1,suprlist)
plt.show()

"""Classical mechanical system with circular motion, yet another version"""

# rough draft 2 for phys 4101, problem 6.18

import math
import numpy as np 
import matplotlib.pyplot as plt

g = 9.8
alpha = math.radians(20)
r_0 = 0.1
t = np.arange(0,20,1)
tstep = 20

def omega(t):
        w = math.sqrt(g/(r_0*math.tan(alpha)))
        if t.any() == 0:
                return w

def r_first(t):
	r1 = 0.999*r_0
	R1 = r1*(omega(t)**2)*(math.sin(alpha))**2 - g*math.cos(alpha)*math.sin(alpha)
	return R1

#Strange result is a consequence of the perfect conditions for circular motion being met
def r_second(t):
	R2 = r_0*(omega(t)**2)*(math.sin(alpha))**2 - g*math.cos(alpha)*math.sin(alpha)
	return R2

def r_third(t):
	r3 = 1.001*r_0
	R3 = r3*(omega(t)**2)*(math.sin(alpha))**2 - g*math.cos(alpha)*math.sin(alpha)
	return R3

def theta(t):
	th = 0 
	if t.any() == 0:
		return th

def r_dot(t):
	rd = 0
	if t.any() == 0:
		return rd

## Verified as working:
#print omega(t)
#print r_first(t), r_second(t), r_third(t)
#print theta(t), r_dot(t)

## From here on out, plots don't show any function.
## Probably need to generalize defintion of r's.

# First plot
plt.plot(t, r_first(t), label='r(t) vs. t for r(0) = 0.999r_0')
plt.ylim([0, 0.005])
plt.xlabel('t')
plt.ylabel('r(t)')
plt.legend()
plt.show()

# Second plot
plt.plot(t, r_second(t), label='r(t) vs. t for r(0) = r_0')
plt.ylim([0, 0.005])
plt.xlabel('t')
plt.ylabel('r(t)')
plt.legend()
plt.show()

# Third plot
plt.plot(t, r_first(t), label='r(t) vs. t for r(0) = 1.001r_0')
plt.ylim([0, 0.005])
plt.xlabel('t')
plt.ylabel('r(t)')
plt.legend()
plt.show()

"""Relative masses from equation of motion"""

# Pr6.27 Draft 3, Full rewrite
# No bugs, but result isn't completely satisfactory
# (Should have a ratio of .208 for m/M = 1/10)

from itertools import *
import numpy as np

m = 1
M = 10
g = 9.8
# Initial value of r
r = 0.5
# Initial velocity 
r1 = 0
# Initial value of theta 
th = 0
# Small interval epsilon
e = 0.01
# Initial and final times
t1 = 0
t2 = 3
#L = 1
L = np.sqrt(M*m*g*(r**2))

# r2 below indicates acceleration
# th1 indicates t-derivative of theta

# Equation of motion for r
rlist = []
for i in np.arange(t1,t2,e):
        r2 = (L/(m*(r**3)) - M*g) / (M + m)
        r = r + e*r1
        r1 = r1 + e*r2
        rlist.append(r)

#Rlist = np.vectorize(R, excluded=['order'])

# Equation of motion for theta (th)
thlist = []
for i in np.arange(t1,t2,e):
    th1 = L/(m*(r**2))
    th = th + e*th1
    thlist.append(th)

rmin = min(rlist)

ratio = rmin/r

print ratio

"""Expanding unknowns"""

# practice example

import numpy as np

x = 2
x1 = 0
e = 0.1
t1 = 0
t2 = 3

#t_step = e

mylist = []
for i in np.arange(t1,t2,e):
        x2 = -5*x
        x = x+e*x1
        x1 = x1+e*x2
        mylist.append(x)
#        nums = list(x)
#        print x

print mylist

"""Visualizing electromagnetic system"""

## EM1 PS4 Problem 3, DRAFT 2
## By Gaith Midani

import numpy as np
import matplotlib.pyplot as plt

##Here is the field point you should use for problem 1 (the computer
##problem) on Problem Set #5.  Each component is in meters.
##
##x               y               z
##
##-1.250  -0.860  -0.680

Q = 0.01
x = -1.250
X = 0 
y = -0.860
Y = 0
z = -0.680
Z = 0
twoL = 0.1
k = 8.988e9 
N = 35

##Elist = []
##for i in np.arange(0,N,1):
##    charge = Q/(N**2)
##    Elist.append(charge)
##
##E_x = x*Elist
##E_y = y*Elist
##E_z = z*Elist
##
##print E_x, E_y, E_z


Chargelist = []
for i in np.arange(0,N,1):
    charge = Q/(N**2)
    Chargelist.append(charge)

Xlist = []
for i in np.arange(0,N,1):
    X += x-N
    Xlist.append(X)

X = np.ndarray(N)
X.fill(x)
Xtot = X - Xlist

Ylist = []
for i in np.arange(0,N,1):
    Y += y-N
    Ylist.append(Y)

Y = np.ndarray(N)
Y.fill(y)
Ytot = Y - Ylist

##Zlist = []
##for i in np.arange(0,N,1):
##    Z += z 
##    Zlist.append(Z)

Ztot = np.ndarray(N)
Ztot.fill(z)

print Chargelist, Xtot, Ytot, Ztot

Eijx = k*(Xtot)/((Xtot**2 + Ytot**2 + Ztot**2)**(3/2))
Eijy = k*(Ytot)/((Xtot**2 + Ytot**2 + Ztot**2)**(3/2))
Eijz = k*(Ztot)/((Xtot**2 + Ytot**2 + Ztot**2)**(3/2))

Enum = (Eijx, Eijy, Eijz)
print Enum

## Calculate on paper!
#Eexact = something

absEdiff = abs(Enum-Eexact)
Err = absEdiff/Eexact

plt.plot(Ztot,Err,label='abs(Eapprox-Eexact)/Eexact vs. z')
plt.xlabel('z')
plt.ylabel('abs(Eapprox-Eexact)/Eexact')
plt.legend()
plt.show()

"""Another electromagnetic system"""

## EM1 PS4 Problem 3, DRAFT 1

import numpy as np

Q = 0.01
x = 1
y = 2
z = 3
twoL = 0.1
N = 35

Elist = []
for i in np.arange(0,N,1):
    charge = Q/(N**2)
    Elist.append(charge)

E_x = x*Elist
E_y = y*Elist
E_z = z*Elist

print E_x, E_y, E_z

"""Momentum probability density v. momentum"""

import numpy as np
import matplotlib.pyplot as plt

L = 1
p = np.linspace(-100, 100)
n = 8
h = 6.582*10**(-16)

fxncoeff = (1/((np.pi*L)*(((n**2)*(np.pi**2)*(h**2) - (L**2)*(p**2))**2)))
fxn = (((n*np.pi*L*(h**2))**2) - 2*((n*np.pi*L*(h**2))**2)*np.cos(n*np.pi)
       *np.cos(L*p/h) + ((n*np.pi*L*(h**2))**2)*((np.cos(n*np.pi))**2) + (h**2)*(L**4)*(p**2)*((np.sin(n*np.pi))**2))

result = fxncoeff * fxn

plt.suptitle('Momentum probability density vs momentum, with L fixed at 1; n = 8')
plt.plot(result)
plt.show()
print result

"""Specific figure"""

import matplotlib.pyplot as plt
import numpy as np

x = np.arange(0,100,100)
L = 7

##fig = (1/(5*L))*((8*((np.sin((np.pi*x)/L))**2)) + 8*(np.sin((np.pi*x)/L))*
##                 (np.sin((2*np.pi*x)/L)) + 2*((np.sin((2*np.pi*x)/L))**2))

fig = (1/5)*(8*((np.sin**2(np.pi*x))) + (8*(np.sin(np.pi*x))*(np.sin(2*np.pi*x))) + (2*((np.sin**2(2*np.pi*x)))))

plt.plot(fig)
plt.show()

"""Another specific figure"""

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0,50,50)

fig = (10**(-(1/2)))*((np.cosh(x/5))**(-1))

plt.plot(fig)
plt.xlabel('x')
plt.ylabel('$\psi_2$(x)')
plt.show()

"""Yet another specific figure"""

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0,100,100)

fig = (4/(np.pi))*(1/(x**2 + 4))

plt.plot(fig)
plt.xlabel('x')
plt.ylabel('$\psi_1$(x)')
plt.show()

"""Incomplete project on potential"""

#Project 3 draft
import numpy as np
from matplotlib import imshow
import math

#a = external length in mm
a = 100
#b = external height in mm
b = 100
#c = radius in mm of internal circle
r = 25
#x, y being locations:
x = 20
y = 20

square = a*b
circle = 3.14*r**2

#Create 2D array of booleans where potential is to be held fixed

#myArray=np.ndarray((100,100))

class Potential:
	def __init__(self, circle):
		self.circle = array(pi*r**2)

class NoPotential:
	def __init__(self, square):
		self.square = array([a,b])

matplotlib.imshow(Potential)
matplotlib.imshow(NoPotential)

#Decide on what range of potential values to work with

#Graphical representation of potential values: create 1D array of 100+ colors, i.e.
#vary from black for lowest potential to white for highest. Use "new" once to create array,
#then write a loop to create each separate Color orbject (either call Color a constructor
#method with "new" or call a static method such as getHSBcolor

#Tack on temporary data into the array of potentials, then write an appropriate
#"paint" method to color each square according to its potential

#Once graphics are working, write a method executing one iteration of the relaxation
#algorithm

#Create a button to call your relaxation method just once


#Loop 1: walking through an array line
#Loop 2: repeat loop 1 throughout array

"""Project on classical mechanical system"""

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

"""Riemann sum plot"""

import numpy as np

def R(t):
    return t**7-2*t**2+0.25*t+0.5

def riemannplot(fct,a,b,n):
    if n<=0:
        print "Error: n needs to be positive"
        return False
    smoothh = (b-a)/100.0
    x = arange(a,b+smoothh,smoothh)
    plot(x,fct(x))
    h = (double(b)-double(a))/double(n)
    riemannx=arange(a,b,h)
    riemanny=fct(riemannx)
    bar(riemannx,riemanny,width=h,alpha=0.5,facecolor='orange')
    xlabel('x')
    ylabel('y')
    title('Riemann Left Sum for R(t)')
    xlim(-1,2)
    show()


np.riemannplot(f,-1,2,8)

"""Malus' law"""

##Make an array which is the ratio of the two measurements,
##at each theta. Calculate your uncertainty in the
##ratio by finding the standard deviation of the 10 ratios
##2
##you took at a fixed angle in the last step of the procedure.
##Plot this ratio versus angle, with error bars. 

import numpy as np
from operator import truediv
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit
from scipy.stats import cosine
import math

# NOTE: All arrays skip the 90 deg data and goes to 92 degrees instead
xaxis = np.linspace(0,92,24)

"""Power through the polarizer: first element
in the pairs that occur in the array"""
Parray = sp.array([0.36, 0.45, 0.45, 0.53, 0.55, 0.56, 0.58, 0.58,
           0.59, 0.59, 0.58, 0.57, 0.56, 0.51, 0.47, 0.47,
           0.44, 0.41, 0.34, 0.30, 0.28, 0.21, 0.19, 0.12])

"""Power without the polarizer; second element
in the pairs that occur in the array"""
NoParray = sp.array([1.69, 1.69, 1.69, 1.87, 1.86, 1.86, 1.87, 1.87,
                     1.86, 1.86, 1.86, 1.87, 1.86, 1.86, 1.86, 1.86,
                     1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 1.86])

RatioArray = map(truediv, Parray, NoParray)
print RatioArray

##Puncertainty = np.std(Parray)
##NoPuncertainty = np.std(NoParray)
##
##print Puncertainty, NoPuncertainty

Parray40 = sp.array([0.58, 0.60, 0.60, 0.60, 0.59])
NoParray40 = sp.array([1.86, 1.86, 1.86, 1.87, 1.86])
RatioArray40 = map(truediv, Parray40, NoParray40)

Puncertainty = sp.std(Parray40)     # Uncertainty of 0.008
NoPuncertainty = sp.std(NoParray40) # Uncertainty of 0.004
Ratiouncertainty = sp.std(RatioArray40)

#print Puncertainty, NoPuncertainty, Ratiouncertainty
print 'ratio uncertainty=', Ratiouncertainty

"""To plot this ratio vs angle with error bars, we first need to
define the error bars"""

# with polarizer
plt.errorbar(xaxis,Parray,yerr=Puncertainty,
         fmt='o',ecolor=None,elinewidth=None,
         capsize=3,barsabove=False,lolims=False,
         uplims=False,xlolims=False,xuplims=False)

plt.plot(xaxis,Parray)
plt.xlabel('Degrees')
plt.ylabel('Power')
plt.title('Input Angle vs. Power (With Polarizer)')
plt.show()

# without polarizer
plt.errorbar(xaxis,NoParray,yerr=NoPuncertainty,
         fmt='o',ecolor=None,elinewidth=None,
         capsize=3,barsabove=False,lolims=False,
         uplims=False,xlolims=False,xuplims=False)

plt.plot(xaxis,NoParray)
plt.xlabel('Degrees')
plt.ylabel('Power')
plt.title('Input Angle vs. Power (Without Polarizer)')
plt.show()

# Ratio
plt.errorbar(xaxis,RatioArray,yerr=Ratiouncertainty,
         fmt='o',ecolor=None,elinewidth=None,
         capsize=3,barsabove=False,lolims=False,
         uplims=False,xlolims=False,xuplims=False)

plt.plot(xaxis,RatioArray)
plt.xlabel('Degrees')
plt.ylabel('Power')
plt.title('Input Angle vs. Power (Ratio of [With Polarizer]/[Without Polarizer])')
plt.show()

"""Model function"""
def Pvec(thetalist = []):
    P_0 = np.math.radians(0.21301775147928995)
#    P_0 = np.math.radians(xaxis)
    theta_0 = 0
    C = 0
    #thetaList = []
    thetalist.append(0)
    thetalist.append(np.math.radians(4))
    thetalist.append(np.math.radians(8))
    thetalist.append(np.math.radians(12))
    thetalist.append(np.math.radians(16))
    thetalist.append(np.math.radians(20))
    thetalist.append(np.math.radians(24))
    thetalist.append(np.math.radians(28))
    thetalist.append(np.math.radians(32))
    thetalist.append(np.math.radians(36))
    thetalist.append(np.math.radians(40))
    thetalist.append(np.math.radians(44))
    thetalist.append(np.math.radians(48))
    thetalist.append(np.math.radians(52))
    thetalist.append(np.math.radians(56))
    thetalist.append(np.math.radians(60))
    thetalist.append(np.math.radians(64))
    thetalist.append(np.math.radians(68))
    thetalist.append(np.math.radians(72))
    thetalist.append(np.math.radians(76))
    thetalist.append(np.math.radians(80))
    thetalist.append(np.math.radians(84))
    thetalist.append(np.math.radians(86))
    thetalist.append(np.math.radians(90))
    return (P_0*((np.math.cos(theta - (theta_0)))*(np.math.cos(theta - (theta_0))))+C)
        #PvecList = list(Pvec(theta))
        #return (PvecList)
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

yunderlying = Pvec(xaxis)

yobserved=yunderlying+ynoise

testmeasurements=np.array(noisegenerator(sp.zeros(10000),'usebias'))

yobservedcorrected=yobserved-sp.mean(testmeasurements)

pguess = [12,5]
(popt, punc, rc, d) = curve_fit_custom(Pvec, xaxis, yobservedcorrected, pguess)
abestfitval=popt[0]
abestfitunc=punc[0]
bbestfitval=popt[1]
bbestfitunc=punc[1]
print "abest= %5.3f +/- %5.3f" % (abestfitval,abestfitunc)
print "bbest= %5.3f +/- %5.3f" % (bbestfitval,bbestfitunc)

xdense=sp.linspace(0,92,1000)

ymodelpoints=Pvec(xaxis,popt[0],popt[1])
ymodelsmooth=Pvec(xdense,popt[0],popt[1])

residuals=ymodelpoints-yobserved

smoothcurve=Pvec(xdense)

#plt.plot(xsamples,residuals,'bo')
plt.plot(xaxis, yobservedcorrected, 'go')
#plt.plot(xdense,ymodelsmooth,'g-')
#plt.plot(xdense,smoothcurve, 'b-')
#plt.xlim([-.1,20.1])
#plt.ylim([-.1,20.1])
plt.show()

print 'residuals=',residuals
print 'mean=               ', sp.mean(residuals)
print 'standard deviation= ', sp.sqrt(sp.var(residuals))

"""plt.autoscale(enable-True, axis=u 'both', tight=False)"""

plt.figure()
#(events, edges, patches)=plt.hist(testmeasurements,bins=np.array([0,20,40,60,80,100]),normed=True)
(events, edges, patches)=plt.hist(testmeasurements,bins=23,normed=True)
#plt.plot(edges, Pvec(edges,0,1), linewidth=2, color='r')
print 'edges = ', edges
print 'Pvec(edges) = ', Pvec(edges)
#plt.plot(edgess,Pvec(edges,0,1),linewidth=2,color='r')
plt.show()

from scipy.stats import chi2
print chi2.sf(114, 98)

"""Malus' law, alt."""

##Make an array which is the ratio of the two measurements,
##at each theta. Calculate your uncertainty in the
##ratio by finding the standard deviation of the 10 ratios
##2
##you took at a fixed angle in the last step of the procedure.
##Plot this ratio versus angle, with error bars. 

import numpy as np
from operator import truediv
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit
from scipy.stats import cosine
import math

# NOTE: All arrays skip the 90 deg data and goes to 92 degrees instead
xaxis = np.linspace(0,92,24)

"""Power through the polarizer: first element
in the pairs that occur in the array"""
Parray = sp.array([0.36, 0.45, 0.45, 0.53, 0.55, 0.56, 0.58, 0.58,
           0.59, 0.59, 0.58, 0.57, 0.56, 0.51, 0.47, 0.47,
           0.44, 0.41, 0.34, 0.30, 0.28, 0.21, 0.19, 0.12])

"""Power without the polarizer; second element
in the pairs that occur in the array"""
NoParray = sp.array([1.69, 1.69, 1.69, 1.87, 1.86, 1.86, 1.87, 1.87,
                     1.86, 1.86, 1.86, 1.87, 1.86, 1.86, 1.86, 1.86,
                     1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 1.86])

RatioArray = map(truediv, Parray, NoParray)
print RatioArray

##Puncertainty = np.std(Parray)
##NoPuncertainty = np.std(NoParray)
##
##print Puncertainty, NoPuncertainty

Parray40 = sp.array([0.58, 0.60, 0.60, 0.60, 0.59])
NoParray40 = sp.array([1.86, 1.86, 1.86, 1.87, 1.86])
RatioArray40 = map(truediv, Parray40, NoParray40)

Puncertainty = sp.std(Parray40)     # Uncertainty of 0.008
NoPuncertainty = sp.std(NoParray40) # Uncertainty of 0.004
Ratiouncertainty = sp.std(RatioArray40)

#print Puncertainty, NoPuncertainty, Ratiouncertainty
print 'ratio uncertainty=', Ratiouncertainty

"""To plot this ratio vs angle with error bars, we first need to
define the error bars"""

# with polarizer
plt.errorbar(xaxis,Parray,yerr=Puncertainty,
         fmt='o',ecolor=None,elinewidth=None,
         capsize=3,barsabove=False,lolims=False,
         uplims=False,xlolims=False,xuplims=False)

plt.plot(xaxis,Parray)
plt.xlabel('Degrees')
plt.ylabel('Power')
plt.title('Input Angle vs. Power (With Polarizer)')
plt.show()

# without polarizer
plt.errorbar(xaxis,NoParray,yerr=NoPuncertainty,
         fmt='o',ecolor=None,elinewidth=None,
         capsize=3,barsabove=False,lolims=False,
         uplims=False,xlolims=False,xuplims=False)

plt.plot(xaxis,NoParray)
plt.xlabel('Degrees')
plt.ylabel('Power')
plt.title('Input Angle vs. Power (Without Polarizer)')
plt.show()

# Ratio
plt.errorbar(xaxis,RatioArray,yerr=Ratiouncertainty,
         fmt='o',ecolor=None,elinewidth=None,
         capsize=3,barsabove=False,lolims=False,
         uplims=False,xlolims=False,xuplims=False)

plt.plot(xaxis,RatioArray)
plt.xlabel('Degrees')
plt.ylabel('Power')
plt.title('Input Angle vs. Power (Ratio of [With Polarizer]/[Without Polarizer])')
plt.show()

"""Model function"""
def Pvec(thetalist = []):
    P_0 = np.math.radians(0.21301775147928995)
#    P_0 = np.math.radians(xaxis)
    theta_0 = 0
    C = 0
    #thetaList = []
    thetalist.append(0)
    thetalist.append(np.math.radians(4))
    thetalist.append(np.math.radians(8))
    thetalist.append(np.math.radians(12))
    thetalist.append(np.math.radians(16))
    thetalist.append(np.math.radians(20))
    thetalist.append(np.math.radians(24))
    thetalist.append(np.math.radians(28))
    thetalist.append(np.math.radians(32))
    thetalist.append(np.math.radians(36))
    thetalist.append(np.math.radians(40))
    thetalist.append(np.math.radians(44))
    thetalist.append(np.math.radians(48))
    thetalist.append(np.math.radians(52))
    thetalist.append(np.math.radians(56))
    thetalist.append(np.math.radians(60))
    thetalist.append(np.math.radians(64))
    thetalist.append(np.math.radians(68))
    thetalist.append(np.math.radians(72))
    thetalist.append(np.math.radians(76))
    thetalist.append(np.math.radians(80))
    thetalist.append(np.math.radians(84))
    thetalist.append(np.math.radians(86))
    thetalist.append(np.math.radians(90))
    return (P_0*((np.math.cos(theta - (theta_0)))*(np.math.cos(theta - (theta_0))))+C)
        #PvecList = list(Pvec(theta))
        #return (PvecList)
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

yunderlying = Pvec(xaxis)

yobserved=yunderlying+ynoise

testmeasurements=np.array(noisegenerator(sp.zeros(10000),'usebias'))

yobservedcorrected=yobserved-sp.mean(testmeasurements)

pguess = [12,5]
(popt, punc, rc, d) = curve_fit_custom(Pvec, xaxis, yobservedcorrected, pguess)
abestfitval=popt[0]
abestfitunc=punc[0]
bbestfitval=popt[1]
bbestfitunc=punc[1]
print "abest= %5.3f +/- %5.3f" % (abestfitval,abestfitunc)
print "bbest= %5.3f +/- %5.3f" % (bbestfitval,bbestfitunc)

xdense=sp.linspace(0,92,1000)

ymodelpoints=Pvec(xaxis,popt[0],popt[1])
ymodelsmooth=Pvec(xdense,popt[0],popt[1])

residuals=ymodelpoints-yobserved

smoothcurve=Pvec(xdense)

#plt.plot(xsamples,residuals,'bo')
plt.plot(xaxis, yobservedcorrected, 'go')
#plt.plot(xdense,ymodelsmooth,'g-')
#plt.plot(xdense,smoothcurve, 'b-')
#plt.xlim([-.1,20.1])
#plt.ylim([-.1,20.1])
plt.show()

print 'residuals=',residuals
print 'mean=               ', sp.mean(residuals)
print 'standard deviation= ', sp.sqrt(sp.var(residuals))

"""plt.autoscale(enable-True, axis=u 'both', tight=False)"""

plt.figure()
#(events, edges, patches)=plt.hist(testmeasurements,bins=np.array([0,20,40,60,80,100]),normed=True)
(events, edges, patches)=plt.hist(testmeasurements,bins=23,normed=True)
#plt.plot(edges, Pvec(edges,0,1), linewidth=2, color='r')
print 'edges = ', edges
print 'Pvec(edges) = ', Pvec(edges)
#plt.plot(edgess,Pvec(edges,0,1),linewidth=2,color='r')
plt.show()

from scipy.stats import chi2
print chi2.sf(114, 98)

"""Malus law, yet another version"""

##Make an array which is the ratio of the two measurements,
##at each theta. Calculate your uncertainty in the
##ratio by finding the standard deviation of the 10 ratios
##2
##you took at a fixed angle in the last step of the procedure.
##Plot this ratio versus angle, with error bars. 

import numpy as np
from operator import truediv
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit
from scipy.stats import cosine
import math

# NOTE: All arrays skip the 90 deg data and goes to 92 degrees instead
xaxis = np.linspace(0,92,24)

"""Power through the polarizer: first element
in the pairs that occur in the array"""
Parray = sp.array([0.36, 0.45, 0.45, 0.53, 0.55, 0.56, 0.58, 0.58,
           0.59, 0.59, 0.58, 0.57, 0.56, 0.51, 0.47, 0.47,
           0.44, 0.41, 0.34, 0.30, 0.28, 0.21, 0.19, 0.12])

"""Power without the polarizer; second element
in the pairs that occur in the array"""
NoParray = sp.array([1.69, 1.69, 1.69, 1.87, 1.86, 1.86, 1.87, 1.87,
                     1.86, 1.86, 1.86, 1.87, 1.86, 1.86, 1.86, 1.86,
                     1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 1.86, 1.86])

RatioArray = map(truediv, Parray, NoParray)
print RatioArray

##Puncertainty = np.std(Parray)
##NoPuncertainty = np.std(NoParray)
##
##print Puncertainty, NoPuncertainty

Parray40 = sp.array([0.58, 0.60, 0.60, 0.60, 0.59])
NoParray40 = sp.array([1.86, 1.86, 1.86, 1.87, 1.86])
RatioArray40 = map(truediv, Parray40, NoParray40)

Puncertainty = sp.std(Parray40)     # Uncertainty of 0.008
NoPuncertainty = sp.std(NoParray40) # Uncertainty of 0.004
Ratiouncertainty = sp.std(RatioArray40)

print Puncertainty, NoPuncertainty, Ratiouncertainty

"""To plot this ratio vs angle with error bars, we first need to
define the error bars"""

# with polarizer
plt.errorbar(xaxis,Parray,yerr=Puncertainty,
         fmt='o',ecolor=None,elinewidth=None,
         capsize=3,barsabove=False,lolims=False,
         uplims=False,xlolims=False,xuplims=False)

plt.plot(xaxis,Parray)
plt.xlabel('Degrees')
plt.ylabel('Power')
plt.title('Input Angle vs. Power (With Polarizer)')
plt.show()

# without polarizer
plt.errorbar(xaxis,NoParray,yerr=NoPuncertainty,
         fmt='o',ecolor=None,elinewidth=None,
         capsize=3,barsabove=False,lolims=False,
         uplims=False,xlolims=False,xuplims=False)

plt.plot(xaxis,NoParray)
plt.xlabel('Degrees')
plt.ylabel('Power')
plt.title('Input Angle vs. Power (Without Polarizer)')
plt.show()

# Ratio
plt.errorbar(xaxis,RatioArray,yerr=Ratiouncertainty,
         fmt='o',ecolor=None,elinewidth=None,
         capsize=3,barsabove=False,lolims=False,
         uplims=False,xlolims=False,xuplims=False)

plt.plot(xaxis,RatioArray)
plt.xlabel('Degrees')
plt.ylabel('Power')
plt.title('Input Angle vs. Power (Ratio of [With Polarizer]/[Without Polarizer])')
plt.show()

"""Model function"""
def Pvec(thetaList = [], *args):
    for theta in thetaList:
        P_0 = np.math.radians(0.21301775147928995)
#    P_0 = np.math.radians(xaxis)
        theta_0 = 0
        C = 0
        return (P_0*((np.math.cos(theta - (theta_0)))*(np.math.cos(theta - (theta_0))))+C)

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

yunderlying = Pvec(xaxis)

yobserved=yunderlying+ynoise

testmeasurements=np.array(noisegenerator(sp.zeros(10000),'usebias'))

yobservedcorrected=yobserved-sp.mean(testmeasurements)

pguess = [12,5]
(popt, punc, rc, d) = curve_fit_custom(Pvec, xaxis, yobservedcorrected, pguess)
abestfitval=popt[0]
abestfitunc=punc[0]
bbestfitval=popt[1]
bbestfitunc=punc[1]
print "abest= %5.3f +/- %5.3f" % (abestfitval,abestfitunc)
print "bbest= %5.3f +/- %5.3f" % (bbestfitval,bbestfitunc)

xdense=sp.linspace(0,92,1000)

ymodelpoints=Pvec(xaxis,popt[0],popt[1])
ymodelsmooth=Pvec(xdense,popt[0],popt[1])

residuals=ymodelpoints-yobserved

smoothcurve=Pvec(xdense)

#plt.plot(xsamples,residuals,'bo')
plt.plot(xaxis, yobservedcorrected, 'go')
#plt.plot(xdense,ymodelsmooth,'g-')
#plt.plot(xdense,smoothcurve, 'b-')
#plt.xlim([-.1,20.1])
#plt.ylim([-.1,20.1])
plt.show()

print 'residuals=',residuals
print 'mean=               ', sp.mean(residuals)
print 'standard deviation= ', sp.sqrt(sp.var(residuals))

"""plt.autoscale(enable-True, axis=u 'both', tight=False)"""

plt.figure()
(events, edges, patches)=plt.hist(testmeasurements,bins=np.array([0,20,40,60,80,100]),normed=True)
#plt.plot(edges, Pvec(edges,0,1), linewidth=2, color='r')
plt.plot(edges,Pvec(edges),linewidth=2,color='r')
plt.show()

from scipy.stats import chi2
print chi2.sf(114, 98)

"""Discrete Fourier transform of a classical system"""

# PHYS4101 HW Set 2, Problem 2
# By Gaith Midani
# Oct. 20, 2015

import numpy as np
import matplotlib.pyplot as plt

w_0 = 10.0
gamma = 2.0
xi = 5.0
F_0 = 1.0
t = 10
t_range = np.arange(0,t,.1)
omega = np.sqrt(w_0**2-gamma**2)

def f(t):
    a = t
    f = a
    flist = []
    for t in t_range:
        a = F_0*np.exp(-xi*t)
        f += a
        flist.append(a)
    return flist

def g(t):
    a = t
    g = a
    glist = []
    for t in t_range:
        a = (1/omega)*np.exp(-gamma*t)*np.sin(omega*t)
        g += a
        glist.append(a)
    return glist

# Pt. c:
xlist = []
for t in t_range:
    x = F_0*np.exp(-t*(gamma+xi))*(np.exp(xi*t)*((xi-gamma)*np.sin(omega*t)-omega*np.cos(omega*t))+omega*np.exp(gamma*t))/(((gamma-xi)**2+omega**2)*omega)
    xlist.append(x)

plt.plot(t_range,xlist,label='x(t) vs. t')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.legend()
plt.show()

# Pt. d
x = np.convolve(f(t),g(t))
plt.xlabel('t')
plt.ylabel('x(t)=f*g')
plt.plot(x,label='Convolution of x(t)')
plt.xlim(0,10)
plt.show()

# Pt. e
f_fft = np.fft.fft(f(t))
g_fft = np.fft.fft(g(t))
xfft = np.fft.ifft(f_fft*g_fft)

plt.plot(xfft,label='Discrete Fourier Transform: IFT(FT(f)*FT(g)) vs. t')
plt.xlabel('t')
plt.ylabel('x(t)=IFT(FT(f)*FT(g))')
plt.show()

## Part e above has an identical result to part d's. Both differ from
## that of part c with respect to amplitude and frequency.

"""Prototype solution, attempted integral class"""

"""I'm being feidistic here. No plot, and hoping everything
else is right, for the sake of time management"""

# Homework Set #3, Problem 5
# Draft proto2

import numpy as np
import matplotlib.pyplot as plt

"""def f(r):
    E = -0.1
    m = 1.
    L = 1.
    k = 1.
    A = (L**2)/(2*m*r**2)
    V = (2*k)/(3*r**(3/2)) 
#    return 1/(np.sqrt(E-((L**2)/(2*m*r**2))-((2*k)/(3*r**(3/2)))))
    denom = np.sqrt(E-A-V)
    return denom**(-1)"""

r_0 = 1
"""Start r"""
r1 = 0.1
"""End r"""
r2 = 5*r_0
"""Steps"""
e = 250

def f(r):
    E = -0.1
    m = 1.
    L = 1.
    k = 1.
    A = (L**2)/(2*m*r**2)
    V = (2*k)/(3*r**(3/2))
    summed = E-A-V
    fxn = summed**(-1/2)
    return fxn

##def g(r):
##   V = (2*k)/(3*r**(3/2))
##   return V

def integral(startr,endr,numbrect):
    width = (float(endr)-float(startr))/numbrect
    runningSum = 0
    integrallist = []
    for i in np.arange(numbrect):
        height = f(startr + i*width)
        area = height*width
        runningSum += area
        integrallist.append(runningSum)
    return integrallist
    
#selectIntegral = integral(0.01,5,250)
#arraySelInt = np.array(selectIntegral)
#IntegralParam = np.arange(r1, r2,e)
#selectIntegral = integral(r1, r2, e)

#print selectIntegral
SpecInt = integral(r1,r2,e)
rmin = min(SpecInt)

print rmin
print SpecInt
"""Vaguely suspicious (but maybe not)
that the last point in the array is the minimum
BUT NOTE that this may very well be right, because the solution to integral
(i.e. singular answer, no array) determines the minimum"""

""""""
### Assuming rmin here is right, here's part c ###

"""Classical mechanical problem"""

## HS3, Problem 5
## Gaith Midani
## PHYS 4101

import numpy as np
import matplotlib.pyplot as plt
import sympy.mpmath

"""def f(r):
    E = -0.1
    m = 1
    L = 1
    k = 1
    sqrt1 = np.sqrt(2/m)
    sqrt2 = np.sqrt(E - (L**2)/(2*m*(float(r)**2)) + (2*k)/(3*(float(r)**(3/2))))
    fxn = sqrt1*sqrt2
    return fxn

print sympy.mpmath.findroot(f,2)"""


rmin = .66708

u = 1/rmin

def ufxn(theta):
    def r(theta):
        for theta in np.arange(0,7*np.pi):
            if theta == 0:
                r = rmin
            else:
                r = theta
        u = 1/r
    return u == np.sqrt(u)

plt.plot(ufxn(theta))
plt.plot()

