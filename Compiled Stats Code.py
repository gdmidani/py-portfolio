# -*- coding: utf-8 -*-
"""COMPILED STATISTICS CODE"""

"""Prints the floor, ceiling, and nearest rounded input of each element in
an inputted 1d array"""

import numpy as np

n = np.array(list(map(float,input().split())))

fl = np.floor(n)
ce = np.ceil(n)
ri = np.rint(n)

print(fl)
print(ce)
print(ri)

"""Sums over axis 0 and takes the product of that result (e.g., sum of
products) for an inputted 2d NxM array. The inputs contain space-separated 
values of N and M, and the next N lines contains M space-separated integers"""

import numpy as np

nm = list(map(int,(input().split())))
n = nm[0]
m = nm[1]

a = [(list(map(int,input().split()))) for i in range(m)]
soo = np.sum(a, axis=0)
proo = np.product(soo)
print(proo)

"""Computes the min over axis 1 and prints the max of that result, for an
inputted NxM array"""

import numpy as np

nm = list(map(int,input().split()))
n = nm[0]
m = nm[0]

a = [list(map(int,input().split())) for i in range(n)]

min1 = np.min(a,axis=1)
max1 = np.max(min1)
print(max1)

"""Calculates the mean (along axis 1), variance (along axis 0)
and standard deviation (along axis None, e.g., the flattened array) for an 
inputted 2d array of size NxM"""

import numpy as np

nm = list(map(int,input().split()))
n = nm[0]
m = nm[1]

a = [list(map(int,input().split())) for i in range(m)]

meana = np.mean(a,axis=1)
vara = np.var(a,axis=0)
stda = np.std(a)

print(meana)
print(vara)
print(stda)

"""Calculates Pearson correlation coefficient between inputted arrays x and y"""

n = input()
x = list(map(float,input().split()))
y = list(map(float,input().split()))

soox = sum(x)
mux = soox/int(n)

sooy = sum(y)
muy = sooy/int(n)

numx = []
numy = []
for i in range(int(n)):
    numx.append(x[i]-mux)
    numy.append(y[i]-muy)

numer = []
for i in range(int(n)):
    numer.append(numx[i]*numy[i])

numer1 = (sum(numer))

numxsq = []
numysq = []
for i in range(int(n)):
    numxsq.append((numx[i])**2)
    numysq.append((numy[i])**2)

sigmax = (sum(numxsq)/int(n))**(1/2)
sigmay = (sum(numysq)/int(n))**(1/2)

denom = int(n)*sigmax*sigmay

rho = (numer1/denom)
rhoround = round(rho,3)

print(rhoround)

"""Calculates Pearson correlation coefficient between array x and array y"""

n = 10
x = [15,  12,  8,   8,   7,   7,   7,   6,   5,   3]
y = [10,  25,  17,  11,  13,  17,  20,  13,  9,   15]

soox = sum(x)
mux = soox/int(n)

sooy = sum(y)
muy = sooy/int(n)

numx = []
numy = []
for i in range(int(n)):
    numx.append(x[i]-mux)
    numy.append(y[i]-muy)

numer = []
for i in range(int(n)):
    numer.append(numx[i]*numy[i])

numer1 = (sum(numer))

numxsq = []
numysq = []
for i in range(int(n)):
    numxsq.append((numx[i])**2)
    numysq.append((numy[i])**2)

sigmax = (sum(numxsq)/int(n))**(1/2)
sigmay = (sum(numysq)/int(n))**(1/2)

denom = int(n)*sigmax*sigmay

rho = (numer1/denom)
rhoround = round(rho,3)

print(rhoround)

"""Computes the slope of the regression line while treating the elements of
array x as the independent variable"""
# Regression line slope = Pearson Correlation Coefficient * (sigma_y)/(sigma_x)... sigma is std dev
# Pearson Correlation Coefficient (rho):

n = 10
x = [15,  12,  8,   8,   7,   7,   7,   6,   5,   3]
y = [10,  25,  17,  11,  13,  17,  20,  13,  9,   15]

soox = sum(x)
mux = soox/int(n)

sooy = sum(y)
muy = sooy/int(n)

numx = []
numy = []
for i in range(int(n)):
    numx.append(x[i]-mux)
    numy.append(y[i]-muy)

numer = []
for i in range(int(n)):
    numer.append(numx[i]*numy[i])

numer1 = (sum(numer))

numxsq = []
numysq = []
for i in range(int(n)):
    numxsq.append((numx[i])**2)
    numysq.append((numy[i])**2)

sigmax = (sum(numxsq)/int(n))**(1/2)
sigmay = (sum(numysq)/int(n))**(1/2)

denom = int(n)*sigmax*sigmay

rho = (numer1/denom)
rhoround = round(rho,3)

# Regression line slope:

regsl = rho*sigmay/sigmax
regslrnd = round(regsl,3)
print(regslrnd)