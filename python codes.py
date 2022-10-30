# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 00:49:45 2022

@author: wangs
"""


'''
This file containes all of my python codes to run the
simulations and calculate some important features.

'''



import math
import matplotlib.pyplot as plt
import numpy as np
import random


'''
We will do the mircrospopic event-driven stochastic
simulation.
Tdt contains the refined time segments.The reason to apply
this is to obtain evenly spaced (rather than random) time
spots so that we can take the mean of many runs.
Here below we set up all the parameters
'''
T=[]
k=[]
Tdt=[]
kdt=[]

k1=0.8
k2=0.6
v0=k1/k2
K=k2/k1
n=3
v=(4**n)*v0
print('the total volume v is ' + str(v))
NAstar=7*(4**n)
NBstar=7*(4**n)
NAB=25
print('At the initial state, NA*,NB* and NAB are ' + str(NAstar)+ ' ' + str(NBstar) + ' ' + str(NAB))
print('The initial concentration of NAB is ' + str(NAB/v))
r1=k1/v
r2=k2
print('r1 and r2 are ' + str(r1) + ' ' + str(r2))
print()

tmax=7
t=0
T1=0
T2=0

ts=0
count=-1
dt=0.001
# how many time spots we take in a microscopic simulation
length_scatter=tmax/dt
# runtime is how many simulations we decide to take in all
countrun=0    
runtime=150

NABsum=[0]*int(length_scatter)


'''
micro deterministic (not event driven) 
There is a problem to this algorithm. When scale is large, one needs to set smaller constant for initial 
P so that python won't go wrong
'''

m=min(NAstar,NBstar)
P=[0.1]*(m+1)
for n in range(1,m+1):
    P[n]=P[n-1]*(r1/r2)*((NAstar-n+1)*(NBstar-n+1))/n

sum0=0
for ele in P:
    sum0=sum0+ele

for i in range(m+1):
    P[i]=P[i]/sum0

EAB=0
for i in range(m+1):
    EAB+=i*P[i]

print('Estimated micro equilibrium EAB is ' + str(EAB))

'''
Macroscopic equilibrium
'''

def calculate(a, b, c):
    discriminant = b ** 2 - 4 * a * c
    if discriminant < 0:
        return None, None
    x1 = (-b - math.sqrt(discriminant)) / 2 * a
    x2 = (-b + math.sqrt(discriminant)) / 2 * a
    return x1, x2

NAini_concentrate = (NAstar)/v
NBini_concentrate = (NBstar)/v
a = float(1)
b = float(-(NAini_concentrate + NBini_concentrate + K))
c = float((NAini_concentrate)*(NBini_concentrate))
x1, x2 = calculate(a, b, c)
print('The estimated macro equilibrium of concentration is ' + str(x1))
if x1<NAB/v:
    print('AB has a decreasing tendency')
print('x1 and NAB concentration are ' + str(x1) + ' ' + str(NAB/v))


'''
Microscopic event driven simulation
'''
#choosse_count is where we select a time spot to calculate mean AB
choose_count=5500
kcollect=[]
#save the initial NAB that we set
NABini=NAB
while countrun<runtime:
    NAB=NABini
    NA = NAstar-NAB
    NB = NBstar-NAB
    t=0
    ts=0
    count=-1
    Tdt=[0]*int(length_scatter)
    kdt=[0]*int(length_scatter)
    while t<tmax:
        NAB_old=NAB
        if NB==0:
            T2=(-math.log(random.random()))/(r2*NAB)
            NA=NA+1
            NB=NB+1
            NAB=NAB-1
            t=t+T2
        elif NA==0:
            T2=(-math.log(random.random()))/(r2*NAB)
            NA=NA+1
            NB=NB+1
            NAB=NAB-1
            t=t+T2
        elif NAB==0:
            T1=(-math.log(random.random()))/(r1*NA*NB)
            NA=NA-1
            NB=NB-1
            NAB=NAB+1
            t=t+T2
        else:
            T1=(-math.log(random.random()))/(r1*NA*NB)
            T2=(-math.log(random.random()))/(r2*NAB)
            if T1<T2:
                NA=NA-1
                NB=NB-1
                NAB=NAB+1
                t=t+T1
            else:
                NA=NA+1
                NB=NB+1
                NAB=NAB-1
                t=t+T2 
        T.append(t)
        k.append(NAB)
            
        while ts<t and ts<tmax and count<=int(length_scatter)-2:
            count+=1
            ts=ts+dt
            Tdt[count]=ts
            kdt[count]=NAB_old
            #pick a time spot to collect data to calculate the std and mean
            if count==choose_count:
                kcollect.append(kdt[count]/v)
            
    for i in range(int(length_scatter)):
        NABsum[i]=NABsum[i]+kdt[i]

    countrun+=1

NAB_ave=[]
NAB_concen_ave=[]
for ele in NABsum:
    NAB_concen_ave.append(ele/runtime/v)
    NAB_ave.append(ele/runtime)
    

def std(a):
    meanval=x1
    stdevi=math.sqrt(sum((i-meanval)**2 for i in a)/(runtime-1))    
    return meanval, stdevi
print()
print('The macroscopic mean concentration at equilibrium'  + ' is ' + str(std(kcollect)[0]))
print('The standard deviation at '  + str(choose_count*dt) + ' sec is '+ str(std(kcollect)[1]))
print('the mean of trials at ' + str(choose_count*dt) + ' is ' + str(sum(i*v/runtime for i in kcollect)) )
print()






'''
Finally, we do the plotting.
There are three graphs in all.
'''

t = np.array(T)
k = np.array(k)

plt.scatter(t, k, label="superposition of {} trials".format(runtime))
plt.legend()
plt.show()

x4 = np.array(Tdt)
y4 = np.array(NAB_ave)

plt.scatter(x4, y4, label="Average of {} trials".format(runtime))
plt.legend()
plt.show()

x = np.array(Tdt)
y = np.array(kdt)

plt.scatter(x, y)
plt.show()


c=x1
a=x2*(x1-NABini/v)/(x2-NABini/v)
b=(x2-x1)*k1
d=(x1-NABini/v)/(x2-NABini/v)

x3 = np.linspace(0,tmax,4000)
#Here, we apply the macroscopic formula for [AB] covered in lecture to see if it coincide with average of even driven simulation
y3 = c + (c*d-a)*np.exp(-b*x3)/(1-d*np.exp(-b*x3))

x2 = np.array(Tdt)
y2 = np.array(NAB_concen_ave)

plt.plot(x3, y3, 'r', label="Predicted [AB](t)")

plt.scatter(x2, y2, label="Average of stochastic simulation [AB](t)")
plt.legend()
plt.show()





















