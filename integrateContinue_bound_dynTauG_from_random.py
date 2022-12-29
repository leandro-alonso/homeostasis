#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Dec. 2022

@author: Leandro Alonso
"""

# Preliminaries and Modules

import os
import numpy as np
from scipy import *
from pylab import *
import math
import sys
import time as ttime
from singlecell_liu_bound_dynTauG import *
import string
import random

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))



#value of gamma 10^-5
boundval=pow(10,-5)*1
#create a random name for the simulation
hashname=id_generator()
name=hashname+'-liu-cubic'+str(boundval)+'.schemeLIU'
# initial condition is random
y0= getRandomInitStateLiuCubic()
y0=list(y0)
# initial condition for deltas, can't be zero if starting from random
y0.extend([1.,1.,1.])
y0.extend([0.0])


dt=0.1
p=getParameters()
#conductance evolution timescale in miliseconds
maxtauG=2000
#parameter for the gaussian functions used in Sf
spread=0.001
#reversal potential of potassium
krev=-80

#this is \rho in eq. (5). 
juicethres=0.075 # THIS IS DYNAMIC. Homeostatic mechanism will be gated on or off
# juicethres=-1 # THIS IS OFF. Homeostatic mechanism is off. conductances won't evolve
# juicethres=10 # THIS IS ALWAYS ON. Homeostatic mechanism never turns off. 


#parameters of the simulation. Two chunks of 30 seconds per file. Each file will contain 1 minute of data
nsecs=30
nchunks=2
timestart = ttime.time()
timeend = ttime.time()
num_minutes=10
totalsecsprevio=0

#sensors targets
targets=[0.10, 0.07, 0.07] #these are the targets for the larger population March 28, 2022


strtargets=str(targets[0])+'-'+str(targets[1])+'-'+str(targets[2])
# pathtostore= '/Users/leandro/4.0.local.workspace/results/homeostasis_2022/move_target_F_'+strtargets+'/init-fullstate-'+name+'-jthres-'+str(juicethres)+'-targets-'+strtargets+'/'
#create a folder to store the simulation data
pathtostore= './simulations/init-fullstate-'+name+'/'
if not os.path.exists(pathtostore):os.makedirs(pathtostore)


#main loop
for min_num in range(num_minutes):
    trace=[]
    fulltrace=[]
    gsevol=[]
    states=[]

    Ftg=[]
    Stg=[]
    Dtg=[]

    error_t=[]
    schemes=[]

    avFt=[]
    avSt=[]
    avDt=[]
    
    alphat=[]
    krevst=[]
    

    for i in range(nchunks):
        timestart = ttime.time()
        sol=[]
        for j in range(nsecs * 10000):
            sol.append(y0)
            y0 = ConductanceBased(y0,0, p,targets,maxtauG,boundval,krev, juicethres,dt,)

        # subsample solutions to save space
        sol=array(sol)[::5,:]
        nexttrace=sol[0,:]
        trace.extend(sol[:,0])
        gs = sol[::100,13:20]
        gsevol.extend(gs)
        states.append(y0)
        taug=int(1000/(0.5))
        F,S,D = getSensorsActivity(sol)
        cvF=convolve(F,ones(taug)/taug, mode = 'valid')
        cvD=convolve(D,ones(taug)/taug, mode = 'valid')
        cvS=convolve(S,ones(taug)/taug, mode = 'valid')
        Ftg.extend(cvF[::100])
        Stg.extend(cvS[::100])
        Dtg.extend(cvD[::100])

        avFt.extend(sol[:,25])
        avSt.extend(sol[:,26])
        avDt.extend(sol[:,27])
        alphat.extend(sol[:,28])

        krevs= ones(len(sol[::40,28:-2]))*krev
        krevst.extend(krevs)
        
        timeend = ttime.time()
        print(hashname, ': ' + str(i) + ' done chunk in: ' + str(timeend-timestart))
    
    sampledt=0.5
    totalsecs=int(len(trace)/(1000/sampledt))
    argsevol=array(gsevol)

    #save chunk and repeat
    fout = pathtostore+'run.'+name+'.nsecs.'+str(totalsecs+totalsecsprevio)+'.'+strtargets
    save(fout+'.data.npy', [trace, argsevol, states, Ftg, Dtg,Stg,avFt,avSt,avDt,alphat,krevst])
    save(fout+'.laststate.npy', y0)
    pathtocontinue = fout+'.laststate.npy'
    totalsecsprevio=totalsecs+totalsecsprevio



