#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 15:51:20 2019

@author: leandro
"""

# Preliminaries and Modules

import os
import numpy as np
from scipy import *
from scipy.integrate import odeint
from pylab import *
import math
import sys
import time as ttime

from singlecell_liu_bound_dynTauG import *
from getCurrents_singlecell_liu_bound_dynTauG import *
# from getRepertoires import *
import string
import random
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def getCurrentscapeData_EK(pathto_initialcondition, krev_highK, boundval, maxtauG, targets):
    pathtocontinue = pathto_initialcondition
    y0= load(pathtocontinue, allow_pickle=True, encoding = 'latin1')    
    
    krev=-80

    dt=0.1
    p=getParameters()    
    spread=0.001

    juicethres=0.075 # THIS IS DYNAMIC
    # juicethres=-1 # THIS IS OFF
    # juicethres=10 # THIS IS ALWAYS ON

    nsecs=5
    nchunks=1
    timestart = ttime.time()
    timeend = ttime.time()
    krev = krev_highK
    sol=[]
    for j in range(nsecs * 10000):
        sol.append(y0)
        y0 = ConductanceBased(y0,0, p,targets,maxtauG,boundval,krev,juicethres,dt,)
    currents = getCurrents_ConductanceBased(array(sol),0, p,targets,maxtauG,boundval,krev,juicethres,dt,)
    voltage = array(sol)[:,0]
    return [voltage, currents]


def getCurrentscapeData_injCurr(pathto_initialcondition, highK, boundval, maxtauG, targets, injCurrent):
    pathtocontinue = pathto_initialcondition
    y0= load(pathtocontinue, allow_pickle=True, encoding = 'latin1')    
    
    krev=-80
    

    dt=0.1
    p=getParameters()
    p[13]=injCurrent

    spread=0.001
    juicethres=0.075 # THIS IS DYNAMIC
    # juicethres=-1 # THIS IS OFF
    # juicethres=10 # THIS IS ALWAYS ON

    nsecs=5
    nchunks=1
    timestart = ttime.time()
    timeend = ttime.time()

    leakval_control=0.01
    rev_leak_control=-50
    rev_leak_high_K= -40

    p[8]=rev_leak_control

    pleak=array(p).copy()
    krev_control=-80
    krev_highK=-40
    krev=krev_control

    
    sol=[]
    for j in range(nsecs * 10000):
        sol.append(y0)
        y0 = ConductanceBased(y0,0, p,targets,maxtauG,boundval,krev,juicethres,dt,)
    currents = getCurrents_ConductanceBased(array(sol),0, p,targets,maxtauG,boundval,krev,juicethres,dt,)
    currents=list(currents)
    injected_current= ones(len(currents[0]))*injCurrent
    print(injected_current)

    currents.append(injected_current)
    voltage = array(sol)[:,0]
    return [voltage, array(currents)]



def getCurrentscapeData_removeCurrent(pathto_initialcondition, control, boundval, maxtauG, targets, whichchannel):
    pathtocontinue = pathto_initialcondition
    y0= load(pathtocontinue, allow_pickle=True, encoding = 'latin1')    
    # y0[whichchannel]=0

    krev=-80
    dt=0.1
    p=getParameters()    
    spread=0.001

    juicethres=0.075 # THIS IS DYNAMIC
    # juicethres=-1 # THIS IS OFF
    # juicethres=10 # THIS IS ALWAYS ON

    nsecs=5
    nchunks=1
    timestart = ttime.time()
    timeend = ttime.time()

    
    sol=[]
    for j in range(nsecs * 10000):
        if(control==False):
            y0[whichchannel]=0
        sol.append(y0)
        y0 = ConductanceBased(y0,0, p,targets,maxtauG,boundval,krev,juicethres,dt,)
    currents = getCurrents_ConductanceBased(array(sol),0, p,targets,maxtauG,boundval,krev,juicethres,dt,)
    voltage = array(sol)[:,0]
    return [voltage, currents]




def getCurrentscapeData_changeLeak(pathto_initialcondition, boundval, maxtauG, targets, newleakval):
    pathtocontinue = pathto_initialcondition
    y0= load(pathtocontinue, allow_pickle=True, encoding = 'latin1')    

    krev=-80
    dt=0.1
    p=getParameters()    
    p[9]=newleakval
    spread=0.001

    juicethres=0.075 # THIS IS DYNAMIC
    # juicethres=-1 # THIS IS OFF
    # juicethres=10 # THIS IS ALWAYS ON

    nsecs=5
    nchunks=1
    timestart = ttime.time()
    timeend = ttime.time()

    # leakval_control=0.01
    # rev_leak_control=-50
    # rev_leak_high_K= -40
    sol=[]
    for j in range(nsecs * 10000):
        # y0[whichchannel]=0
        sol.append(y0)
        y0 = ConductanceBased(y0,0, p,targets,maxtauG,boundval,krev,juicethres,dt,)
    currents = getCurrents_ConductanceBased(array(sol),0, p,targets,maxtauG,boundval,krev,juicethres,dt,)
    voltage = array(sol)[:,0]
    return [voltage, currents]

