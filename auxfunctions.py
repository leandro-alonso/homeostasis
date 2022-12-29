#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 14:56:04 2019

@author: leandro
"""
from pylab import *


def getCAcurrents(cellsol):
    s=cellsol
    Temp=10
    V = s[:,0]
    IntCa = s[:,12]
    CaTM = s[:,3]
    CaTH = s[:,4]
    CaSM = s[:,5]
    CaSH = s[:,6]

    gCaT = s[:,14]  # Low Threshold Calcium Maximal Conductance
    gCaS = s[:,15]  # Slow Calcium Maximal Conductance

    CaRev = CaNernst(IntCa, Temp)
    ECaT = CaRev  # Low Threshold Calcium Reversal Potential
    ECaS = CaRev  # Slow Calcium Reversal Potential
    area=1
    iCaT = iIonic(gCaT, CaTM, CaTH, 3, V, ECaT)*area
    iCaS = iIonic(gCaS, CaSM, CaSH, 3, V, ECaS)*area
    # Calcium Current Sensor Dependence
    ICa = iCaT + iCaS
    return ICa

def getRandomInitStateLiuCubic():
    # Initial Conditions
    V0 = -50.0
    IntCa0 = 0.05
    # Initial Maximal Conductances
    gNa0 = rand()
    gCaT0 = rand()
    gCaS0 = rand()
    gH0 = rand()
    gKd0 = rand()
    gKCa0 = rand()
    gA0 = rand()
    
    initV=-50
    initstate = [initV, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0.05, gNa0, gCaT0, gCaS0, gH0, gKd0, gKCa0, gA0, 0.0, 1.0, 0.0, 1.0, 0.0]#, 0, 0 ,0]
    return initstate


def getRandomInitStateLiuAdaptive():
    # Initial Conditions
    V0 = 0.0
    IntCa0 = 0.05
    # Initial Maximal Conductances
    gNa0 = rand()
    gCaT0 = rand()
    gCaS0 = rand()
    gH0 = rand()
    gKd0 = rand()
    gKCa0 = rand()
    gA0 = rand()
    
    initV=-50
    initstate = [initV, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0.05, gNa0, gCaT0, gCaS0, gH0, gKd0, gKCa0, gA0, 0.0, 1.0, 0.0, 1.0, 0.0]#, 0, 0 ,0]
    return initstate


def getRandomInitStateLiuAdaptive_multiple_readouts():
    # Initial Conditions
    V0 = 0.0
    IntCa0 = 0.05
    # Initial Maximal Conductances
    gNa0 = rand()
    gCaT0 = rand()
    gCaS0 = rand()
    gH0 = rand()
    gKd0 = rand()
    gKCa0 = rand()
    gA0 = rand()
    
    initV=-50
    initstate = [initV, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0.05, gNa0, gCaT0, gCaS0, gH0, gKd0, gKCa0, gA0, 0.0, 1.0, 0.0, 1.0, 0.0,0,1,0,1]#, 0, 0 ,0]
    return initstate


def getRandomInitStateOleary():
    # Initial Conditions
    V0 = -50.0
    IntCa0 = 0.05
    # Initial Maximal Conductances
    gNa0 = rand()
    gCaT0 = rand()
    gCaS0 = rand()
    gH0 = rand()
    gKd0 = rand()
    gKCa0 = rand()
    gA0 = rand()
    
    initV=-50
    initstate = [initV, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0.05, gNa0, gCaT0, gCaS0, gH0, gKd0, gKCa0, gA0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]#, 0, 0 ,0]
    return initstate




    
# Synaptic Functions
# Synaptic Current
# 'Vpost' is the postsynaptic potential
# 'sact' is the synaptic activation
def iSyn(gsyn, sact, Vpost, Esyn):
    isyn = gsyn*sact*(Vpost - Esyn)
    return isyn


# Steady State Synaptic Activation
# 'Vth' is the synaptic threshold
# 'Vpre' is the presynaptic potential
# 'Delta' is the half activation potential
def synSS(Vth, Vpre, Delta):
    synSS = 1/(1 + np.exp((Vth - Vpre)/Delta))
    return synSS


# Synaptic Activation Time Constant
# 'k' is the rate constant for receptor-transmitter dissociation
def tauSyn(synSS, k):
    tauSyn = 1  + (1 - synSS)/k
    return tauSyn

def tauSynBound(synSS, k,bound):
    tauSyn = bound  + (1 - synSS)/k
    return tauSyn


# Synaptic Constants
# Glutamatergic reversal potential
Eglut = -70.0
# Cholinergic reversal potential
Echol = -80.0
# K for glutamateric synapses
kglut = 1./40.
# K for cholinergic synapses
kchol = 1./100.
# Threshold for all synapse types
Vthresh = -35.0
# Delta for all synapse types
DSyn = 5.00


def getParameters():
    # Temperature of Simulation
    Temp = 10.0  # degrees Celsius
    
    # Initialize Parameters
    p = [1.,                 # 0 Capacitance
         30.0,                # 1 ENa
         80.0,                # 2 ECaT
         80.0,                # 3 ECaS
         -20.0,               # 4 EH
         -80.0,               # 5 EK
         -80.0,               # 6 EKCa
         -80.0,               # 7 EA
         -50.0,               # 8 EL
    # =============================================================================
    #      1.0,               # 9 gNa
    #      2.00,                # 10 gCaT
    #      4.00,                # 11 gCaS
    #      0.00,                # 12 gH
    #      1.0,                # 13 gKd
    #      1.0,                # 14 gKCa
    #      1.0,                # 15 gA
    #      0.01,                # 16 gL
    #      0.94,                 # 17 F
    #      1.0,               # 9 gNa
    #      2.00,                # 10 gCaT
    #      4.00,                # 11 gCaS
    #      0.00,                # 12 gH
    #      70.0,                # 13 gKd
    #      40.0,                # 14 gKCa
    #      50.0,                # 15 gA
    # =============================================================================
         0.01,                # 9 gL
         0.94,                 # 10 F
         0.05,                # 11 Ca0 (uM)
         20.0,                # 12 tauIntCa
    # =============================================================================
    #      0.1,                 # 20 FNa
    #      0.1,                 # 21 SCaT
    #      0.1,                 # 22 SCaS
    #      0.1,                 # 23 FKd
    #      0.1,                 # 24 SKd
    #      0.1,                 # 25 SKCa
    #      0.1,                 # 26 DKCa
    #      0.1,                 # 27 SA
    #      0.1,                 # 28 DA
    #      0.1,                 # 29 SH
    #      0.1,                 # 30 DH
    # =============================================================================
         0.0]                 # 13 Iapp (nA)
    return p


# Define the steady state activation/inactivation functions
def boltzSS(Volt, A, B):
    act = 1/(1 + np.exp((Volt + A)/B))
    return act


# Define the time constant function
def tauX(Volt, CT, DT, AT, BT):
    timeconst = CT - DT/(1 + np.exp((Volt + AT)/BT))
    return timeconst


def differenceWithGap(x, y, gap): 
    if(abs(x-y)<gap):
        return 0
    else:
        return x-y


# Some Models require a special time constant function
def spectau(Volt, CT, DT, AT, BT, AT2, BT2):
    spec = CT + DT/(np.exp((Volt + AT)/BT) + np.exp((Volt + AT2)/BT2))
    return spec


# Define the ionic currents
# q is the exponent of the activation variable m
def iIonic(g, m, h, q, Volt, Erev):
    flux = g*pow(m, q)*h*(Volt - Erev)
    return flux


# Define the regulation of the maximal conductances
def GEvol(A, Fbar, F, B, Sbar, S, C, Dbar, D, g):
    tauG = 5000.0 * 2
    gdot = ((A*(Fbar - F) + B*(Sbar - S) + C*(Dbar - D))*g)/tauG
    return gdot


def altGEvol(A, Fbar1,Fbar2, F, B, Sbar1,Sbar2, S, C, Dbar1,Dbar2, D, g):
    tauG = 5000.0 * 2
    gdot = (((A*(Fbar1 - F) + B*(Sbar1 - S) + C*(Dbar1 - D))*(A*(Fbar2 - F) + B*(Sbar2 - S) + C*(Dbar2 - D)))*g)/tauG
    return gdot

def altGEvolBis(A, Fbar1, F, B, Sbar1, S, C, Dbar1, D, g):
    tauG = 1000.0 * 2
    gdot = ((A*(Fbar1 - F)**3 + B*(Sbar1 - S)**3 + C*(Dbar1 - D)**3)*g)/tauG
    return gdot

def altGEvol3(A, Fbar1, F, B, Sbar1, S, C, Dbar1, D, g):
    tauG = 5000.0 * 2
    gdot = ((A* differenceWithGap(Fbar1, F, 0.01) + B*differenceWithGap(Sbar1, S, 0.01) + C*differenceWithGap(Dbar1, D, 0.01))*g)/tauG
    return gdot

# =============================================================================
# def altGEvol4(A, mrnaF, B, mrnaS, C, mrnaD, g):
#     tauG = 1000.0 * 1
#     gdot = ((A* mrnaF + B*mrnaS + C*mrnaD)*g)/tauG
#     return gdot
# =============================================================================

def altGEvol4(A, mrnaF, B, mrnaS, C, mrnaD, g):
    tauG = 1000.0 * 1
    gdot = (((A* mrnaF + B*mrnaS + C*mrnaD)*g)-(10**(-5))*g**3)/tauG
    return gdot

def altGEvol5(A, Fbar, F, B, Sbar, S, C, Dbar, D, g):
    tauG = 1000.0 * 1
    gdot = (((A*(Fbar - F) + B*(Sbar - S) + C*(Dbar - D))*g)-(10**(-6))*g**3)/tauG
    return gdot

def altGEvol6(A, Fbar, F, B, Sbar, S, C, Dbar, D, g, gamma):
    tauG = 1000.0 * 1
    gdot = (((A*(Fbar - F) + B*(Sbar - S) + C*(Dbar - D))*g)-(gamma)*g**3)/tauG
    return gdot


def altGEvol7(A, Fbar, F, B, Sbar, S, C, Dbar, D, g, gamma, tauG):
#    tauG = 1000.0 * 1
    gdot = (((A*(Fbar - F) + B*(Sbar - S) + C*(Dbar - D))*g)-(gamma)*g**3)/tauG
    return gdot


# Update the sensor values
def Sensor(G, M, H):
    sensor = G*pow(M, 2)*H
    return sensor


# Calcium current dependence of sensor activation
# Z is the sensor threshold
def SigM(Z, ICa):
    M = 1/(1 + np.exp(Z + ICa))
    return M


# Calcium current dependence of sensor inactivation
def SigH(Z, ICa):
    H = 1/(1 + np.exp(-Z - ICa))
    return H


# Define the concentration dependent Ca reversal potential
def CaNernst(CaIn, temp):
    R = 8.314*pow(10, 3)  # Ideal Gas Constant (*10^3 to put into mV)
    T = 273.15 + temp  # Temperature in Kelvin
    z = 2.0  # Valence of Calcium Ions
    Far = 96485.33  # Faraday's Constant
    CaOut = 3000.0  # Outer Ca Concentration (uM)
    CalRev = ((R*T)/(z*Far))*np.log(CaOut/CaIn)
    #print 'calrev ', CalRev
    return CalRev




def getSensorsActivity(sol):
    MF = sol[:,20]
    HF = sol[:,21]
    MS = sol[:,22]
    HS = sol[:,23]
    MD = sol[:,24]
    Gf=10
    Gs=3
    Gd=1
    F = Gf * MF*MF *HF
    S = Gs * MS*MS *HS
    D = Gd * MD    *MD
    return [F,S,D]