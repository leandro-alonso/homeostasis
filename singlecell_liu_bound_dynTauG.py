#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Dec 2022

@author: Leandro Alonso
"""

from pylab import *
from auxfunctions import *

def getLiuScheme():
    A1=[
    [1.,0.,0.],
    [0.,1.,0.],
    [0.,1.,0.],
    [0.,1.,1.],
    [1.,-1.,0.],
    [0.,-1.,-1.],
    [0.,-1.,-1.]
    ]
    return array(A1)


def ConductanceBased(s,t, p,targets,maxTauG,boundval,krev, juicethres,dt,):

    # Fixed parameters
    Temp=10
    
    # State Variables and Parameters
    # Voltage and Gating Variables
    V =     s[0]
    NaM =   s[1]
    NaH =   s[2]
    CaTM =  s[3]
    CaTH =  s[4]
    CaSM =  s[5]
    CaSH =  s[6]
    HM =    s[7]
    KdM =   s[8]
    KCaM =  s[9]
    AM =    s[10]
    AH =    s[11]
    IntCa = s[12]
    
    # Dynamic Maximal Conductances
    gNa =   s[13]      # Transient Sodium Maximal Conductance
    gCaT =  s[14]      # Low Threshold Calcium Maximal Conductance
    gCaS =  s[15]      # Slow Calcium Maximal Conductance
    gH =    s[16]      # Hyperpolarization Activated Cation Maximal Conductance
    gKd =   s[17]      # Potassium Maximal Conductance
    gKCa =  s[18]      # Calcium Dependent Potassium Maximal Conductance
    gA =    s[19]      # Transient Potassium Maximal Conductance
    
    # Activation/Inactivation of Sensors
    MF = s[20]
    HF = s[21]
    MS = s[22]
    HS = s[23]
    MD = s[24]
    
    # average value of sensors 
    avF = s[25] 
    avS = s[26]
    avD = s[27]

    alpha = s[28]
    conductances=array(s[13:20])

    # Update Ca Reversal Potential
    CaRev = CaNernst(IntCa, Temp)
    
    # ====================
    # Parameters
    # ====================
    C = p[0]  # Capacitance (uF / cm^2)

    # Reversal Potentials (mV)
    ENa =       p[1]  # Sodium Reversal Potential
    ECaT = CaRev  # Low Threshold Calcium Reversal Potential
    ECaS = CaRev  # Slow Calcium Reversal Potential
    EH =        p[4]  # Hyperpolarization Activated Cation Reversal Potential
    
    EKd = krev #p[5]  # Potassium Reversal Potential
    EKCa =krev # p[6]  # Calcium Dependent Potassium Reversal Potential
    EA = krev #p[7]  # Transient Potassium Reversal Potential    
    EL =        p[8]  # Leak Reversal Potential
    # Fixed Maximal Conductances
    gL =        p[9]  # Leak Maximal Conductance    
    # Intracellular Calcium Constants
    caF =0.94# p[10]  # Current (nA) to Concentration (uM) Conversion Factor (uM / nA)
    #ESTE HARDCODEO ES PROBLEMATICO.
#    caF =9.4# p[10]  # Current (nA) to Concentration (uM) Conversion Factor (uM / nA)
    Ca0 =       p[11]  # Background Intracellular Calcium Concentration (uM)
    tauIntCa =  p[12]  # Calcium buffer time constant (ms)
    # Applied Current
    Iapp =      p[13]
    
    # ====================
    # Equations
    # ====================

    # Sensor Update    
    Gf=10
    Gs=3
    Gd=1
    F = Gf * MF*MF *HF
    S = Gs * MS*MS *HS
    D = Gd * MD    *MD
    
    #F = Sensor(10.0, MF, HF)  # Fast sensor
    #S = Sensor(3.0, MS, HS)   # Slow sensor
    #D = Sensor(1.0, MD, 1.0)  # DC sensor
    

    # Steady State Gating Variables
    NaMinf = boltzSS(V, 25.5, -5.29)  # m^3
    NaHinf = boltzSS(V, 48.9, 5.18)  # h
    CaTMinf = boltzSS(V, 27.1, -7.20)  # m^3
    CaTHinf = boltzSS(V, 32.1, 5.50)  # h
    CaSMinf = boltzSS(V, 33.0, -8.1)  # m^3
    CaSHinf = boltzSS(V, 60.0, 6.20)  # h
    HMinf = boltzSS(V, 70.0, 6.0)  # m
    KdMinf = boltzSS(V, 12.3, -11.8)                         # m^4
    KCaMinf = (IntCa/(IntCa + 3.0))*boltzSS(V, 28.3, -12.6)  # m^4
    AMinf = boltzSS(V, 27.2, -8.70)                          # m^3
    AHinf = boltzSS(V, 56.9, 4.90)                           # h
    
    # Time Constants (ms)
    tauNaM = tauX(V, 1.32, 1.26, 120.0, -25.0)
    tauNaH = tauX(V, 0.0, -0.67, 62.9, -10.0)*tauX(V, 1.50, -1.00, 34.9, 3.60)
    tauCaTM = tauX(V, 21.7, 21.3, 68.1, -20.5)
    tauCaTH = tauX(V, 105.0, 89.8, 55.0, -16.9)
    tauCaSM = spectau(V, 1.40, 7.00, 27.0, 10.0, 70.0, -13.0)
    tauCaSH = spectau(V, 60.0, 150.0, 55.0, 9.00, 65.0, -16.0)
    tauHM = tauX(V, 272.0, -1499.0, 42.2, -8.73)
    tauKdM = tauX(V, 7.20, 6.40, 28.3, -19.2)
    tauKCaM = tauX(V, 90.3, 75.1, 46.0, -22.7)
    tauAM = tauX(V, 11.6, 10.4, 32.9, -15.2)
    tauAH = tauX(V, 38.6, 29.2, 38.9, -26.5)
    
    # Ionic Currents (mV / ms)
    area=1.#0.628
    iNa = iIonic(gNa, NaM, NaH, 3, V, ENa)
    iCaT = iIonic(gCaT, CaTM, CaTH, 3, V, ECaT)
    iCaS = iIonic(gCaS, CaSM, CaSH, 3, V, ECaS)
    iH = iIonic(gH,  HM, 1, 1, V, EH)
    iKd = iIonic(gKd, KdM, 1, 4, V, EKd)
    iKCa = iIonic(gKCa, KCaM, 1, 4, V, EKCa)
    iA = iIonic(gA, AM, AH, 3, V, EA)
    iL = iIonic(gL, 1, 1, 1, V, EL)
    # Calcium Current Sensor Dependence
    ICa = iCaT + iCaS


#   Evolution of V using Exponential Euler
    giNa = get_g(gNa, NaM, NaH, 3, V, ENa)
    giCaT = get_g(gCaT, CaTM, CaTH, 3, V, ECaT)
    giCaS = get_g(gCaS, CaSM, CaSH, 3, V, ECaS)
    giH = get_g(gH,  HM, 1, 1, V, EH)
    giKd = get_g(gKd, KdM, 1, 4, V, EKd)
    giKCa = get_g(gKCa, KCaM, 1, 4, V, EKCa)
    giA = get_g(gA, AM, AH, 3, V, EA)
    giL = get_g(gL, 1, 1, 1, V, EL)


    EiNa  = giNa * ENa   
    EiCaT = giCaT* ECaT
    EiCaS = giCaS* ECaS
    EiH   = giH  * EH
    EiKd  = giKd * EKd
    EiKCa = giKCa* EKCa
    EiA   = giA  * EA
    EiL   = giL  * EL

    gE = (EiNa + EiCaT + EiCaS+ EiH+ EiKd+ EiKCa+ EiA+ EiL)
    g = giNa + giCaT + giCaS+ giH+ giKd+ giKCa+ giA+ giL
    vinf = (gE + (Iapp/1.))/g
    tauV = C / g
    #print(tauV,vinf,g,gE, ECaT)
    dV = vinf + (V -vinf) * exp(-dt / tauV)

#ALTERNATIVE SENSORS
    MbarF = SigM(14.8, ICa)  # Fast sensor activation
    HbarF = SigH(9.8, ICa)  # Fast sensor inactivation
    MbarS = SigM(7.20, ICa)  # Slow sensor activation
    HbarS = SigH(2.80, ICa)  # Slow sensor inactivation
    MbarD = SigM(3.0, ICa)  # DC sensor activation
    
    
 
    
    # State Equations
    # Voltage Time Evolution: C*dV/dt = -I_ionic + I_applied; I_ionic = sum(g_i*m^q*h*(V - E_i))
#    iH=0
    #dV = (-(iNa + iCaT + iCaS + iH + iKd + iKCa + iA + iL) + Iapp)/C
    # Gating Variable Time Evolution: dX/dt = (X_inf - X)/tau_X
# =============================================================================
#     dNaM = (NaMinf - NaM)/tauNaM
#     dNaH = (NaHinf - NaH)/tauNaH
#     dCaTM = (CaTMinf - CaTM)/tauCaTM
#     dCaTH = (CaTHinf - CaTH)/tauCaTH
#     dCaSM = (CaSMinf - CaSM)/tauCaSM
#     dCaSH = (CaSHinf - CaSH)/tauCaSH
#     dHM = (HMinf - HM)/tauHM
#     dKdM = (KdMinf - KdM)/tauKdM
#     dKCaM = (KCaMinf - KCaM)/tauKCaM
#     dAM = (AMinf - AM)/tauAM
#     dAH = (AHinf - AH)/tauAH
# =============================================================================

#   Evolution of gating variables using Exponential Euler    
    dNaM   = NaMinf + (-NaMinf + NaM) * exp(-dt / tauNaM)
    dNaH   = NaHinf + (-NaHinf + NaH) * exp(-dt / tauNaH)
    dCaTM   = CaTMinf + (-CaTMinf + CaTM) * exp(-dt / tauCaTM)
    dCaTH   = CaTHinf + (-CaTHinf + CaTH) * exp(-dt / tauCaTH)
    dCaSM   = CaSMinf + (-CaSMinf + CaSM) * exp(-dt / tauCaSM)
    dCaSH   = CaSHinf + (-CaSHinf + CaSH) * exp(-dt / tauCaSH)
    dHM   = HMinf + (-HMinf + HM) * exp(-dt / tauHM)
    dKdM   = KdMinf + (-KdMinf + KdM) * exp(-dt / tauKdM)
    dKCaM   = KCaMinf + (-KCaMinf + KCaM) * exp(-dt / tauKCaM)
    dAM   = AMinf + (-AMinf + AM) * exp(-dt / tauAM)
    dAH   = AHinf + (-AHinf + AH) * exp(-dt / tauAH)

    # Conductance Time Evolution    
# =============================================================================
#     dgNa =  altGEvol7(1.0, Fbar1, F, 0.0, Sbar1, S, 0.0, Dbar1, D, gNa, gamma, tauG)
#     dgCaS = altGEvol7(0.0, Fbar1, F, 1.0, Sbar1, S, 0.0, Dbar1, D, gCaS, gamma, tauG)
#     dgCaT = altGEvol7(0.0, Fbar1, F, 1.0, Sbar1, S, 0.0, Dbar1, D, gCaT, gamma, tauG)
#     dgH =   altGEvol7(0.0, Fbar1, F, 1.0, Sbar1, S, 1.0, Dbar1, D, gH, gamma, tauG)
#     dgKd =  altGEvol7(1.0, Fbar1, F, -1.0, Sbar1, S, 0.0,Dbar1, D, gKd, gamma, tauG)
#     dgKCa = altGEvol7(0.0, Fbar1, F, -1.0, Sbar1, S, -1.0, Dbar1, D, gKCa, gamma, tauG)
#     dgA =   altGEvol7(0.0, Fbar1, F, -1.0, Sbar1, S, -1.0, Dbar1, D, gA, gamma, tauG)
# 
# =============================================================================

    # dgNa =  gNa + altGEvol7(1.0, Fbar1, F, 0.0, Sbar1, S, 0.0, Dbar1, D, gNa, boundval, tauG) * dt
    # dgCaS = gCaS + altGEvol7(0.0, Fbar1, F, 1.0, Sbar1, S, 0.0, Dbar1, D, gCaS, boundval, tauG)* dt
    # dgCaT = gCaT + altGEvol7(0.0, Fbar1, F, 1.0, Sbar1, S, 0.0, Dbar1, D, gCaT, boundval, tauG)* dt    
    # dgH = gH +  altGEvol7(0.0, Fbar1, F, 1.0, Sbar1, S, 1.0, Dbar1, D, gH, boundval, tauG)* dt
    # dgKd = gKd+ altGEvol7(1.0, Fbar1, F, -1.0, Sbar1, S, 0.0,Dbar1, D, gKd, boundval, tauG)* dt
    # dgKCa =gKCa+ altGEvol7(0.0, Fbar1, F, -1.0, Sbar1, S, -1.0, Dbar1, D, gKCa, boundval, tauG)* dt
    # dgA =gA+   altGEvol7(0.0, Fbar1, F, -1.0, Sbar1, S, -1.0, Dbar1, D, gA, boundval, tauG)* dt

    tauG=maxTauG

#   Conductance Time Evolution  --- Euler  

    Fbar = targets[0]
    Sbar = targets[1]
    Dbar = targets[2]

    A=getLiuScheme()
    delta = array([(Fbar-F),(Sbar-S),(Dbar-D)])
    R=dot(A, delta)
    
    # this is \gamma
    boundvals = ones(7)*boundval

    dconductances = (1./tauG) * (R * conductances - boundvals * conductances**3) * alpha
    conductances = conductances + (dconductances * dt)

#   Sensors Time Evolution  --- Exponential Euler  
    dMF = MbarF + (-MbarF + MF)*exp(-dt/0.5)
    dHF = HbarF + (-HbarF + HF)*exp(-dt/1.50)
    dMS = MbarS + (-MbarS + MS)*exp(-dt/50.0)
    dHS = HbarS + (-HbarS + HS)*exp(-dt/60.0)
    dMD = MbarD + (-MbarD + MD)*exp(-dt/500.0)


    taualpha=2000
    gaussian_spread=0.001

    # reward is Sf. See eqn. (7)
    reward = exp(-(avF**2)/gaussian_spread) * exp(-(avS**2)/gaussian_spread) * exp(-(avD**2)/gaussian_spread) 
    
    #this is \rho    
    offset = juicethres

    #alpha evolution
    dalpha = (alpha_inf(reward, offset) - alpha) * (1./taualpha)

    # this is eq. (6). Average errors
    davF = (1./maxTauG) * ((Fbar-F)-avF)
    davS = (1./maxTauG) * ((Sbar-S)-avS)
    davD = (1./maxTauG) * ((Dbar-D)-avD)

#   Sensors Average Time Evolution  --- Euler  
    avF = avF + (davF*dt)
    avS = avS + (davS*dt)
    avD = avD + (davD*dt)
    
    alpha = alpha + (dalpha * dt)
##   Sensor Time Evolution -- explicit ODE
# =============================================================================
#     dMF = (MbarF - MF)/0.5
#     dHF = (HbarF - HF)/1.50
#     dMS = (MbarS - MS)/50.0
#     dHS = (HbarS - HS)/60.0
#     dMD = (MbarD - MD)/500.0
#     
# =============================================================================
    #print F,S,D , 'mD', MD, 'MbarD', MbarD, 'ICa ', ICa , 'dMD ', dMD

#   Intracellular Calcium Time Evolution -- Exponential Euler
    Ca_inf = (-caF*(iCaT + iCaS) + Ca0)
    dIntCa   = Ca_inf + (IntCa-Ca_inf) * exp(-dt / tauIntCa)
    
    #Pack all evolved state variables in a vector and return 
    activations = [dNaM, dNaH, dCaTM, dCaTH, dCaSM, dCaSH, dHM, dKdM,dKCaM, dAM, dAH]
    sensors_activations=[dMF, dHF, dMS, dHS, dMD] 
    av_sensors=[avF, avS, avD , alpha] 

    y_dot = [dV]
    y_dot.extend(activations)
    y_dot.extend([dIntCa])
    y_dot.extend(conductances)
    y_dot.extend(sensors_activations)
    y_dot.extend(av_sensors)
    return y_dot

def get_g(g,m,h,q,Volt,Erev):
    gd = g*pow(m, q)*h;
    return gd

def alpha_inf(reward, offset):
    delta=100
    # offset=0.075 # THIS IS DYNAMIC
    # offset=-1 # THIS IS OFF
    # offset=10 # THIS IS ALWAYS ON
    ainf = 1/(1 + np.exp(-delta*(-reward +offset)))
    return ainf
