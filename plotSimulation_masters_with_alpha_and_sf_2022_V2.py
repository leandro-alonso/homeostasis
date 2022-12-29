from pylab import *

import os
from auxfunctions import *
import time as ttime
import numpy as np
from scipy import *
from scipy.integrate import odeint
from pylab import *
import math
import sys

import string
import random
import glob
from analizeSolution import *

close('all')

#path to folder with simulations data
pathtosimulations='./simulations/'
spread=0.001

print('plotting')

#create a folder for the plots
pathtoplots='./simulation-plots/'
if not os.path.exists(pathtoplots):os.makedirs(pathtoplots)

modelsfolders= glob.glob(pathtosimulations+'init*/')


# plot all models in folder ./simulations/
for folder in modelsfolders:
	fnames= glob.glob(folder+'*.data.npy')
	fnames.sort(key=os.path.getmtime)

	name=folder.split('/')[-2]
	print("name is", name, folder)
	hashname=name.split('-')[2]
	boundval=folder.split('/')[-2].split('-')[-1]
	print(name)
	delta=60
	sampledt=0.5	
	auxname=name.split('_')[0]
	
	for index in arange(0, len(fnames), delta):
		Krev=-80
		trace=[]
		gsevol=[]
		states=[]
		targets=[]
		Ftg=[]
		Stg=[]
		Dtg=[]

		scores_t=[]
		


		avFt=[]
		avSt=[]
		avDt=[]
		krevs=[]
		
		targets=[0.1, 0.07, 0.07]
		numfiles= len(fnames)
		isbursting=[]
		alphat=[]

		intCa=[]
		Vdist=[]
		for fname in fnames[index:index+delta]:	
			print('loading: ', fname)
			trace_i, argsevol_i, states_i, Ftg_i, Dtg_i,Stg_i, avFt_i,avSt_i,avDt_i,alphat_i ,krevs_i,=load(fname, allow_pickle=True, encoding = 'latin1')
			trace_i=array(trace_i)
			tausecs=10
			tau=tausecs*1000*sampledt
			stamps=arange(tau,len(trace_i)-tau-1,tau/100)
			
			ymin=-80
			ymax=40

			for i in stamps:
				localdata1=trace_i[int(i-tau):int(i+tau)]
				h,xed = histogram(localdata1, range=(ymin,ymax),bins=250)
				Vdist.append(h)
			
			trace.extend(trace_i)
			Ftg.extend(Ftg_i)
			Stg.extend(Stg_i)
			Dtg.extend(Dtg_i)
			avFt.extend(avFt_i)
			avSt.extend(avSt_i)
			avDt.extend(avDt_i)
			alphat.extend(alphat_i)			
			gsevol.extend(argsevol_i)
		avFt=array(avFt)
		avSt=array(avSt)
		avDt=array(avDt)
		reward = exp(-(avFt**2)/spread) * exp(-(avSt**2)/spread) * exp(-(avDt**2)/spread) 

		totalmins = index +len(fnames[index:index+delta])
		fout = pathtoplots+name+'.boundval.'+boundval+'.minstart.'+str(index)+'.minend.'+str(index+delta)

		sampledt=0.5
		spikesthres=-15
		fig,spkts, r= analizeActivity(trace, sampledt,spikesthres)

		Vdist=array(Vdist)
		import matplotlib.pyplot as plt	
		plt.rcParams["font.family"] = "Helvetica"
		



#Actual plot
########################################################################################################
		
		close('all')
		figure(figsize=(6,8))
		# suptitle('model ID: ' +  name)
		plottau=int(600/sampledt)

		numminutes=len(fnames)
		timee = linspace(0, numminutes , len(trace))

		# rr =[plottau,1*(len(trace)/60.),2*(len(trace)/60.), 5*(len(trace)/60.), 10*(len(trace)/60.),30*(len(trace)/60.) , len(trace)-plottau]
		# rr =[plottau,1*(len(trace)/20.),2*(len(trace)/20.), 3*(len(trace)/20.), 5*(len(trace)/20.),10*(len(trace)/20.)]
		rr =[plottau,1*(len(trace)/numminutes),3*(len(trace)/numminutes), 4*(len(trace)/numminutes), 5*(len(trace)/numminutes),9*(len(trace)/numminutes)]
		# rr =[3*(len(trace)/numminutes),5*(len(trace)/numminutes), 7*(len(trace)/numminutes), 10*(len(trace)/numminutes),19*(len(trace)/numminutes),21*(len(trace)/numminutes)]
		# [3,5,7,10,19,21]
		# howmanytraces=len(rr)
		howmanytraces=3
		tracesstamps=[]
		for nn,i in enumerate(rr[0:3]):
		# for nn,i in enumerate(logspace(log(plottau), log(len(trace))-plottau, 5)):		
			subplot2grid((6,howmanytraces),(0,nn),rowspan=1,colspan=1)
			title(str(int(round(timee[int(i)]))),color= 'red')
			# title( str(round(i/(1000/sampledt))) )#+ ' T:' + str(plottau))
			tracesstamps.append(timee[int(i)])
			plot(trace[int(i-plottau):int(i+plottau)], color='black', lw=0.25)
			plot(ones(len(trace[int(i-plottau):int(i+plottau)]))*(-50), ls='dashed', color='red', lw=0.5)            
			ylabel('mV')
			xlabel('time [msec]')
			xlim(0,2*plottau)
			ylim(ymin,ymax)
			# if(nn==9):plot(ones(int(250/sampledt))*(-70), lw=2,color='red')
			print(nn)
			axis('off')

		for nn,i in enumerate(rr[3:6]):
		# for nn,i in enumerate(logspace(log(plottau), log(len(trace))-plottau, 5)):		
			subplot2grid((6,howmanytraces),(1,nn),rowspan=1,colspan=1)
			title(str(int(round(timee[int(i)]))),color= 'red')
			# title( str(round(i/(1000/sampledt))) )#+ ' T:' + str(plottau))
			tracesstamps.append(timee[int(i)])
			plot(trace[int(i-plottau):int(i+plottau)], color='black', lw=0.25)
			plot(ones(len(trace[int(i-plottau):int(i+plottau)]))*(-50), ls='dashed', color='red', lw=0.5)            
			ylabel('mV')
			xlabel('time [msec]')
			xlim(0,2*plottau)
			ylim(ymin,ymax)
			if(nn==0):plot(ones(int(250/sampledt))*(-70), lw=2,color='black')
			print(nn)
			axis('off')


		ax = subplot2grid((6,howmanytraces),(2,0), colspan=howmanytraces)
		# elcolormap='gray_r'
		# totalstamps=shape(Vdist)[0]
		# imshow(flipud(log(Vdist+1).T), aspect='auto',extent=[0,timee[-1]/1.,ymin,ymax], interpolation='spline36',cmap=elcolormap)     
		# plot(ones(totalstamps)*(-50), ls='dashed', color='red',lw=1)
		# clim(0.5, 3)
		# xlim(-0.1,timee[-1]/2.)

		# elcolormap='gray_r'
		totalstamps=shape(Vdist)[0]
		# imshow(flipud(log(Vdist+1).T), aspect='auto',extent=[0,timee[-1]/1.,ymin,ymax], interpolation='spline36',cmap=elcolormap)     
		# plot(ones(totalstamps)*(-50), ls='dashed', color='red',lw=1)
		# clim(0.5, 3)
		# xlim(-0.1,timee[-1]/1.)
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ylabel('Sf')	
		ts=trace
		timee = linspace(index, totalmins, len(ts))
		plot(timee, ts,color='black',lw=2, zorder=10)
		plot(ones(totalstamps)*(-50), ls='dashed', color='red',lw=1)
		# legend(loc="upper left")
		# ylim(-0.1,1.1)
		ylim(ymin,ymax)
		xlim(-0.1,timee[-1]/1.)

		tracesstamps=[1,3,4,5,6,9]
		for tstamp in tracesstamps:
			vlines(tstamp, ymin,ymax,color='red',ls=':',lw=1)
			pad=0
			text(int(tstamp-pad), ymax+2, str(int(tstamp)),color='red')
		plt.gca().xaxis.set_major_locator(plt.NullLocator())
		ylabel('V [mV]')
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ylim(ymin,ymax)
		# xlim(-0.1,timee[-1]+0.1)
		xlim(-0.1,timee[-1]/1.)


		ax = subplot2grid((6,howmanytraces),(3,0), colspan=howmanytraces)
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		glabels = ['gNa', 'gCaT', 'gCaS', 'gH', 'gKd', 'gKCa', 'gA']
		totalsecs=len(trace)*sampledt
		# localtime=linspace(0, totalsecs, len(trace))

		for i,ts in enumerate(array(gsevol).T):
			timee = linspace(0, len(fnames), len(ts))
			plot(timee, ts, label=glabels[i])
		ylabel('uS')
		legend(loc='right')
		yscale('log')
		xlabel('time [min]')
		xlim(-0.1,timee[-1]/1.)

		# ax = subplot2grid((7,howmanytraces),(4,0), colspan=howmanytraces)
		# ts=intCa
		# ax.spines['right'].set_visible(False)
		# ax.spines['top'].set_visible(False)
		# ylabel('intCA')	
		# timee = linspace(index, totalmins, len(ts))

		# meanintCa=convolve(intCa, ones(10000)/10000, mode='same')
		# plot(timee, ts,color='blue',lw=2, zorder=10)
		# plot(timee, meanintCa,color='red',lw=2, zorder=10, ls='dashed')
		# # legend(loc="upper left")
		# # ylim(-0.1,1.1)
		# xlim(-0.1,timee[-1]/1.)


		ax = subplot2grid((6,howmanytraces),(4,0), colspan=howmanytraces)
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ylabel('Sf')	
		ts=reward
		timee = linspace(index, totalmins, len(ts))
		plot(timee, ts,color='blue',lw=2, zorder=10)
		# legend(loc="upper left")
		ylim(-0.1,1.1)
		xlim(-0.1,timee[-1]/1.)

		ax = subplot2grid((6,howmanytraces),(5,0), colspan=howmanytraces)
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		# ylabel('$\alpha$')
		ylabel('alpha')	
		ts=alphat
		timee = linspace(index, totalmins, len(ts))
		plot(timee, ts,color='blue',lw=2, zorder=10)
		# legend(loc="upper left")
		ylim(-0.1,1.1)
		xlim(-0.1,timee[-1]/1.)
		# ylabel('a.u.')
		xlabel('time [min]')

		subplots_adjust(left=0.12, right=0.98, top=0.95,bottom=0.1, hspace=0.3, wspace=0.1)
		savefig(fout+'.gevols-Sf-alpha-MASTER.png',dpi=400)
		savefig(fout+'.gevols-Sf-alpha-MASTER.eps',dpi=400)




