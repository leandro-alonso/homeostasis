from pylab import *

#THIS IS THE MOST UPDATED VERSION, Y AXIS IS INVERTED, MARGINS ARE FIXED. YLIMITS ARE CHOSEN FOR LIU ADAPTIVE
def plotCurrentscape(voltage, currents):	
	# make a copy of currents
	# CURRENTSCAPE CALCULATION STARTS HERE. 
	curr=array(currents)	
	cpos= curr.copy()
	cpos[curr<0]=0
	cneg= curr.copy()
	cneg[curr>0]=0

	normapos = sum(abs(array(cpos)),axis=0)
	normaneg = sum(abs(array(cneg)),axis=0)
	npPD=normapos
	nnPD=normaneg
	cnorm=curr.copy()
	cnorm[curr>0]=(abs(curr)/normapos)[curr>0]
	cnorm[curr<0]=-(abs(curr)/normaneg)[curr<0]

	resy=1000
	impos=zeros((resy,shape(cnorm)[-1])) 
	imneg=zeros((resy,shape(cnorm)[-1])) 

	times=arange(0,shape(cnorm)[-1])
	for t in times:
	    lastpercent=0
	    for numcurr, curr in enumerate(cnorm):
	        if(curr[t]>0):
	            percent = int(curr[t]*(resy))   
	            impos[lastpercent:lastpercent+percent,t]=numcurr
	            lastpercent=lastpercent+percent        
	for t in times:
	    lastpercent=0
	    for numcurr, curr in enumerate(cnorm):
	        if(curr[t]<0):
	            percent = int(abs(curr[t])*(resy))   
	            imneg[lastpercent:lastpercent+percent,t]=numcurr
	            lastpercent=lastpercent+percent        
	im0= vstack((impos,imneg))
	# CURRENTSCAPE CALCULATION ENDS HERE. 

	#PLOT CURRENTSCAPE
	# fig = figure(figsize=(3,4))
	fig = figure(figsize=(6,4))

	#PLOT VOLTAGE TRACE
	xmax=len(voltage)
	swthres=-50        
	ax=subplot2grid((7,1),(0,0),rowspan=2)	
	t=arange(0,len(voltage))
	plot(t, voltage, color='black',lw=1.)
	plot(t,ones(len(t))*swthres,ls='dashed',color='black',lw=0.75)
	vlines(1,-50,-20,lw=1)
	ylim(-80,30)
	xlim(0,xmax)
	axis('off')         

	#PLOT TOTAL INWARD CURRENT IN LOG SCALE
	ax=subplot2grid((7,1),(2,0),rowspan=1)
	fill_between(arange(len((npPD))),(npPD),color='black')
	plot(1.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(10.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(100.*ones(len(nnPD)),color='black', ls=':',lw=1)
	yscale('log')
	ylim(0.001,1000)
	xlim(0,xmax)
	axis('off') 

	#PLOT CURRENT SHARES
	elcolormap='Set1'
	ax=subplot2grid((7,1),(3,0),rowspan=3)
	imshow(im0[::1,::1],interpolation='nearest',aspect='auto',cmap=elcolormap)
	ylim(2*resy,0)
	plot(resy*ones(len(npPD)),color='black',lw=2)
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	plt.gca().yaxis.set_major_locator(plt.NullLocator())
	xlim(0,xmax)
	clim(0,8)
	axis('off') 

	#PLOT TOTAL OUTWARD CURRENT IN LOG SCALE
	ax=subplot2grid((7,1),(6,0),rowspan=1)
	fill_between(arange(len((nnPD))),(nnPD),color='black')
	plot(1.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(10.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(100.*ones(len(nnPD)),color='black', ls=':',lw=1)
	yscale('log')
	ylim(1000,0.001)
	# plt.gca().invert_yaxis()
	xlim(0,xmax)
	axis('off') 
	subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
	return fig

#THIS IS THE MOST UPDATED VERSION, Y AXIS IS INVERTED, MARGINS ARE FIXED. YLIMITS ARE CHOSEN FOR LIU ADAPTIVE
def plotCurrentscape_AltSize(voltage, currents, xsize):	
	# make a copy of currents
	# CURRENTSCAPE CALCULATION STARTS HERE. 
	curr=array(currents)	
	cpos= curr.copy()
	cpos[curr<0]=0
	cneg= curr.copy()
	cneg[curr>0]=0

	normapos = sum(abs(array(cpos)),axis=0)
	normaneg = sum(abs(array(cneg)),axis=0)
	npPD=normapos
	nnPD=normaneg
	cnorm=curr.copy()
	cnorm[curr>0]=(abs(curr)/normapos)[curr>0]
	cnorm[curr<0]=-(abs(curr)/normaneg)[curr<0]

	resy=1000
	impos=zeros((resy,shape(cnorm)[-1])) 
	imneg=zeros((resy,shape(cnorm)[-1])) 

	times=arange(0,shape(cnorm)[-1])
	for t in times:
	    lastpercent=0
	    for numcurr, curr in enumerate(cnorm):
	        if(curr[t]>0):
	            percent = int(curr[t]*(resy))   
	            impos[lastpercent:lastpercent+percent,t]=numcurr
	            lastpercent=lastpercent+percent        
	for t in times:
	    lastpercent=0
	    for numcurr, curr in enumerate(cnorm):
	        if(curr[t]<0):
	            percent = int(abs(curr[t])*(resy))   
	            imneg[lastpercent:lastpercent+percent,t]=numcurr
	            lastpercent=lastpercent+percent        
	im0= vstack((impos,imneg))
	# CURRENTSCAPE CALCULATION ENDS HERE. 

	#PLOT CURRENTSCAPE
	fig = figure(figsize=(xsize,4))

	#PLOT VOLTAGE TRACE
	xmax=len(voltage)
	swthres=-50        
	ax=subplot2grid((7,1),(0,0),rowspan=2)	
	t=arange(0,len(voltage))
	plot(t, voltage, color='black',lw=1.)
	plot(t,ones(len(t))*swthres,ls='dashed',color='black',lw=0.75)
	vlines(1,-50,-20,lw=1)
	ylim(-70,30)
	xlim(0,xmax)
	axis('off')         

	#PLOT TOTAL INWARD CURRENT IN LOG SCALE
	ax=subplot2grid((7,1),(2,0),rowspan=1)
	fill_between(arange(len((npPD))),(npPD),color='black')
	plot(1.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(10.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(100.*ones(len(nnPD)),color='black', ls=':',lw=1)
	yscale('log')
	ylim(0.001,1000)
	xlim(0,xmax)
	axis('off') 

	#PLOT CURRENT SHARES
	elcolormap='Set1'
	ax=subplot2grid((7,1),(3,0),rowspan=3)
	imshow(im0[::1,::1],interpolation='nearest',aspect='auto',cmap=elcolormap)
	ylim(2*resy,0)
	plot(resy*ones(len(npPD)),color='black',lw=2)
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	plt.gca().yaxis.set_major_locator(plt.NullLocator())
	xlim(0,xmax)
	clim(0,8)
	axis('off') 

	#PLOT TOTAL OUTWARD CURRENT IN LOG SCALE
	ax=subplot2grid((7,1),(6,0),rowspan=1)
	fill_between(arange(len((nnPD))),(nnPD),color='black')
	plot(1.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(10.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(100.*ones(len(nnPD)),color='black', ls=':',lw=1)
	yscale('log')
	ylim(1000,0.001)
	plt.gca().invert_yaxis()
	xlim(0,xmax)
	axis('off') 
	subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
	return fig



#THIS IS THE MOST UPDATED VERSION, Y AXIS IS INVERTED, MARGINS ARE FIXED. YLIMITS ARE CHOSEN FOR LIU ADAPTIVE
def plotCurrentscape_AltSize_moretrace(voltage, currents, xsize):	
	# make a copy of currents
	# CURRENTSCAPE CALCULATION STARTS HERE. 
	curr=array(currents)	
	cpos= curr.copy()
	cpos[curr<0]=0
	cneg= curr.copy()
	cneg[curr>0]=0

	normapos = sum(abs(array(cpos)),axis=0)
	normaneg = sum(abs(array(cneg)),axis=0)
	npPD=normapos
	nnPD=normaneg
	cnorm=curr.copy()
	cnorm[curr>0]=(abs(curr)/normapos)[curr>0]
	cnorm[curr<0]=-(abs(curr)/normaneg)[curr<0]

	resy=1000
	impos=zeros((resy,shape(cnorm)[-1])) 
	imneg=zeros((resy,shape(cnorm)[-1])) 

	times=arange(0,shape(cnorm)[-1])
	for t in times:
	    lastpercent=0
	    for numcurr, curr in enumerate(cnorm):
	        if(curr[t]>0):
	            percent = int(curr[t]*(resy))   
	            impos[lastpercent:lastpercent+percent,t]=numcurr
	            lastpercent=lastpercent+percent        
	for t in times:
	    lastpercent=0
	    for numcurr, curr in enumerate(cnorm):
	        if(curr[t]<0):
	            percent = int(abs(curr[t])*(resy))   
	            imneg[lastpercent:lastpercent+percent,t]=numcurr
	            lastpercent=lastpercent+percent        
	im0= vstack((impos,imneg))
	# CURRENTSCAPE CALCULATION ENDS HERE. 

	#PLOT CURRENTSCAPE
	fig = figure(figsize=(xsize,5))

	#PLOT VOLTAGE TRACE
	xmax=len(voltage)
	swthres=-50        
	ax=subplot2grid((10,1),(0,0),rowspan=5)	
	t=arange(0,len(voltage))
	plot(t, voltage, color='black',lw=1.)
	plot(t,ones(len(t))*swthres,ls='dashed',color='black',lw=0.75)
	vlines(1,-50,-20,lw=1)
	ylim(-60,20)
	xlim(0,xmax)
	axis('off')         

	#PLOT TOTAL INWARD CURRENT IN LOG SCALE
	ax=subplot2grid((10,1),(5,0),rowspan=1)
	fill_between(arange(len((npPD))),(npPD),color='black')
	plot(1.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(10.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(100.*ones(len(nnPD)),color='black', ls=':',lw=1)
	yscale('log')
	ylim(0.001,1000)
	xlim(0,xmax)
	axis('off') 

	#PLOT CURRENT SHARES
	elcolormap='Set1'
	ax=subplot2grid((10,1),(6,0),rowspan=3)
	imshow(im0[::1,::1],interpolation='nearest',aspect='auto',cmap=elcolormap)
	ylim(2*resy,0)
	plot(resy*ones(len(npPD)),color='black',lw=2)
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	plt.gca().yaxis.set_major_locator(plt.NullLocator())
	xlim(0,xmax)
	clim(0,8)
	axis('off') 

	#PLOT TOTAL OUTWARD CURRENT IN LOG SCALE
	ax=subplot2grid((10,1),(9,0),rowspan=1)
	fill_between(arange(len((nnPD))),(nnPD),color='black')
	plot(1.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(10.*ones(len(nnPD)),color='black', ls=':',lw=1)
	plot(100.*ones(len(nnPD)),color='black', ls=':',lw=1)
	yscale('log')
	ylim(1000,0.001)
	plt.gca().invert_yaxis()
	xlim(0,xmax)
	axis('off') 
	subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
	return fig



def plotVoltageDistributions(Vdist): 
    im=Vdist
    fig = figure()
    cmmap='Greys'		
    #proper visualization requires that we take the logarithm (we add +1 to avoid NaN)       
    imshow((log10(im+1)),aspect='auto', cmap=cmmap, extent=[1,0,-75,30], interpolation='nearest')
    xlim(1,0)
    ylim(-75,30)
    clim(1,5)
    axis('off')     
    subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)    
    return fig

def plotVoltageDistributionsEnhanceEdges(Vdist): 
    im=Vdist
    fig = figure(figsize=(3,6))    
    #choose a colormap (see also gnuplot1, gnuplot2, helix)
    cmmap='hot'
    #proper visualization requires that we take the logarithm     
    a=log10(im+1)
    #We smooth the distribution using a convolution with the identity. This does not have a big effect. 
    filt=ones(3)/3.
    r=np.apply_along_axis(lambda m: np.convolve(m, filt, mode='same'), axis=0, arr=a)
    #We enhance edges by performing a derivative (along the V-axis). 
    imm=abs(diff(r,axis=0))
    imshow(imm,aspect='auto', cmap=cmmap, extent=[1,0,-75,30], interpolation='nearest')     
    xlim(1,0)
    ylim(-75,30)
    axis('off')     
    subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)    
    clim(0.0,0.3)
    return fig

def plotCurrentSharesDistributions(current_share_dist): 
    im=current_share_dist
    fig = figure()
    cmmap='gnuplot2'		
    #proper visualization requires that we take the logarithm (we add +1 to avoid NaN)       
    imshow(log10(im+1),aspect='auto', cmap=cmmap, extent=[1,0,0,1], interpolation='nearest')
    percents=linspace(1,0,101)
    plot(percents,ones(len(percents))*(0.25), color = 'white',ls=':')    
    plot(percents,ones(len(percents))*(0.5), color = 'white',ls=':')    
    plot(percents,ones(len(percents))*(0.75), color = 'white',ls=':')   
    xlim(1,0)
    ylim(0,1)
    clim(1,5)
    axis('off')     
    subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)    
    return fig    
   



	