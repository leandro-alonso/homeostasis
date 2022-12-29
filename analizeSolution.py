from pylab import *
import numpy

def getSpiketimes(trace, spikethres):
	# print('finding spiketimes fancy (convolution is expensive)') 
	#FIND SPIKE TIMES FANCY 
	# sstrace= convolve(trace , ones(avtau)/avtau , mode='same') 
	sstrace=array(trace).copy()
	all_loc_max = numpy.r_[True, sstrace[1:] >= sstrace[:-1]] & numpy.r_[sstrace[:-1] >= sstrace[1:], True]
	candidates=argwhere(all_loc_max)
	# spikethres=-30
	#of all local maxima keep those that are bigger than their two neighbors (probably useless) and are above spikethreshold
	spkts=[]
	if(len(candidates)>2):
		for lm in candidates[2:-2]: 
		    #if(len(sstrace)<lm+2 and len(sstrace)>lm-2): 
		    lm=lm[0]
		    if(sstrace[lm-2]<=sstrace[lm] and sstrace[lm+2]<sstrace[lm] and sstrace[lm] > spikethres):
		        spkts.append(lm)	
	return spkts

def getInstantStats(st,isithres, dt):
    dc=[]
    ibi=[]
    sb=[]
    eb=[]
    val=isithres/dt
    for i in arange(0,len(st)-1):        
        if(i==0 or (st[i]-st[i-1] > val)):
            #print 'burst starts here ', st[i]
            starts= st[i]
            sb.append(st[i])
        if(st[i+1]-st[i] > val):
            #print 'burst ends here ', st[i]
            eb.append(st[i])
            ends = st[i]
            period = float(ends-starts+ st[i+1]-st[i])
            dc.append((ends-starts)/period)
            ibi.append(period)
            #temps.append(spktemps[i])
    return array(dc), array(ibi),array(sb),array(eb)

def getStats(trace, sampledt):
	spikethres=-30
	spkts=getSpiketimes(trace, spikethres)
	r = getInstantStats(spkts,50, sampledt) 	
	return r, spkts
# figure()
# plot(trace)
# scatter(spkts,array(trace)[spkts])

def analizeActivity(trace, sampledt,spikethres ):
	# spikethres=-30
	spkts=getSpiketimes(trace, spikethres)
	r = getInstantStats(spkts,50, sampledt) 

	fig = figure(figsize=(15,8))
	subplot(221)
	title('< duty cycle > ' + str(mean(r[0])))
	hist(r[0])

	subplot(222)
	title('< burst freq > ' + str(mean(1000./r[1])))
	hist(1000./r[1])

	subplot(223)
	scatter(r[3],r[0])
	ylabel('duty cycle')
	subplot(224)
	ylabel('period [msec]')
	scatter(r[3],r[1])
	subplots_adjust(left=0.1, right=0.95, top=0.90,bottom=0.05, hspace=0.35, wspace=0.25)
	return fig,spkts, r



