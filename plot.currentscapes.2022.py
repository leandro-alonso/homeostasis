from pylab import *
from saveCurrents_bound_dynTauG_V2 import *
import os
from currents_visualization import * 
import glob

pathtostore='/Users/leandro/4.0.local.workspace/results/homeostasis_2022/currentscapes_multiple_injections/'

targets=[0.10, 0.07, 0.07]
boundval=pow(10,-5)*1
maxtauG=2000

#PLOT CURRENTSCAPE CONTROL CONDITION -- MULTIPLE INJECTIONS
pathto_initialcondition = '/Users/leandro/4.0.local.workspace/results/homeostasis_2022/MULTIPLE_DEPOL_INJCURR-taug-2sec-singlebound-e-5-0.1-0.07-0.07-INITCONDITIONS-intCA-REDO-10mins/init-fullstate-LL82XM-liu-cubic1e-05.schemeLIU-jthres-0.075-Inj_curr-0.4-tauG-2000/run.LL82XM-liu-cubic1e-05.schemeLIU.nsecs.60.laststate.npy'
minute=0
injCurrent=0
# getCurrentscapeData_EK(pathto_initialcondition, highK, boundval, maxtauG, targets)
voltage, currents = getCurrentscapeData_injCurr(pathto_initialcondition, False, boundval, maxtauG, targets, injCurrent)
ini=0
fig = plotCurrentscape(array(voltage)[ini:ini+10000], array(currents)[ini:ini+10000,:])
fig.savefig(pathtostore+'currentscape.minute.'+str(minute)+'.png',dpi=250)

#MINUTE 10
pathto_initialcondition = '/Users/leandro/4.0.local.workspace/results/homeostasis_2022/MULTIPLE_DEPOL_INJCURR-taug-2sec-singlebound-e-5-0.1-0.07-0.07-INITCONDITIONS-intCA-REDO-10mins/init-fullstate-LL82XM-liu-cubic1e-05.schemeLIU-jthres-0.075-Inj_curr-0.4-tauG-2000/run.LL82XM-liu-cubic1e-05.schemeLIU.nsecs.600.laststate.npy'
minute=10
injCurrent=0.4
# getCurrentscapeData_EK(pathto_initialcondition, highK, boundval, maxtauG, targets)
voltage, currents = getCurrentscapeData_injCurr(pathto_initialcondition, False, boundval, maxtauG, targets, injCurrent)
ini=0
fig = plotCurrentscape(array(voltage)[ini:ini+10000], array(currents)[ini:ini+10000,:])
fig.savefig(pathtostore+'currentscape.minute.'+str(minute)+'.png',dpi=250)

# fig.title('10 mins')


#MINUTE 20
pathto_initialcondition = '/Users/leandro/4.0.local.workspace/results/homeostasis_2022/MULTIPLE_DEPOL_INJCURR-taug-2sec-singlebound-e-5-0.1-0.07-0.07-INITCONDITIONS-intCA-REDO-10mins/init-fullstate-LL82XM-liu-cubic1e-05.schemeLIU-jthres-0.075-Inj_curr-0.4-tauG-2000/run.LL82XM-liu-cubic1e-05.schemeLIU.nsecs.1200.laststate.npy'
injCurrent=0
# getCurrentscapeData_EK(pathto_initialcondition, highK, boundval, maxtauG, targets)
voltage, currents = getCurrentscapeData_injCurr(pathto_initialcondition, False, boundval, maxtauG, targets, injCurrent)
ini=0
fig = plotCurrentscape(array(voltage)[ini:ini+10000], array(currents)[ini:ini+10000,:])
minute=20
fig.savefig(pathtostore+'currentscape.minute.'+str(minute)+'.png',dpi=250)



#MINUTE 30
pathto_initialcondition ='/Users/leandro/4.0.local.workspace/results/homeostasis_2022/MULTIPLE_DEPOL_INJCURR-taug-2sec-singlebound-e-5-0.1-0.07-0.07-INITCONDITIONS-intCA-REDO-10mins/init-fullstate-LL82XM-liu-cubic1e-05.schemeLIU-jthres-0.075-Inj_curr-0.4-tauG-2000/run.LL82XM-liu-cubic1e-05.schemeLIU.nsecs.1800.laststate.npy'
injCurrent=0.4
# getCurrentscapeData_EK(pathto_initialcondition, highK, boundval, maxtauG, targets)
voltage, currents = getCurrentscapeData_injCurr(pathto_initialcondition, False, boundval, maxtauG, targets, injCurrent)
ini=0
fig = plotCurrentscape(array(voltage)[ini:ini+10000], array(currents)[ini:ini+10000,:])
minute=30
fig.savefig(pathtostore+'currentscape.minute.'+str(minute)+'.png',dpi=250)


#MINUTE 40
pathto_initialcondition ='/Users/leandro/4.0.local.workspace/results/homeostasis_2022/MULTIPLE_DEPOL_INJCURR-taug-2sec-singlebound-e-5-0.1-0.07-0.07-INITCONDITIONS-intCA-REDO-10mins/init-fullstate-LL82XM-liu-cubic1e-05.schemeLIU-jthres-0.075-Inj_curr-0.4-tauG-2000/run.LL82XM-liu-cubic1e-05.schemeLIU.nsecs.2400.laststate.npy'
injCurrent=0.
# getCurrentscapeData_EK(pathto_initialcondition, highK, boundval, maxtauG, targets)
voltage, currents = getCurrentscapeData_injCurr(pathto_initialcondition, False, boundval, maxtauG, targets, injCurrent)
ini=0
fig = plotCurrentscape(array(voltage)[ini:ini+10000], array(currents)[ini:ini+10000,:])
minute=40
fig.savefig(pathtostore+'currentscape.minute.'+str(minute)+'.png',dpi=250)


#MINUTE 50
pathto_initialcondition ='/Users/leandro/4.0.local.workspace/results/homeostasis_2022/MULTIPLE_DEPOL_INJCURR-taug-2sec-singlebound-e-5-0.1-0.07-0.07-INITCONDITIONS-intCA-REDO-10mins/init-fullstate-LL82XM-liu-cubic1e-05.schemeLIU-jthres-0.075-Inj_curr-0.4-tauG-2000/run.LL82XM-liu-cubic1e-05.schemeLIU.nsecs.3000.laststate.npy'
injCurrent=0.4
# getCurrentscapeData_EK(pathto_initialcondition, highK, boundval, maxtauG, targets)
voltage, currents = getCurrentscapeData_injCurr(pathto_initialcondition, False, boundval, maxtauG, targets, injCurrent)
ini=0
fig = plotCurrentscape(array(voltage)[ini:ini+10000], array(currents)[ini:ini+10000,:])
minute=50
fig.savefig(pathtostore+'currentscape.minute.'+str(minute)+'.png',dpi=250)


#MINUTE 60
pathto_initialcondition ='/Users/leandro/4.0.local.workspace/results/homeostasis_2022/MULTIPLE_DEPOL_INJCURR-taug-2sec-singlebound-e-5-0.1-0.07-0.07-INITCONDITIONS-intCA-REDO-10mins/init-fullstate-LL82XM-liu-cubic1e-05.schemeLIU-jthres-0.075-Inj_curr-0.4-tauG-2000/run.LL82XM-liu-cubic1e-05.schemeLIU.nsecs.3600.laststate.npy'
injCurrent=0
# getCurrentscapeData_EK(pathto_initialcondition, highK, boundval, maxtauG, targets)
voltage, currents = getCurrentscapeData_injCurr(pathto_initialcondition, False, boundval, maxtauG, targets, injCurrent)
ini=0
fig = plotCurrentscape(array(voltage)[ini:ini+10000], array(currents)[ini:ini+10000,:])
minute=60
fig.savefig(pathtostore+'currentscape.minute.'+str(minute)+'.png',dpi=250)



# lag=10000
# V = array(voltage[ini:ini+lag])
# C = array(currents)[:, ini:ini+lag]
# fig = plotCurrentscape(V,C)
