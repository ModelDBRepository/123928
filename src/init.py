# This file contains all the parameters.
# Additional comments are provided in the header of main.py

print '********'
print '* Init *'
print '********'

import os
isunix = lambda: os.name == 'posix'
import matplotlib   
if isunix():
    matplotlib.use('Agg') # no graphic link with Unix cluster for me

import brian_no_units
from brian import * 
from scipy.io import *
from scipy import weave
from time import time
from time import localtime
from customrefractoriness import *
import pickle
import glob


set_global_preferences(useweave=True) # to compile C code

globalStartTime=time()*second


#*************************
# COMPUTATION PARAMETERS *
#*************************
neuronTimeCompression = 2.0**0 # to compress all the neuronal time constants
pbTimeCompression = 2.0**0 # to compress the problem and oscillation time constants

get_default_clock().set_dt(.1*ms/neuronTimeCompression) # set time step

# random state
randState = 28 # use this to specify it
# # use this for batches:
#rl = glob.glob('../data/rand*.mat')
#if size(rl)==0:
#    randState = 0
#else:
#    randState = int(rl[-1][12:15])+1
#savemat('../data/rand'+'%03d' % (randState)+'.mat',{'tmp':[]})
seed(randState)

imposedEnd = 3*second/pbTimeCompression # imposed end time
N = 2000 # number of presynaptic neurons
nG = 1 # number of gmax values
nR = 1 # number of ratio LTD/LTP
M = nG*nR # number of postsynaptic neurons, numbered like that [ (r_0,g_0)...(r_0,g_nG),(r_1,g_0)...(r_1,g_nG),...,(r_nR,g_0)...(r_nR,g_nG)]

recomputeSpikeList = True # recompute input spike trains as opposed to load appropriate mat files
dumpSpikeList = False # save input spike list in mat files
computeOutput = True # compute output (STDP) layer

useSavedWeight = False # load previously dumped weights
timeOffset = 0*second # simulation starts at t=timeOffset

useReset = False # impose resets on input layer at dates specified in reset.###.mat file

graph = True # graph output
monitorInput = True # monitor input spikes
monitorOutput = True # monitor output spikes
monitorPot = True # monitor potential in output layer
monitorCurrent = False # monitor potential in output layer
monitorInputPot = False # monitor potential in input layer
monitorRate = False # monitor rates in output layer
isMonitoring = False # flag saying if currently monitoring
monitorTime = (imposedEnd-6/pbTimeCompression)*second # start monitoring only at that time (to save memory)
analyzePeriod = Inf # periodically launches analyze.py

if not recomputeSpikeList and dumpSpikeList:
    print 'Warning: dumping a spike list which is not re-computed makes no sense. Setting dumpSpikeList to False'
    dumpSpikeList = False

if not computeOutput and ( monitorOutput or monitorPot or monitorCurrent or monitorRate):
    print 'Warning: can not monitor output, which is not computed. Setting monitorOutput, monitorPot and monitorCurrent to False'
    monitorOutput = False
    monitorPot = False
    monitorCurrent = False
    monitorRate = False

# load pattern values (the values are only useful for plotting), but the length of the vector is used to scale gmax
if os.path.exists(os.path.join('..','data','realValuedPattern.'+'%03d' % (randState)+'.mat')):
    realValuedPattern=loadmat(os.path.join('..','data','realValuedPattern.'+'%03d' % (randState)+'.mat'))
    realValuedPattern=realValuedPattern['realValuedPattern']
else:
    realValuedPattern=zeros(round(.5*N))
    
#**************************
# NEURON MODEL PARAMETERS *
#**************************

# Types of input and output neurons
poissonInput = False # poisson input neurons (as opposed to LIF)
conductanceOutput = False # conductance-base output neurons (as opposed to LIF). Tim: 9/2008: has never been critical so far
poissonOutput = False # stochastic (Poisson) output neurons (as opposed to deterministic). Note that their differential equations are the same as the LIF, but firing is stochastic.

#neurons (Dayan&Abbott 2001)
refractoryPeriod = 1*ms/neuronTimeCompression
R=10*Mohm
if poissonOutput:
    vt=400 # not used (just for graph scaling)
    vr=0 # not used
    El=-450 # resting
#    El=-280 # resting (use that value to have more false alarms)
else:
    vt=-54*mV # threshold
    vr=-60*mV # reset
    El=-70*mV # resting
    Ee=0*mV # for excitatory conductance
taum = 20*ms/neuronTimeCompression # membrane time constant
taue=taum/4 # synapse time constant

sigma = 0.015*(vt-vr) # white Gaussian noise. Applies to input and output neurons

# array of max conductance values
unitaryEffect = taue/(taum-taue)*((taum/taue)**(-taue/(taum-taue))-(taum/taue)**(-taum/(taum-taue))) # this corresponds to the maximum of the kernel (exp(-t/taum) - exp(-t/taue)) taue / (taum-taue)
if poissonOutput:
#    gmax=7*100.0/(1.0*N)*1.0/unitaryEffect*exp(log(2)*(array(nR*range(nG))-nG/2)) # appropriate for (1000/2000 afferents, LIF+reset)
    gmax=1.05**-10*2.5*1500/(1.0*size(realValuedPattern))*1.0/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # appropriate for (1000/2000 afferents, PLoS patten)
#    gmax=1.05**-16*2.5*1500/(1.0*size(realValuedPattern))*1.0/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # appropriate for (1000/2000 afferents, N1N2)
else:
    
#    LIF+reset every 250ms, all spikes
#     gmax = 1.05**12/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=2.5%
#     gmax = 1.05**18/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=5%
#     gmax = 1.05**24/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=10%
#     gmax = 1.05**30/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=20%
#     gmax = 1.05**36/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=40%
     
#    LIF+reset every 125ms, all spikes (Tim 03/09) 
#    gmax = 1.05**2/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=2.5%
#    gmax = 1.05**8/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=5%
#    gmax = 1.05**14/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=10%
#    gmax = 1.05**20/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=20%
#    gmax = 1.05**26/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=40%

#    8Hz sin oscillations with 1-3 spikes/cycle, nearest spike
#    gmax = 1.05**4/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) ## time compression x2^-2
#    gmax = 1.05**4/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # time compression x2^-1.5
#    gmax = 1.05**2/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # time compression x2^-1
#    gmax = 1.05**-2/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # time compression x1.0
#    gmax = 1.05**4/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) ## time compression x2^1.25
#    gmax = 1.05**4/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # time compression x2^1.5
#    gmax = 1.05**-2/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*2*log(1.05)*(array(nR*range(nG))-nG/2)) ## time compression x2^1.75
#    gmax = 1.05**6/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*2*log(1.05)*(array(nR*range(nG))-nG/2)) ## time compression x2^2.0
#    gmax = 1.05**6/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*2*log(1.05)*(array(nR*range(nG))-nG/2)) ## time compression x2^2.5

#    8Hz sin oscillations with 1-3 spikes/cycle, all spikes
#    gmax = 1.05**-16/(size(realValuedPattern))*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=2.5%
#    gmax = 1.05**-8/(size(realValuedPattern))*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=5%
    gmax = 1.05**-2/(size(realValuedPattern))*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=10%
#    gmax = 1.05**8/(size(realValuedPattern))*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=20%
#    gmax = 1.05**12/(size(realValuedPattern))*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # x=40%

#    gmax = 1.05**4/(size(realValuedPattern))*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) ## time compression x 2^-2
#    gmax = 1.05**2/(size(realValuedPattern))*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # time compression x 2^-1
#    gmax = 1.05**-2/(size(realValuedPattern))*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) ## time compression x 1.0
#    gmax = 1.05**-6/(size(realValuedPattern))*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) # time compression x 2
#    gmax = 1.05**-10/(size(realValuedPattern))*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) ## time compression x 2^2
#    gmax = 1.05**-14/(size(realValuedPattern))*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) ## time compression x 2^2.5
#    gmax = 1.05**-22/(size(realValuedPattern))*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) ## time compression x 2^3

    # academic patterns (Matthieu Gilson)
#    gmax = 1.05**30/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) ## x=25% LIF
#    gmax = 1.05**38/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) ## x=10% LIF
#    gmax = 1.05**34/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) ## x=6.25% LIF
#    gmax = 1.05**50/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) ## x=100% of 200 afferents, LIF
#    gmax = 1.05**38/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(2*log(1.05)*(array(nR*range(nG))-nG/2)) ## x=100% of 20 afferents, LIF
#    gmax = 1.05**30/size(realValuedPattern)*(vt-El)/unitaryEffect*exp(3*log(1.05)*(array(nR*range(nG))-nG/2)) ## x=50% LIF

if conductanceOutput:
    gmax /= (Ee-vt)

print 'gmax=' + str(gmax)

#********************
# WEIGHT PARAMETERS *
#********************
# initial synaptic weight are randomly picked (uniformly) between those two bounds
if poissonOutput:
    initialWeight_min = 0
    initialWeight_max = 10    
else:
    initialWeight_min = 0*volt
    initialWeight_max = 2*8.6e-5*volt
burstingCriterion = .5 # unplug neurons whose mean normalized synaptic weight is above this value. This allows to save memory by not computing neurons whose normalized synaptic weights all go to 1.

#***************************
# INPUT CURRENT PARAMETERS *
#***************************    
# used only if recomputeSpikeList = True

# input current (normalization with respect to raw values in input file)
if poissonInput: # in this case I values are in fact firing rates
    Imin = 0
    Imax = 100
else:

#    Imax = (1.05)*(vt-El)/R*ones(N) # Tim 12/08: appropriate for (LIF+reset)
#    Imin = (1.0)*(vt-El)/R*ones(N) # Tim 12/08: appropriate for (LIF+reset)  
  
#    Imax = (1.140)*(vt-El)/R # Tim 12/08: appropriate for (8Hz sawtooth oscillations). up to 3 spikes/cycle
#    Imin = (1.025)*(vt-El)/R  # Tim 12/08: appropriate for (8Hz sawtooth oscillations). up to 3 spikes/cycle

    Imax = (1.07)*(vt-El)/R*ones(N) # Tim 12/08: appropriate for (8Hz sin oscillations). 1-3 spikes/cycle
    Imin = (0.95)*(vt-El)/R*ones(N) # Tim 12/08: appropriate for (8Hz sin oscillations). 1-3 spikes/cycle

#    Imax = (1.012)*(vt-El)/R*ones(N) # Tim 12/08: appropriate for (8Hz sin oscillations). 1-3 spikes/cycle. taum=10ms
#    Imin = (0.936)*(vt-El)/R*ones(N) # Tim 12/08: appropriate for (8Hz sin oscillations). 1-3 spikes/cycle. taum=10ms

#    Imax = (1.2)*(vt-El)/R*ones(N) # 
#    Imin = (0.95)*(vt-El)/R*ones(N) #

#    Imax = (vt-El)/R * (1.07)*ones(N) # partial drive
#    Imin = (vt-El)/R * (0.95)*ones(N) # partial drive

#    Imax = (0.992)*(vt-El)/R # Tim 12/08: appropriate for (8Hz sin oscillations). 1 spike/cycle
#    Imin = (0.95)*(vt-El)/R # Tim 12/08: appropriate for (8Hz sin oscillations). 1 spike/cycle



# oscillation amplitude (use 0 not to have any). An array is used so that different amplitudes may be used for the input neurons
#a = .075*(vt-El)/R * array(12*range(9))/8.0
#for i in range(9):
#    a[i*12:(i+1)*12] = .075*(vt-El)/R * (i/8.0)
a = .075*(vt-El)/R*ones(N)
# oscillation fequence
oscilFreq = 8*Hz*pbTimeCompression

#******************
# STDP PARAMETERS *
#******************    
nearestSpike = False # Tim 08/2008: not critical with low input spike rates (10Hz)
tau_post=33.7*ms #(source: Bi & Poo 2001)
tau_pre=16.8*ms #(source: Bi & Poo 2001)
#tau_post=20*ms #(source: Song, Miller & Abbott 2000)
#tau_pre=20*ms #(source: Song, Miller & Abbott 2000)
a_pre= .005 # Tim 12/08: no more that .005
w_out = - 0*.005/3 # see Matt Gilson
w_in=zeros(M) # array of w_in
for i in range(nR): # here nR corresponds to number of w_in/w_out ratios (and not number of LTD/LTP ratios)
    # see Matt Gilson
    w_in[i*nG:(i+1)*nG] = -w_out*1.05**0*exp((i-nR/2)*.5*log(2))
a_post=zeros(M) # array of LTD/LTP ratios
for i in range(nR):
    
#    LIF+reset every 250ms, all spikes
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**-5*exp((i-nR/2)*log(1.05)) # x=2.5%, 5% 10%
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**-2*exp((i-nR/2)*log(1.05)) # x=20%
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**0*exp((i-nR/2)*log(1.05)) # x=40%

#    LIF+reset every 125ms, all spikes
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**-1*exp((i-nR/2)*1*log(1.05)) #
   
#    8Hz sin oscillations with 1-3 spikes/cycle, nearest spike
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**3*exp((i-nR/2)*1*log(1.05)) # time compression x2^-2
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**5*exp((i-nR/2)*1*log(1.05)) # time compression x2^-1.5
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**8*exp((i-nR/2)*1*log(1.05)) # time compression x2^-1
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**9*exp((i-nR/2)*1*log(1.05)) # time compression x1.0
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**7*exp((i-nR/2)*1*log(1.05)) # time compression x2^1
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**7*exp((i-nR/2)*1*log(1.05)) ## time compression x2^1.25
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**5*exp((i-nR/2)*1*log(1.05)) # time compression x2^1.5
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**5*exp((i-nR/2)*1*log(1.05)) ## time compression x2^1.75
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**5*exp((i-nR/2)*1*log(1.05)) ## time compression x2^2
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**3*exp((i-nR/2)*1*log(1.05)) # time compression x2^2.5

#    8Hz sin oscillations with 1-3 spikes/cycle, all spikes
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**9*exp((i-nR/2)*1*log(1.05)) # x=2.5%
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**9*exp((i-nR/2)*1*log(1.05)) # x=5%
    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**8*exp((i-nR/2)*1*log(1.05)) # x=10%
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**8*exp((i-nR/2)*1*log(1.05)) # x=20%
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**8*exp((i-nR/2)*1*log(1.05)) # x=40%

#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**3*exp((i-nR/2)*.5*log(1.05)) ## time compression x2^-2
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**7*exp((i-nR/2)*.5*log(1.05)) # time compression x2^-1
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**8*exp((i-nR/2)*.5*log(1.05)) ## time compression x1.0    
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**2*exp((i-nR/2)*.5*log(1.05)) # time compression x2
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**-5*exp((i-nR/2)*.5*log(1.05)) ## time compression x2^2
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**-8*exp((i-nR/2)*.5*log(1.05)) ## time compression x2^2.5
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**-10*exp((i-nR/2)*.5*log(1.05)) ## time compression x2^3

    # academic patterns (Matthieu Gilson)
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**-11*exp((i-nR/2)*1*log(1.05)) # Tim 05/09: x=25% LIF 
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**-7*exp((i-nR/2)*1*log(1.05)) # Tim 05/09: x=10% LIF 
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**-10*exp((i-nR/2)*1*log(1.05)) # Tim 06/09: x=100% of 200 or 20 afferents, LIF 
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**-12*exp((i-nR/2)*1*log(1.05)) # Tim 06/09: x=50% Poisson, PLoS pattern
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**-9*exp((i-nR/2)*2*log(1.05)) # Tim 06/09: x=50% Poisson, N1N2 pattern
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**-9*exp((i-nR/2)*1*log(1.05)) # Tim 06/09: x=50% Poisson, PLoS pattern, mu = 0.1
#    a_post[i*nG:(i+1)*nG] = -a_pre*1.05**-11*exp((i-nR/2)*1*log(1.05)) # Tim 07/09: x=50% LIF 
    
    
mu = 0.0 # see Gutig et al. 2003.
print 'normalized a_post/a_pre ratios' + str(a_post/a_pre)



#***********
# FUNCIONS *
#***********

def printtime(mess):
    t = localtime()
    print  '%02d' % t[3] + ':' + '%02d' % t[4] + ' ' + mess

def listMatFile(directory,randState):                                        
    "get list of spike lists" 
    fileList = os.listdir(directory)
    fileList.sort()
    fileList = [f 
               for f in fileList
                # e.g. spikeList.200.010.mat
                #      01234567890123456789
                if f[0:9] in ['spikeList'] and f[10:13] in ['%03d' % (randState)] and f[-4:] in ['.mat'] ]
    return fileList


# Called whenever output neurons fire. Resets the potential and trigger STDP updates.
# Includes C code, will be compiled the first time it is called
# Param:
# P: NeuronGroup (output neurons)
# spikes: list of the indexes of the neuron that fire.
def neurons_reset(P,spikes):
    if size(spikes):
        if not poissonOutput:
            P.v_[spikes]=vr # reset pot
        if nearestSpike:
            nspikes = size(spikes)
            A_pre = mirror.A_pre_
            code = '''
            for(int si=0;si<nspikes;si++)
            {
                int i = spikes(si);
                for(int j=0;j<N;j++)
                {
                    if(!_alreadyPotentiated(j,i))
                    {
                        double wnew;
                        if(mu==0) { /* additive. requires hard bound */
                            wnew = _synW(j,i)+_gmax(i)*(w_out+A_pre(j));
                            if(wnew>_gmax(i)) wnew = _gmax(i);
                            if(wnew<0) wnew = 0.0;
                        }
                        else { /* soft bound */
                            wnew = _synW(j,i)+_gmax(i)*(w_out+A_pre(j)*exp(mu*log(1-_synW(j,i)/_gmax(i))));
                            if(wnew>_gmax(i)) wnew = _gmax(i);
                            if(wnew<0) wnew = 0.0;
						}
                        _synW(j,i) = wnew;
                        _alreadyPotentiated(j,i) = true;
                    }
                }
            }
            '''
            weave.inline(code,
                        ['spikes', 'nspikes', 'N', '_alreadyPotentiated', '_synW', '_gmax', 'A_pre', 'mu','w_out'],
                        compiler='gcc',
                        type_converters=weave.converters.blitz,
                        extra_compile_args=['-O3'])
            _alreadyDepressed[:,spikes] = False
            P.A_post_[spikes]=a_post[spikes] # reset A_post (~start a timer for future LTD)

        else: # all spikes
            nspikes = size(spikes)
            A_pre = mirror.A_pre_
            code = '''
            for(int si=0;si<nspikes;si++)
            {
                int i = spikes(si);
                for(int j=0;j<N;j++)
                {
                        double wnew;
                        if(mu==0){ /* additive. requires hard bound */
                            wnew = _synW(j,i)+_gmax(i)*(w_out+A_pre(j));
                            if(wnew>_gmax(i)) wnew = _gmax(i);
                            if(wnew<0) wnew = 0.0;
                        }
                        else { /* soft bound */
                            wnew = _synW(j,i)+_gmax(i)*(w_out+A_pre(j)*exp(mu*log(1-_synW(j,i)/_gmax(i))));
                            if(wnew>_gmax(i)) wnew = _gmax(i);
                            if(wnew<0) wnew = 0.0;
						}
                        _synW(j,i) = wnew;
                }
            }
            '''
            weave.inline(code,
                        ['spikes', 'nspikes', 'N', '_synW', '_gmax', 'A_pre', 'mu','w_out'],
                        compiler='gcc',
                        type_converters=weave.converters.blitz,
                        extra_compile_args=['-O3'])
            P.A_post_[spikes]+=a_post[spikes] # reset A_post (~start a timer for future LTD)

# Called whenever input neurons fire. Resets the potentials and trigger STDP updates.
# Note that mirror is a fake group that mirrors input neuron, only used for implementation issues.
# Includes C code, will be compiled the first time it is called
# Param:
# P: NeuronGroup (input neurons)
# spikes: list of the indexes of the neuron that fire.
def mirror_reset(P,spikes):
    if size(spikes):
        P.v_[spikes] = 0
        if nearestSpike:
            nspikes = size(spikes)
            A_post = neurons.A_post_
            code = '''
            for(int si=0;si<nspikes;si++)
            {
                int i = spikes(si);
                for(int j=0;j<M;j++)
                {
                    if(!_alreadyDepressed(i,j))
                    {
                        double wnew;
                        if(mu==0) { /* additive. requires hard bound */
                            wnew = _synW(i,j)+_gmax(j)*(w_in(j)+A_post(j)); 
                            if(wnew>_gmax(j)) wnew = _gmax(j);
                            if(wnew<0.0) wnew = 0.0;
                        }
                        else { /* soft bound */
                            wnew = _synW(i,j)+_gmax(j)*(w_in(j)+A_post(j)*exp(mu*log(_synW(i,j)/_gmax(j))));
                            if(wnew>_gmax(j)) wnew = _gmax(j);
                            if(wnew<0.0) wnew = 0.0;
						}
                        _synW(i,j) = wnew;
                        _alreadyDepressed(i,j) = true;
                    }
                }
            }
            '''
            weave.inline(code,
                        ['spikes', 'nspikes', 'M', '_alreadyDepressed', '_synW', '_gmax', 'A_post', 'mu','w_in'],
                        compiler='gcc',
                        type_converters=weave.converters.blitz,
                        extra_compile_args=['-O3'])
            _alreadyPotentiated[spikes,:]=False
            P.A_pre_[spikes]=a_pre  # reset A_pre (~start a timer for future LTP)
            
        else: # all spikes
            nspikes = size(spikes)
            A_post = neurons.A_post_
            code = '''
            for(int si=0;si<nspikes;si++)
            {
                int i = spikes(si);
                for(int j=0;j<M;j++)
                {
                        double wnew;
                        if(mu==0) { /* additive. requires hard bound*/
                            wnew = _synW(i,j)+_gmax(j)*(w_in(j)+A_post(j));
                            if(wnew>_gmax(j)) wnew = _gmax(j);
                            if(wnew<0.0) wnew = 0.0;
                        }
                        else { /* soft bound */
                            wnew = _synW(i,j)+_gmax(j)*(w_in(j)+A_post(j)*exp(mu*log(_synW(i,j)/_gmax(j))));
                            if(wnew>_gmax(j)) wnew = _gmax(j);
                            if(wnew<0.0) wnew = 0.0;
						}
                        _synW(i,j) = wnew;
                }
            }
            '''
            weave.inline(code,
                        ['spikes', 'nspikes', 'M', '_synW', '_gmax', 'A_post', 'mu','w_in'],
                        compiler='gcc',
                        type_converters=weave.converters.blitz,
                        extra_compile_args=['-O3'])
            P.A_pre_[spikes]+=a_pre  # reset A_pre (~start a timer for future LTP)
