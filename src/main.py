# Main file.

execfile('init.py') # sets the parameters

#****************************************************
printtime('************************')
printtime('* Computing (rand=' + '%03d' % randState + ') *')
printtime('************************')


if nearestSpike: # flags useful for nearest spike mode
    alreadyDepressed = zeros([N,M],dtype=int)
    alreadyPotentiated = zeros([N,M],dtype=int)

# for memory issues the simulation is divided in multiple periods
# count is the current period number
count = 0

if computeOutput: # instantiate output neurons
    
    _gmax = asarray(gmax) # faster implementation (does not verify dimensions)
    
    if nearestSpike: # flags useful for nearest spike mode
        _alreadyPotentiated = asarray(alreadyPotentiated)
        _alreadyDepressed = asarray(alreadyDepressed)
        
    # mirrors
    mirror_eqs=''' 
    v:1
    dA_pre/dt=-A_pre/tau_pre : 1
    ''' 
    mirror = NeuronGroup(N, model=mirror_eqs, threshold=0.5, reset=mirror_reset)
      
    #STDP neuron
    if conductanceOutput:
        eqs_neurons='''
        dv/dt=(ge*(Ee-v)+El-v)/taum + sigma*xi/taum**.5 : volt
        dge/dt=-ge/taue : 1
        dA_post/dt=-A_post/tau_post : 1
        '''        
    else:
        eqs_neurons='''
        dv/dt=(ge+El-v)/taum + sigma*xi/taum**.5 : volt
        dge/dt=-ge/taue : volt
        dA_post/dt=-A_post/tau_post : 1
        '''
    
    neurons_cr = CustomRefractoriness(neurons_reset,period=refractoryPeriod,state='v')
    
    if poissonOutput: # stochastic spike generation
        neurons=NeuronGroup(M,model=eqs_neurons,threshold=PoissonThreshold(),reset=neurons_cr)
    else: # deterministic spike generation
        neurons=NeuronGroup(M,model=eqs_neurons,threshold=vt,reset=neurons_cr)    
   
    #connections
    synapses=Connection(mirror,neurons,'ge',structure='dense')
    seed(randState)
    if useSavedWeight and os.path.exists(os.path.join('..','data','weight.'+'%03d' % (randState)+'.mat')):
        print 'Loading previously dumped weight'    
        tmp=loadmat(os.path.join('..','data','weight.'+'%03d' % (randState)+'.mat'))
        tmp=tmp['weight']     
        initialWeight = zeros([N,M])
        if M>1:
            for j in range(M):
                initialWeight[:,j] = tmp[:,j]*gmax[j]
        else:
            initialWeight[:,0] = tmp[:]*gmax[0]
        del tmp
    else: # start from random synaptic weights
        initialWeight = zeros([N,M])
        for i in range(N):
            initialWeight[i,:] = initialWeight_min + rand(1)*(initialWeight_max-initialWeight_min)
        if initialWeight.max() > min(gmax):
            print '***********************************************************'
            print '* WARNING: Initial weight > gmax. This should not happen. *'
            print '***********************************************************'
    synapses.connect(mirror,neurons,initialWeight)
    synapses.compress() 
    _synW = asarray(synapses.W)

    # affect initial values
    neurons.v_ = vr+rand(1)*ones(len(neurons))*(vt-vr)
    
    neurons.A_post_=0*volt
    neurons.ge_=0*volt
    mirror.A_pre_=0*volt
     
startTime = timeOffset;

if recomputeSpikeList: # input layer is computed. Need to instantiate corresponding neurons
    if poissonInput:
        input=PoissonGroup(N,zeros(N)) 
    else:    
        if sum(a)==0*namp:
            input_eqs=''' 
            dv/dt = (-v + El + R*I )/taum + sigma*xi/taum**.5 : volt
            I : amp
            ''' 
        else:
            # sinusoidal oscillation
            input_eqs=''' 
            dv/dt = (-v + El + R * ( I + aa*(sin(t*2*pi*oscilFreq-pi))))/taum + sigma*xi/taum**.5 : volt
            I : amp
            aa : amp
            ''' 
#            # sawtooth oscillation
#            input_eqs=''' 
#            dv/dt = (-v + El + R * ( I + a*((t-floor(t*oscilFreq)*1/oscilFreq)*oscilFreq-1.0)))/taum + sigma*xi/taum**.5 : volt
#            I : amp
#            ''' 
        input=NeuronGroup(N,model=input_eqs,threshold=vt,reset=vr,refractory=refractoryPeriod) 
        # affect initial values
        seed(randState)
        input.v_ = vr+rand(N)*(vt-vr)
#        input.v_ = El
        input.aa_ = a # amplitude of the oscillatory input current

    if computeOutput: # connect input layer to its mirror
        C_input_mirror = IdentityConnection(input, mirror)
        
        
    if useReset: # load reset times from reset.###.mat
        reset=loadmat(os.path.join('..','data','reset.'+'%03d' % (randState)+'.mat'))
        reset=reset['resetTimes']
#        # to artificially double the reset frequency:
#        reset=.5*concatenate([reset,reset[-1]+250e-3+reset])
    else:
        reset=[Inf]
    rcursor=0 # cursor for reset array
    
    
    I = zeros(N)*namp # just to initialize the array
    
    
    # open input current file
    f = open(os.path.join('..','data','inputValues.'+'%03d' % (randState)+'.txt'),'r')
    
    printtime('Starting (recompute spike list)')
    
    localStartTime = time()*second

    inputSpike = SpikeMonitor(input,True) # input spikes need to be monitored all the time with this kind of computation

    for l in f: # iterate on file lines
    
        # read time
        endTime = double(l[0:9])*second/pbTimeCompression

        # imposed end time
        endTime = min(imposedEnd,endTime)
    
        # monitors
        if endTime>=monitorTime and not isMonitoring:
            print '********************'            
            print '* Start monitoring *'            
            print '********************'            
            isMonitoring = True
#            if monitorInput:
#                inputSpike = SpikeMonitor(input,True)
            if monitorInputPot:
                inputPot = StateMonitor(input,'v',record=True)
            if monitorOutput:
                outputSpike = SpikeMonitor(neurons,True)
            if monitorPot:
                pot = StateMonitor(neurons,'v',record=True)
            if monitorCurrent:
                current = StateMonitor(neurons,'ge',record=True)
            if monitorRate:
                rate = []
                for i in range(M):
                    rate.append(PopulationRateMonitor(neurons[i],bin=10000*ms))
    
        # read input currents
        for i in range(N):
            I[i]= Imin[i] + double(l[11+i*8:11+(i+1)*8-2]) * (Imax[i]-Imin[i])
#            I[i]= Imin[i] + i*1.0/(N-1) * (Imax[i]-Imin[i])
#            I[i]= Imin[i] + mod(i,12)*1.0/(N/9-1) * (Imax[i]-Imin[i])
                     
        # affect currents
        if poissonInput: # with poisson currents in fact correspond to rates
            input._S[0,:] = I
        else:
            input.I_ = I
    
        if reset[rcursor]> endTime: # no reset during period
            defaultclock.reinit(startTime) # make sure end time is exactly the one we want, to avoid drifting
            run(endTime-startTime) # run Brian simulator until endTime
        else:
            fromTime = startTime
            while reset[rcursor]<= endTime: # iterate on resets
                defaultclock.reinit(fromTime) # make sure end time is exactly the one we want, to avoid drifting
                run(reset[rcursor]-fromTime) # run Brian until reset
#                input.v_ = El # reset to resting potential
                input.v_ = vr # reset to reset potential
                
                fromTime = reset[rcursor]
                rcursor = rcursor+1
                
            # last bit    
            defaultclock.reinit(reset[rcursor-1]) # make sure end time is exactly the one we want, to avoid drifting
            run(endTime-reset[rcursor-1]) # run Brian until end time
            
    
        if endTime>=imposedEnd: # exit condition
            break
    
        if inputSpike.nspikes > 1000000: # periodic log output
            printtime('Period # '+ str(count+1) +' - simulated time: '+str(endTime)+' - computation time: ' + str(time()*second-localStartTime))
            localStartTime = time()*second
            if dumpSpikeList: # dump spike list of this period
                print 'Dumping #' + str(count+1) + ' (nspikes=' + str(inputSpike.nspikes) + ')'
                savemat(os.path.join('..','data','spikeList.'+'%03d' % (randState)+'.'+ '%03d' % (count+1) +'.mat'),{'sl':inputSpike.spikes})
                print 'Dumping time: '+ str(time()*second-localStartTime)                        
                localStartTime = time()*second
            inputSpike.reinit()
            count += 1
                
        # periodic graphic plot output
        if floor(endTime/analyzePeriod)!=floor(startTime/analyzePeriod):
            # compute final normalized weight (there's probably a smarter way to do that...)
            if computeOutput:  
                finalWeight = zeros([N,M])
                for i in range(N):
                    for j in range(M):
                        finalWeight[i,j] = _synW[i,j]/gmax[j]
            execfile('analyze.py')
                                        
        # start := end
        startTime = endTime
            
    f.close() # close input current file
    
    # last dump
    if inputSpike.nspikes>0: # something to dump
        printtime('Period # '+ str(count+1) +' - simulated time: '+str(endTime)+' - computation time: ' + str(time()*second-localStartTime))
        localStartTime = time()*second
        if dumpSpikeList:
            print 'Dumping #' + str(count+1) + ' (nspikes=' + str(inputSpike.nspikes) + ')'
            savemat(os.path.join('..','data','spikeList.'+'%03d' % (randState)+'.'+ '%03d' % (count+1) +'.mat'),{'sl':inputSpike.spikes})
            print 'Dumping time: ' + str(time()*second-localStartTime)                        
            localStartTime = time()*second
            
else: # not in recomputeSpikeList mode. Thus spikes are read from files spikeList.###.###.mat (first number: random seed, second number: file number) 
    
    if not computeOutput:
        print 'Warning: bad configuration: compute neither output nor output...'
    
                
    
    printtime('Starting (use saved spike list)')
    
    # look for spike list files
    fileList = listMatFile('../data/',randState)
    print str(len(fileList)) + ' spike list files found'
    
   
    for fl in fileList: # iterate on spile list files

        # read spike list
        localStartTime = time()*second
        print 'Reading '+ fl
        spikeList=loadmat(os.path.join('..','data',fl))
        spikeList=spikeList['sl']
        spikeList[:,1]+=timeOffset
        spikeList[:,1]/=pbTimeCompression
        print str(size(spikeList,0)) + ' spikes read (in ' + str(time()*second-localStartTime) + ')'
        
        input = SpikeGeneratorGroup(N, spikeList) # special Brian NeuronGroup that fire at specified dates
        endTime = spikeList[-1][1]       
        del spikeList
        
        # monitors
        if endTime>=monitorTime and not isMonitoring:
            print '********************'            
            print '* Start monitoring *'            
            print '********************'            
            isMonitoring = True
            if monitorInput:
                inputSpike = SpikeMonitor(input,True)
            if monitorOutput:
                outputSpike = SpikeMonitor(neurons,True)
            if monitorPot:
                pot = StateMonitor(neurons,'v',record=True)
            if monitorCurrent:
                current = StateMonitor(neurons,'ge',record=True)
            if monitorRate:
                rate = []
                for i in range(M):
                    rate.append(PopulationRateMonitor(neurons[i],bin=2000*ms))
                                
        # imposed end time
        endTime = min(imposedEnd,endTime)
        
        # connect new spike generator
        C_input_mirror = IdentityConnection(input, mirror)
        

        # run
        print 'Running from t=' + str(startTime) + ' to t=' + str(endTime)
        defaultclock.reinit(startTime) # make sure end time is exactly the one we want, to avoid drifting
        run(endTime-startTime) # run Brian simulator until endTime
    
        # periodic graphic plot output
        if floor(endTime/analyzePeriod)!=floor(startTime/analyzePeriod):
            # compute final normalized weight (there's probably a smarter way to do that...)
            if computeOutput:  
                finalWeight = zeros([N,M])
                for i in range(N):
                    for j in range(M):
                        finalWeight[i,j] = _synW[i,j]/gmax[j]
            execfile('analyze.py')

        # start := end
        startTime = endTime

        # explicitly free memory
        del input

        printtime('Period # '+ str(count+1) +': computation time: ' + str(time()*second-localStartTime))
        localStartTime = time()*second
        count += 1
        
        if endTime>=imposedEnd:
            break
        
        for j in range(M):
            if mean(_synW[:,j])/gmax[j]>burstingCriterion:
                print 'WARNING: neuron # ' + str(j) + ' is bursting. Disconnecting it.'
                _synW[:,j] = 0*mV
                gmax[j] = 0*mV


        

print 'Total computation time: ' + str(time()*second-globalStartTime)

# compute final normalized weight (there's probably a smarter way to do that...)
if computeOutput:  
    finalWeight = zeros([N,M])
    for i in range(N):
        for j in range(M):
            finalWeight[i,j] = _synW[i,j]/gmax[j]


#execfile('pickleAll.py') # pickle all variable (under development...)
if imposedEnd>6: # don't dump short simulations, probably done for display only
    execfile('saveWeight.py') # dump final weights 
    
#execfile('toMatlab.py') # dump variables in a mat file (under development)

if graph: # graphical plot output
    execfile('analyze.py')
    show()

if imposedEnd>6 and monitorOutput: # mutual info (stimulus, response)
    execfile('mutualInfo.py')
    show()


