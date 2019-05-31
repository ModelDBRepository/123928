import matplotlib.pyplot as plt


# Graphical plot output

printtime('*************')
printtime('* Analyzing *')
printtime('*************')

outputSelection = range(M)
#outputSelection = arange(3, 6)
#outputSelection = [83] # 42 and 82 stays at the end

# zoom inserts
nmin=194.5
nmax=205.5
tmin=-20
tmax=100

nbPlot = sum([ monitorInput, monitorInputPot, (monitorOutput or monitorPot) , monitorCurrent, monitorRate, computeOutput ])
# additional plot w = f(patternValue)
if len(outputSelection)==1 and computeOutput and os.path.exists(os.path.join('..','data','realValuedPattern.'+'%03d' % (randState)+'.mat')):
    nbPlot+=1
## additional plot w = f(index) for 'academic' patterns
#if len(outputSelection)==1 and computeOutput:
#    nbPlot+=1    
if monitorOutput and sum(a)>0 and not monitorPot:
    nbPlot += 1 # phase plot

useSubplot = True
nCol = 1 # 1
nRow = int(ceil(1.0*nbPlot/nCol))-1

## end
#minTime = max(0*second,(endTime-100/oscilFreq)*1000)
minTime = max(0*second,(endTime-3/pbTimeCompression))
maxTime = endTime

### beginning
#minTime = 0
#maxTime = 3/pbTimeCompression

## all
#minTime = 0*second
#maxTime = endTime*1000

#figure()
if useSubplot:
    fig = plt.figure()



idxPlot = 1

if monitorInput:
    # raster plots
    if useSubplot:
        idxPlotTransp = mod(idxPlot-1,nRow)*nCol + 1 + (idxPlot-1)/nRow
        labels = []
        ticks = arange(minTime,maxTime+.5,.5)
        for i in range(len(ticks)):
            labels.append(str(ticks[i]))
        ax = fig.add_subplot(nRow,nCol,idxPlotTransp,xticks=1000*ticks,xticklabels=labels,frame_on=False)
        plot([1000*minTime, 1000*maxTime],[150, 150],'k',linewidth=1.5)
        plot([1000*minTime, 1000*minTime],[150, 250],'k',linewidth=1.5)
        
#        xticks(1000*ticks,labels,size=3)

#        t=matplotlib.axis.Tick(ax, 1000*ticks, labels, size=None, gridOn=None, tick1On=True, tick2On=True, label1On=True, label2On=False, major=True)
#        t= ax.get_xticks();
#        ax.xaxis.set_major_formatter( NullFormatter() )
        
#        ax = fig.add_subplot(nRow,nCol,idxPlotTransp)
    else:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1,title='Input spikes')
        
#    raster_plot(inputSpike)

    for i in range(inputSpike.nspikes):
        if inputSpike.spikes[i][1]<minTime:
            continue
        if inputSpike.spikes[i][0]>=150 and inputSpike.spikes[i][0]<=250:
            plot([1000*inputSpike.spikes[i][1]],[inputSpike.spikes[i][0]],'.b') 
        if inputSpike.spikes[i][1]>maxTime:
            break
        
#    if useReset and maxTime - minTime < 100000:
#        plot([1000*reset,1000*reset],[-1*ones(len(reset)),N*ones(len(reset))],'-r')

#    axis([minTime, maxTime, -1, N])
    axis([1000*minTime, 1000*maxTime, 150, 250])
#    axis([1000*reset[4]-20, 1000*reset[4]+100, 180, 220])
    xlabel('Time (s)')
    ylabel('Afferent #')
#    title('A')
    text(1.02, 0.5,'A', fontsize=20, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
    idxPlot+=1
    if inputSpike.nspikes>0:
        print 'Mean input firing rate = ' + '%.1f' % (inputSpike.nspikes/(endTime-inputSpike.spikes[0][1])/N)
        
    # pattern periods     
    if os.path.exists(os.path.join('..','data','patternPeriod.'+'%03d' % (randState)+'.mat')):
        patternPeriod=loadmat(os.path.join('..','data','patternPeriod.'+'%03d' % (randState)+'.mat'))
        patternPeriod=patternPeriod['patternPeriod']/pbTimeCompression
        for i in range(size(patternPeriod,0)):
            if patternPeriod[i,1]<minTime:
                continue
            rect = Rectangle( [1000*patternPeriod[i,0], -.5], 1000*(patternPeriod[i,1]-patternPeriod[i,0]), len(realValuedPattern), facecolor='grey', edgecolor='none')
            ax.add_patch(rect)
            if patternPeriod[i,0]>maxTime:
                break
            
    print 'Input spikes done'
    
    #zoom inserts
    if useReset:
        rect = Rectangle( [1000*reset[4]+tmin, nmin], tmax, nmax-nmin, facecolor='none', edgecolor='black',linewidth=1.5)
        ax.add_patch(rect)
        rect = Rectangle( [1000*reset[7]+tmin, nmin], tmax, nmax-nmin, facecolor='none', edgecolor='black',linewidth=1.5)
        ax.add_patch(rect)
        rect = Rectangle( [1000*reset[8]+tmin, nmin], tmax, nmax-nmin, facecolor='none', edgecolor='black',linewidth=1.5)
        ax.add_patch(rect)
    else:
        rect = Rectangle( [1000*1.43+tmin, nmin], tmax, nmax-nmin, facecolor='none', edgecolor='black',linewidth=1.5)
        ax.add_patch(rect)
        rect = Rectangle( [1000*2.18+tmin, nmin], tmax, nmax-nmin, facecolor='none', edgecolor='black',linewidth=1.5)
        ax.add_patch(rect)
        rect = Rectangle( [1000*2.305+tmin, nmin], tmax, nmax-nmin, facecolor='none', edgecolor='black',linewidth=1.5)
        ax.add_patch(rect)
        
    
if monitorInputPot:
    # membrane potential
    subplot(nRow,nCol,idxPlot)
    if useSubplot:
        idxPlotTransp = mod(idxPlot-1,nRow)*nCol + 1 + (idxPlot-1)/nRow
        ax = fig.add_subplot(nRow,nCol,idxPlotTransp,xticks=[minTime,maxTime],xticklabels=[str(minTime/1000),str(maxTime/1000)])
    else:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1,xticks=[minTime,maxTime],xticklabels=[str(minTime/1000),str(maxTime/1000)])
    for j in range(N):
        plot(inputPot.times/(1*msecond),inputPot.values[j,:]/(1*mvolt))
    axis([minTime, maxTime, 1000*(El-.1*(vt-El)), 1000*(vt+.1*(vt-El))])
    xlabel('Time (s)')
    ylabel('Input - Membrane potential (in mV)')
    idxPlot+=1
    print 'Input pot. done'
if monitorCurrent:    
    # current
    if useSubplot:
        idxPlotTransp = mod(idxPlot-1,nRow)*nCol + 1 + (idxPlot-1)/nRow
        ax = fig.add_subplot(nRow,nCol,idxPlotTransp,xticks=[minTime,maxTime],xticklabels=[str(minTime/1000),str(maxTime/1000)])
    else:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1,xticks=[minTime,maxTime],xticklabels=[str(minTime/1000),str(maxTime/1000)])
    for j in outputSelection:
        plot(current.times/(1*msecond),current.values[j,:]/(1*mvolt))
    plot([0, endTime*1000],[1000*(vt-El),1000*(vt-El)],':r',linewidth=2)    
##    axis([0, 300*1000, -1*1000*(vt-El), 5.0*1000*(vt-El)])
    axis([minTime, maxTime, .0*1000*(vt-El), 2.0*1000*(vt-El)])
    xlabel('Time (s)')
    ylabel('R x input current (in mV)')
    idxPlot+=1
    print 'Current done'
if monitorPot:    
    # membrane potential
    if useSubplot:
        idxPlotTransp = mod(idxPlot-1,nRow)*nCol + 1 + (idxPlot-1)/nRow
#        ax = fig.add_subplot(nRow,nCol,idxPlotTransp,xticks=[minTime,maxTime],xticklabels=[str(minTime/1000),str(maxTime/1000)])
        ax = fig.add_subplot(nRow,nCol,idxPlotTransp,frame_on=False)
        plot([1000*minTime, 1000*maxTime],[1000*(El-.1*(vt-El)), 1000*(El-.1*(vt-El))],'k',linewidth=1.5)
        plot([1000*minTime, 1000*minTime],[1000*(El-.1*(vt-El)), 1000*(vt+.1*(vt-El))],'k',linewidth=1.5)
        
        text(1.02, 0.5,'B', fontsize=20, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
    else:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1,title='Postsynaptic potential',xticks=[minTime,maxTime],xticklabels=[str(minTime/1000),str(maxTime/1000)])
    for j in outputSelection:
        plot(pot.times/(1*second),pot.values[j,:]/(1*mvolt))
    plot([0, endTime],[1000*vt,1000*vt],':r',linewidth=1.5)    

    # pattern periods     
    if os.path.exists(os.path.join('..','data','patternPeriod.'+'%03d' % (randState)+'.mat')):
        patternPeriod=loadmat(os.path.join('..','data','patternPeriod.'+'%03d' % (randState)+'.mat'))
        patternPeriod=patternPeriod['patternPeriod']/pbTimeCompression
        for i in range(size(patternPeriod,0)):
            if patternPeriod[i,1]<=minTime-10:
                continue
            rect = Rectangle( [patternPeriod[i,0],1000*(El-.1*(vt-El))], (patternPeriod[i,1]-patternPeriod[i,0]), 1000*(1.2*(vt-El)), facecolor='grey', edgecolor='none')
            ax.add_patch(rect)
            if patternPeriod[i,0]>=maxTime+10:
                break
            
    # add vertical lines for postsynaptic spikes
    if monitorOutput: 
        tmp = reshape(outputSpike.spikes,[outputSpike.nspikes,2])[:,1]  
        plot([tmp,tmp],1000*double([1.1*El*ones(len(tmp)),0.9*vt*ones(len(tmp))]),':k')

    
    axis([minTime, maxTime, 1000*(El-.1*(vt-El)), 1000*(vt+.1*(vt-El))])
    xlabel('Time (s)')
    ylabel('Membrane potential (in mV)')
#    title('B')
    idxPlot+=1
else:
    if monitorOutput:    
        # raster plots
        if useSubplot:
            idxPlotTransp = mod(idxPlot-1,nRow)*nCol + 1 + (idxPlot-1)/nRow
            ax = fig.add_subplot(nRow,nCol,idxPlotTransp,title='B',xticks=[minTime,maxTime],xticklabels=[str(minTime/1000),str(maxTime/1000)])
        else:
            fig = plt.figure()
    #        ax = fig.add_subplot(1,1,1,title='Output spikes',xticks=[minTime,maxTime],xticklabels=[str(minTime/1000),str(maxTime/1000)])
            ax = fig.add_subplot(1,1,1,title='Output spikes')
        raster_plot(outputSpike) 
        if useReset and maxTime - minTime < 100000:
            plot([1000*reset,1000*reset],[-1*ones(len(reset)),M*ones(len(reset))],'-r')
    
       # pattern periods     
        if os.path.exists(os.path.join('..','data','patternPeriod.'+'%03d' % (randState)+'.mat')):
            patternPeriod=loadmat(os.path.join('..','data','patternPeriod.'+'%03d' % (randState)+'.mat'))
            patternPeriod=patternPeriod['patternPeriod']/pbTimeCompression
            for i in range(size(patternPeriod,0)):
                if 1000*patternPeriod[i,1]<=minTime-10000:
                    continue
                rect = Rectangle( [1000*patternPeriod[i,0], -1], 1000*(patternPeriod[i,1]-patternPeriod[i,0]), M+1, facecolor='grey', edgecolor='none')
                ax.add_patch(rect)
                if 1000*patternPeriod[i,0]>=maxTime+10000:
                    break
    
    
        axis([minTime, maxTime, min(outputSelection)-.5, max(outputSelection)+.5])
    #    axis([minTime, maxTime, 19.5, 21.5])
        xlabel('Time (s)')
        ylabel('Neuron #')
        title('B')
        idxPlot+=1
        print str(outputSpike.nspikes) + ' post synaptic spikes'
            
        print 'Output spikes done'
    #    print 'First output spikes: ' + str( outputSpike.spikes[0:min(5,outputSpike.nspikes)] )
    
        if sum(a)>0: # phases
            if useSubplot:
                idxPlotTransp = mod(idxPlot-1,nRow)*nCol + 1 + (idxPlot-1)/nRow
                ax = fig.add_subplot(nRow,nCol,idxPlotTransp)
            else:
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
            # init array
            phase = []
            times = []
            for i in range(M):
                phase.append([])
                times.append([])
            # gather data
            for i in range(outputSpike.nspikes):
                phase[outputSpike.spikes[i][0]].append( outputSpike.spikes[i][1] - floor(outputSpike.spikes[i][1]*oscilFreq)*1.0/oscilFreq )
                times[outputSpike.spikes[i][0]].append( outputSpike.spikes[i][1] )
            # plot
            for i in outputSelection:
                plot(times[i],phase[i],'.',markersize=1.5,label= str(i)+' (gmax=' + '%2.1e' % gmax[i] + ', r=' + '%.2f' % (a_post[i]/a_pre) + ')' )
    #        axis([minTime, maxTime, 0, 1/oscilFreq])
            xlabel('Time (in s)')
            xlabel('Phase (in s)')
            leg = legend(loc='lower left')
            # matplotlib.text.Text instances
            for t in leg.get_texts():
                t.set_fontsize(8)    # the legend text fontsize
            # matplotlib.lines.Line2D instances
            for l in leg.get_lines():
                l.set_markersize(3)  # the legend line width
                
            # draw pattern periods
            if os.path.exists(os.path.join('..','data','patternPeriod.'+'%03d' % (randState)+'.mat')):
                for i in range(size(patternPeriod,0)):
                    rect = Rectangle( [patternPeriod[i,0], 0], patternPeriod[i,1]-patternPeriod[i,0], 1/oscilFreq, facecolor='grey', edgecolor='none')
                    ax.add_patch(rect)
                
                
            idxPlot+=1
            print 'Output phases done'

if monitorRate:    
    # firing rates
    if useSubplot:
        idxPlotTransp = mod(idxPlot-1,nRow)*nCol + 1 + (idxPlot-1)/nRow
        ax = fig.add_subplot(nRow,nCol,idxPlotTransp)
    else:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    for i in outputSelection:
        f_pre = 10
        dwdt = f_pre*rate[i]._rate[0] * (a_pre*tau_pre/(1+tau_pre*(f_pre+rate[i]._rate[0]))+a_post[i]*tau_post/(1+tau_post*(f_pre+rate[i]._rate[0])) )
#        plot(rate[0].times/ms,rate[i].smooth_rate(6000*ms),label= str(i)+' (gmax=' + '%2.1e' % gmax[i] + ', r=' + '%.2f' % (a_post[i]/a_pre) + ', dW/dt='+'%2.1e' % dwdt)
        plot(rate[0].times/ms,rate[i]._rate,label= str(i)+' (gmax=' + '%2.1e' % gmax[i] + ', r=' + '%.2f' % (a_post[i]/a_pre) + ', dW/dt='+'%2.1e' % dwdt)
#    axis([minTime, 400*1000, 0, 90])
    xlabel('Time (in ms)')
    ylabel('Rates in Hz')
    leg = legend(loc='upper center')
    # matplotlib.text.Text instances
    for t in leg.get_texts():
        t.set_fontsize(8)    # the legend text fontsize

    # matplotlib.lines.Line2D instances
    for l in leg.get_lines():
        l.set_linewidth(1.5)  # the legend line width
    idxPlot+=1
    print 'Rates done'
if computeOutput:  
    # synaptic weights
    if useSubplot:
        idxPlotTransp = mod(idxPlot-1,nRow)*nCol + 1 + (idxPlot-1)/nRow
#        ax = fig.add_subplot(nRow,nCol,idxPlotTransp)
        ax = fig.add_subplot(3,2,5,frame_on=False)
#        text(1.02, 0.5,'C', fontsize=20, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
        text(1.04, 0.5,'C', fontsize=20, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
        plot([0,1],[0,0],'k',linewidth=1.5)
        plot([0, 0],[0, N],'k',linewidth=1.5)
        xticks([],[])
        yticks([],[])            
        text(0, -.04,'0', fontsize=12, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
        text(1, -.04,'1', fontsize=12, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
        text(-.02, 0,'0', fontsize=12, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
        text(-.02, N,'2000', fontsize=12, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
    else:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        
    # always show final synatpic weights
#    bar(histc(finalWeight,arange(0, 1.1, .1)))
    for i in outputSelection:
#        hist(finalWeight[:,i],bins=20,label= str(i)+' (gmax=' + '%2.1e' % gmax[i] + ', r=' + '%.2f' % (a_post[i]/a_pre) )
        hist(finalWeight[:,i],label= str(i)+' (gmax=' + '%2.1e' % gmax[i] + ', r=' + '%.2f' % (a_post[i]/a_pre) )
    axis([0, 1, 0, N])
#        bar(hist(finalWeight[:,i],bins=arange(0, 1.05, .05)),color=None,label= str(i)+' (gmax=' + '%2.1e' % gmax[i] + ', r=' + '%.1f' % (a_post[i]/a_pre) )
        
#    leg = legend()    
#    # matplotlib.text.Text instances
#    for t in leg.get_texts():
#        t.set_fontsize(8)    # the legend text fontsize
#    # matplotlib.lines.Line2D instances
#    for l in leg.get_lines():
#        l.set_linewidth(1.5)  # the legend line width
        
#    hist(x, bins=10, range=[0.,1.])
    xlabel('Normalized weight')
    ylabel('#')
    

    idxPlot+=1
    
    print "# of selected synapses=",str(sum(finalWeight>.5,0))    
    print "weight sum=",str(sum(finalWeight,0))    
    print 'Weight hist. done'
    
    if len(outputSelection)==1 and os.path.exists(os.path.join('..','data','realValuedPattern.'+'%03d' % (randState)+'.mat')):        
        if useSubplot:
            idxPlotTransp = mod(idxPlot-1,nRow)*nCol + 1 + (idxPlot-1)/nRow
#            ax = fig.add_subplot(nRow,nCol,idxPlotTransp)
            ax = fig.add_subplot(3,2,6,frame_on=False)
#            text(1.02, 0.5,'D', fontsize=20, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
            text(1.04, 0.5,'D', fontsize=20, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
            plot([0,1],[0,0],'k',linewidth=1.5)
            plot([0,0],[0,1],'k',linewidth=1.5)
            xticks([],[])
            yticks([],[])            
            text(0, -.04,'0', fontsize=12, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
            text(1, -.04,'1', fontsize=12, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
            text(-.02, 0,'0', fontsize=12, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
            text(-.02, 1,'1', fontsize=12, horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
        else:
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
        realValuedPattern=loadmat(os.path.join('..','data','realValuedPattern.'+'%03d' % (randState)+'.mat'))
        realValuedPattern=realValuedPattern['realValuedPattern']
        plot(realValuedPattern,finalWeight[range(len(realValuedPattern)),outputSelection],'.')
        xlabel('Pattern activation level')
        ylabel('Normalized weight')
        axis([0, 1, 0, 1])
        idxPlot+=1
    
#    if len(outputSelection)==1: # w=f(index) for academic patterns        
#        if useSubplot:
#            idxPlotTransp = mod(idxPlot-1,nRow)*nCol + 1 + (idxPlot-1)/nRow
#            ax = fig.add_subplot(nRow,nCol,idxPlotTransp,title='E',xticks=[0,1])
#        else:
#            fig = plt.figure()
#            ax = fig.add_subplot(1,1,1)
#        plot(range(N),finalWeight[range(len(realValuedPattern)),outputSelection],'.')
#        xlabel('Afferent #')
#        ylabel('Normalized weight')
#        axis([0, N-1, -.1, 1.1])
#        idxPlot+=1

#phase distrib of afferents
if sum(a)>0 and N<=108:
    # gather data
    phase = [];
    for i in range(N):
        phase.append([])
    for i in range(inputSpike.nspikes):
        phase[inputSpike.spikes[i][0]].append( inputSpike.spikes[i][1] - floor(inputSpike.spikes[i][1]*oscilFreq)*1.0/oscilFreq )
           
    fig2 = plt.figure()
    nCol2 = round(2*(3*N)**.5/3)
    nRow2 = ceil(N/nCol2)
    for i in range(N):
        ax = fig2.add_subplot(nRow2,nCol2,i+1)
        if len(phase[i])>0:
    #        b= hist(phase[i],bins=arange(0, 1/oscilFreq, 1/oscilFreq/25))
    #        print len(b[0])
    #        print len(b[1])
    #        print '-'
    #        bar(b[1],b[0])
            hist(phase[i],bins=arange(0*second, 1/oscilFreq+1/oscilFreq/25, 1/oscilFreq/25))
            axis([0*second, 1/oscilFreq, 0, 1000])
            title(str(len(phase[i])))
        if i<N-1:
            for ticklabel in ax.get_xticklabels():
                ticklabel.set_visible(False)
        if i>0:
            for ticklabel in ax.get_yticklabels():
                ticklabel.set_visible(False)

print 'Input phases done'


#for i in range(M):
#for i in []:
for i in []:
    figure()
    hist(finalWeight[:,i],bins=20)
    title( str(i)+' (gmax=' + '%2.1e' % gmax[i] + ', r=' + '%.2f' % (a_post[i]/a_pre) )
    print 'Individual hist.' + str(i) +  ' done'
    
    
# zoom inserts
# zoom 1
ax = axes([.35,.75,.1,.2],xticks=[],yticks=[])
rect = Rectangle( [1000*patternPeriod[0,0], -.5], 1000*(patternPeriod[0,1]-patternPeriod[0,0]), len(realValuedPattern), facecolor='grey', edgecolor='none')
ax.add_patch(rect)
if useReset:
    center = 1000*reset[4]
else:
    center = 1000*1.43

# grid
tmp = arange(center+tmin,center+tmax,20)
plot([tmp,tmp],double([nmin*ones(len(tmp)),nmax*ones(len(tmp))]),':k',linewidth=.75)
tmp = arange(nmin,nmax)
plot(double([(center+tmin)*ones(len(tmp)),(center+tmax)*ones(len(tmp))]),[tmp,tmp],':k',linewidth=.75)

# reset line
if useReset:
    plot([center,center],[-1,N],'-r')
    
raster_plot(inputSpike)
axis([center+tmin, center+tmax, nmin, nmax])
xlabel('')
ylabel('')

# zoom 2
ax = axes([.56,.75,.1,.2],xticks=[],yticks=[])
rect = Rectangle( [1000*patternPeriod[1,0], -.5], 1000*(patternPeriod[1,1]-patternPeriod[1,0]), len(realValuedPattern), facecolor='grey', edgecolor='none')
ax.add_patch(rect)
if useReset:
    center = 1000*reset[7]
else:
    center = 1000*2.18

# grid
tmp = arange(center+tmin,center+tmax,20)
plot([tmp,tmp],double([nmin*ones(len(tmp)),nmax*ones(len(tmp))]),':k',linewidth=.75)
tmp = arange(nmin,nmax)
plot(double([(center+tmin)*ones(len(tmp)),(center+tmax)*ones(len(tmp))]),[tmp,tmp],':k',linewidth=.75)

# reset line
if useReset:
    plot([center,center],[-1,N],'-r')
    
raster_plot(inputSpike)
axis([center+tmin, center+tmax, nmin, nmax])
xlabel('')
ylabel('')
    
# zoom 3
ax = axes([.79,.75,.1,.2],xticks=[],yticks=[])
if useReset:
    center = 1000*reset[8]
    rect = Rectangle( [1000*patternPeriod[1,0], -.5], 1000*(patternPeriod[1,1]-patternPeriod[1,0]), len(realValuedPattern), facecolor='grey', edgecolor='none')
else:
    center = 1000*2.305
    rect = Rectangle( [1000*patternPeriod[1,0], -.5], 1000*(patternPeriod[2,1]-patternPeriod[1,0]), len(realValuedPattern), facecolor='grey', edgecolor='none')
ax.add_patch(rect)

# grid
tmp = arange(center+tmin,center+tmax,20)
plot([tmp,tmp],double([nmin*ones(len(tmp)),nmax*ones(len(tmp))]),':k',linewidth=.75)
tmp = arange(nmin,nmax)
plot(double([(center+tmin)*ones(len(tmp)),(center+tmax)*ones(len(tmp))]),[tmp,tmp],':k',linewidth=.75)

# reset line
if useReset:
    plot([center,center],[-1,N],'-r')
    
raster_plot(inputSpike)
axis([center+tmin, center+tmax, nmin, nmax])
xlabel('')
ylabel('')

show()

# for i in range(N-1):
#    print pot.mean[i]
#     print pot[i][len(pot[i])-1]
#     print pot[i][0]
#     print len(pot[i])
# print pot.times[len(pot.times)-1]

# write output in a file
# f = open("spikeList.txt", 'w')
# f.write("neuron spikeTime\n")
# for i in range(spike.nspikes-1):
#     f.write(str(spike.spikes[i][0])+" "+str(spike.spikes[i][1])+"\n")
# f.close()
