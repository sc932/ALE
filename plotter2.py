#!/usr/bin/python

import matplotlib.pylab as plt
import numpy
from Bio.Seq import Seq
from Bio import SeqIO
import scipy
import time, sys, math

__usage__ = "(C) Scott Clark 2010.\nUsage:\npython plotter2.py [placementFile.ale]\nWhere this is the file outputed by ALE\n"

#class Library():
    #def __init__(self):
        #self.inwardInsertLen = {0:0}
        #self.outwardInsertLen = {0:0}
        #self.placementsPerRead = {0:0,1:0}
    #def printMe(self, start = 0, end = 0):
        #maxValIn = max(libraries[i].inwardInsertLen.keys())
        #maxValOut = max(libraries[i].outwardInsertLen.keys())
        #if end == 0:
            #end = int(numpy.max(maxValIn, maxValOut))
        #distIn = numpy.zeros(maxValIn+1)
        #distOut = numpy.zeros(maxValOut+1)
        #for j in libraries[i].inwardInsertLen.keys():
            #distIn[j] = libraries[i].inwardInsertLen[j]
        #for j in libraries[i].outwardInsertLen.keys():
            #distOut[j] = libraries[i].outwardInsertLen[j]
        #fig = plt.figure()
        #plt.plot(distIn[start:end])
        #plt.plot(distOut[start:end])

class Contig():
    def __init__(self, length, name = "name"):
        self.name = name
        self.length = length
        self.depth = numpy.ones(length)
        self.depthProb = numpy.zeros(length)
        self.placementProb = numpy.zeros(length)
        self.kmerProb = numpy.zeros(length)
        self.totalProb = numpy.zeros(length)
        self.GCcorrection = numpy.ones(length)
        self.depthProbNormalizer = 1.0
        self.placementProbNormalizer = 1.0
        self.kmerProbNormalizer = 1.0
        self.weightsDepth = numpy.ones(length)
        self.weightsPlace = numpy.ones(length)
        self.weightsKmer = numpy.ones(length)
        self.avgProbNorm = 10.0
        self.GCcont = numpy.zeros(length)
    def prePlot(self, infile, window=77):
        self.getGCcont(infile, window)
        #self.GCcorrect()
        self.normalizeMe()
        #self.GCnormalizeMe()
    def getGCcont(self, infile, window=77):
        #infile is the fasta for the assembly
        self.GCcont = GCcontentFunc(infile, window)
    def GCcorrect(self):
        # this is broken right now
        self.GCcorrection = GCvsDepthPlot(self.GCcont, self.depth)
    def plotDepth(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.depth*self.GCcorrection)
    def plotGCandDepthHistogram(self):
        numLevels = 10
        GClevels = numpy.zeros((numLevels, numpy.max(self.depth)))
        for i in range(len(self.GCcont)):
            GClevels[math.floor(self.GCcont[i]*10.0)][self.depth[i]] += 1.0
        leg = []
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i in range(numLevels):
            ax.plot(GClevels[i][:])
            leg.append(str(i*10) + "% - " + str((i+1)*10))
        ax.legend(leg, 'upper right')
        plt.show()
    def GCnormalizeMe(self):
        self.GCcorrection = smooth(self.GCcorrection, 5000)
        for i in range(self.length):
            self.weightsDepth[i] = 1.0/self.GCcorrection[i]
            self.weightsKmer[i] += 1.0 - 0.5/self.GCcorrection[i]
            self.weightsPlace[i] += 1.0 - 0.5/self.GCcorrection[i]
    def normalizeMe(self):
        self.placementProb = (self.placementProb > -60)*self.placementProb - 60*(self.placementProb <= -60)
        self.depthProb = (self.depthProb > -60)*self.depthProb - 60*(self.depthProb <= -60)
        self.kmerProb = (self.kmerProb > -60)*self.kmerProb - 60*(self.kmerProb <= -60)
        print "type, min, max, mean, std"
        print "place", numpy.min(self.placementProb), numpy.max(self.placementProb), numpy.mean(self.placementProb), numpy.std(self.placementProb)
        print "place", numpy.min(self.depthProb), numpy.max(self.depthProb), numpy.mean(self.depthProb), numpy.std(self.depthProb)
        print "place", numpy.min(self.kmerProb), numpy.max(self.kmerProb), numpy.mean(self.kmerProb), numpy.std(self.kmerProb)

        depthTotal = numpy.sum(self.depthProb)
        placementTotal = numpy.sum(self.placementProb)
        kmerTotal = numpy.sum(self.kmerProb)
        self.depthProbNormalizer = -self.avgProbNorm*float(self.length)/depthTotal
        self.placementProbNormalizer = -self.avgProbNorm*float(self.length)/placementTotal
        self.kmerProbNormalizer = -self.avgProbNorm*float(self.length)/kmerTotal


    def thresholdMe(self, threshold = -60, minLen = 1000):
        print "Computing thresholds..."
        ALEcontigs = []
        lens = []
        passed = self.totalProb > threshold
        start = 0
        bar = progressBar(0, self.length, 42)
        for i in range(len(passed)):
            if (i)%(self.length/40) == 0:
                bar(i)
            if passed[i] == True and i > 0 and passed[i-1] == False:
                start = i
            elif passed[i] == False and i > 0 and passed[i-1] == True:
                end = i
                if end - start >= minLen:
                    ALEcontigs.append([start,end])
                    lens.append(end-start)
                start = end
        print ""
        return ALEcontigs, lens
    def plotMe(self, start = 0, end = 0, plotType = "dpkt", saver = False, smoothLevelD = 10000, smoothLevelP = 1000, smoothLevelK = 1000):
        if end < 1 or end < start or end > self.length:
            end = self.length
        fig = plt.figure()
        ax2 = fig.add_subplot(111)
        legendList = []
        if "d" in plotType:
            ax2.plot(smooth(self.depthProbNormalizer*self.depthProb[numpy.max(0,start-smoothLevelD):numpy.min(self.length-1,end+smoothLevelD)], smoothLevelD)[start:end], 'r')
            legendList.append('Average Depth Likelihood')
        if "p" in plotType:
            ax2.plot(smooth(self.placementProbNormalizer*self.placementProb[numpy.max(0,start-smoothLevelP):numpy.min(self.length-1,end+smoothLevelP)], smoothLevelP)[start:end], 'b')
            legendList.append('Average Placement Likelihood')
        if "k" in plotType:
            ax2.plot(smooth(self.kmerProbNormalizer*self.kmerProb[numpy.max(0,start-smoothLevelK):numpy.min(self.length-1,end+smoothLevelK)], smoothLevelK)[start:end], 'g')
            legendList.append('Average k-mer Likelihood')
        if "t" in plotType:
            ax2.plot(smooth(self.depthProbNormalizer*self.depthProb[numpy.max(0,start-smoothLevelD):numpy.min(self.length-1,end+smoothLevelD)]*self.weightsDepth[numpy.max(0,start-smoothLevelD):numpy.min(self.length-1,end+smoothLevelD)], smoothLevelD)[start:end] + smooth(self.placementProbNormalizer*self.placementProb[numpy.max(0,start-smoothLevelP):numpy.min(self.length-1,end+smoothLevelP)]*self.weightsPlace[numpy.max(0,start-smoothLevelP):numpy.min(self.length-1,end+smoothLevelP)], smoothLevelP)[start:end] + smooth(self.kmerProbNormalizer*self.kmerProb[numpy.max(0,start-smoothLevelK):numpy.min(self.length-1,end+smoothLevelK)]*self.weightsKmer[numpy.max(0,start-smoothLevelK):numpy.min(self.length-1,end+smoothLevelK)], smoothLevelK)[start:end], 'm')
            legendList.append('Average total Likelihood')
        ax2.legend(legendList, 'center left')
        ax2.plot([0, end-start], [-10, -10], 'b--')
        ax2.plot([0, end-start], [-30, -30], 'b--')
        ax2.plot([0, end-start], [-12.5, -12.5], 'r--')
        ax2.plot([0, end-start], [-35, -35], 'r--')
        ax2.plot([0, end-start], [-7.5, -7.5], 'r--')
        ax2.plot([0, end-start], [-25, -25], 'r--')
        ax2.set_title('Average Likelihoods')
        ax2.set_xlabel('Position')
        ax2.set_ylabel('Avg (log) Likelihood')
        ax2.set_xlim((0, end-start))
        #newXtick = []
        #partitions = len(ax2.get_xticklabels())
        #for i in range(partitions):
            #newXtick.append(start + i*(end-start)/(partitions-1))
        #ax2.set_xticklabels(newXtick)
        #plt.plot(numpy.log(contigs[i].placementProb/(contigs[i].inDepth + contigs[i].outDepth)))
        plt.show()
        if saver == True:
            plt.savefig(self.name + '.png')
    def gifPlotMe(self, start = 174600, end = 175300, plotType = "dpkt", saver = False, smoothLevelD = 1, smoothLevelP = 1, smoothLevelK = 100, namer = "gif1"):
        self.normalizeMe()
        if end < 1 or end < start or end > self.length:
            end = self.length
        fig = plt.figure()
        ax2 = fig.add_subplot(111)
        legendList = []
        if "d" in plotType:
            ax2.plot(smooth(self.depthProbNormalizer*self.depthProb[numpy.max(0,start-smoothLevelD):numpy.min(self.length-1,end+smoothLevelD)], smoothLevelD)[start:end], 'r')
            legendList.append('Depth Like.')
        if "p" in plotType:
            ax2.plot(smooth(self.placementProbNormalizer*self.placementProb[numpy.max(0,start-smoothLevelP):numpy.min(self.length-1,end+smoothLevelP)], smoothLevelP)[start:end], 'b')
            legendList.append('Placement Like.')
        if "k" in plotType:
            ax2.plot(smooth(self.kmerProbNormalizer*self.kmerProb[numpy.max(0,start-smoothLevelK):numpy.min(self.length-1,end+smoothLevelK)], smoothLevelK)[start:end], 'g')
            legendList.append('k-mer Like.')
        if "t" in plotType:
            ax2.plot(smooth(self.depthProbNormalizer*self.depthProb[numpy.max(0,start-smoothLevelD):numpy.min(self.length-1,end+smoothLevelD)]*self.weightsDepth[numpy.max(0,start-smoothLevelD):numpy.min(self.length-1,end+smoothLevelD)], smoothLevelD)[start:end] + smooth(self.placementProbNormalizer*self.placementProb[numpy.max(0,start-smoothLevelP):numpy.min(self.length-1,end+smoothLevelP)]*self.weightsPlace[numpy.max(0,start-smoothLevelP):numpy.min(self.length-1,end+smoothLevelP)], smoothLevelP)[start:end] + smooth(self.kmerProbNormalizer*self.kmerProb[numpy.max(0,start-smoothLevelK):numpy.min(self.length-1,end+smoothLevelK)]*self.weightsKmer[numpy.max(0,start-smoothLevelK):numpy.min(self.length-1,end+smoothLevelK)], smoothLevelK)[start:end], 'm')
            legendList.append('Total Like.')
        ax2.legend(legendList, 'lower left')
        ax2.plot([0, end-start], [-10, -10], 'b--')
        ax2.plot([0, end-start], [-30, -30], 'b--')
        ax2.plot([0, end-start], [-12.5, -12.5], 'r--')
        ax2.plot([0, end-start], [-35, -35], 'r--')
        ax2.plot([0, end-start], [-7.5, -7.5], 'r--')
        ax2.plot([0, end-start], [-25, -25], 'r--')
        ax2.set_title('Average Likelihoods')
        ax2.set_xlabel('Position')
        ax2.set_ylabel('Avg (log) Likelihood')
        ax2.set_xlim((0, end-start))
        ax2.set_ylim((-350, 0))
        #newXtick = []
        #partitions = len(ax2.get_xticklabels())
        #for i in range(partitions):
            #newXtick.append(start + i*(end-start)/(partitions-1))
        #ax2.set_xticklabels(newXtick)
        #plt.plot(numpy.log(contigs[i].placementProb/(contigs[i].inDepth + contigs[i].outDepth)))
        #plt.show()
        plt.savefig(namer + '.png')
    def scoreMe(self, start = 0, end = 0):
        if end < 1 or end < start or end > self.length:
            end = self.length
        scoreVec = numpy.zeros(end-start)
        print "Scoring contig " + self.name + "..."
        bar = progressBar(0, end-start, 42)
        for i in range(start,end):
            if (i - start)%(end/40) == 0:
                bar(i - start)
            if not numpy.isnan(self.depthProb[i]):
                scoreVec[i-start] += self.depthProb[i]
            else:
                scoreVec[i-start] += -20.0
            if not numpy.isnan(self.placementProb[i]):
                scoreVec[i-start] += self.placementProb[i]
            else:
                scoreVec[i-start] += -20.0
            if not numpy.isnan(self.kmerProb[i]):
                scoreVec[i-start] += self.kmerProb[i]
            else:
                scoreVec[i-start] += -20.0
            if numpy.isnan(scoreVec[i-start]) or scoreVec[i-start] < -60:
                scoreVec[i-start] = -60
        score = (1.0/(float(end - start)))*numpy.sum(scoreVec)
        bar(end-start)
        print ""
        return score
    
def gifScore():
    scores = []
    for i in range(11):
        filein = "gifMaker/Ecoli_first350k_e" + str(i) + ".ale"
        ecoliE0 = ReadInInfo(filein)
        ecoliE0[0].normalizeMe()
        scores.append(ecoliE0[0].scoreMe())
        print scores[-1]
    return scores

def smooth2(x, y):
    return x

def plotGCandDepthHistogram(depth, GCcont):
    numLevels = 10
    chunkSize = 0.001
    DepthVsGC = numpy.zeros(int(1/chunkSize))
    GClevels = numpy.zeros((numLevels, numpy.max(depth)+1))
    for i in range(len(GCcont)):
        GClevels[math.floor(GCcont[i]*10.0)][depth[i]] += 1.0
        DepthVsGC[math.floor(GCcont[i]/chunkSize)] += 1.0
    leg = []
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(1, numLevels-2):
        ax.plot(GClevels[i][1:])
        leg.append(str(i*10) + "% - " + str((i+1)*10) + "%")
    ax.legend(leg, 'upper right')
    plt.xlabel("Depth")
    plt.ylabel("Number of observations of that depth")
    plt.title("Distribution of depths, for different 'bins' of average GC content")
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #xvals = numpy.arange(0,100, 0.1)
    #print numpy.shape(xvals)
    #print numpy.shape(DepthVsGC)
    #ax.plot(xvals, DepthVsGC, 'r*')
    #plt.xlabel("GC%")
    #plt.ylabel("Depth")
    #plt.title("Depth vs avg GC content over windows of size 5000")
    plt.show()

def ScoreContigsSection(c, start, end):
    print "Scoring all contigs in section (" + str(start) + "," + str(end) + "):"
    for i in range(len(c)):
        print "Contig " + str(i+1) + " score in section (" + str(start) + "," + str(end) + "):" + str(c[i].scoreMe(start=start, end=end))
        c[i].plotMe(start=start, end=end)

def ScoreContigs(c):
    print "Scoring all contigs:"
    for i in range(len(c)):
        print "Contig " + str(i+1) + " score:" + str(c[i].scoreMe())
        c[i].plotMe()

def ReadInInfo(placementFile):
    fPl = open(placementFile, 'r')

    contigs = []

    for line in fPl:
        if line[0] == '#':
            print ""
            try:
              tName = line.split(':')[1].split(' ')[0]
            except IndexError:
              tName = ""
            if tName != "":
              tLen = int(line.split(' ')[3])
              contigs.append(Contig(tLen, name = tName))
              place = 0
              print "Reading in contig: " + tName + " len " + str(tLen)
              bar = progressBar(0, tLen, 42)
              fPl.readline() # get past info
        else:
            data = line.split(' ')
            contigs[-1].depth[place] = numpy.double(data[0])
            for i in range(1,5):
                if "-nan"==data[i] or "nan"==data[i] or "inf"==data[i] or "-inf"==data[i] or numpy.double(data[i]) != numpy.double(data[i]):
                    data[i] = -60.0
            contigs[-1].depthProb[place] = numpy.double(data[1])
            contigs[-1].placementProb[place] = numpy.double(data[2])
            contigs[-1].kmerProb[place] = numpy.double(data[3])
            contigs[-1].totalProb[place] = numpy.double(data[4])
            place += 1
            if (place)%(tLen/40) == 0:
                bar(place)
    print ""
    return contigs

#//uses Stirlings approximation to high precision
#double lnfact(double input){
    #return (input - 0.5)*log(input) - input + lnfactconst - 1.0/(12.0*input) - 1.0/(360.0*input*input*input) - 1.0/(1260.0*input*input*input*input*input);
#}

def poissonPMF(val, mean):
    P = -mean + val*numpy.log(mean) - lnfact(val)
    return P

def lnfact(ins):
    return ins*numpy.log(ins) - ins + 0.5*numpy.log(2*numpy.pi*ins) + 1.0/(12.0*ins) - 1.0/(360.0*ins*ins*ins) + 1.0/(1260.0*ins*ins*ins*ins*ins)
    #lnfactconst = 0.918938533204672741780329
    #return (ins - 0.5)*numpy.log(ins) - ins + lnfactconst - 1.0/(12.0*ins) - 1.0/(360.0*ins*ins*ins) - 1.0/(1260.0*ins*ins*ins*ins*ins)

# from: http://code.activestate.com/recipes/168639-progress-bar-class/
class progressBar:
    """ Creates a text-based progress bar. Call the object with the `print'
        command to see the progress bar, which looks something like this:

        [=======>        22%                  ]

        You may specify the progress bar's width, min and max values on init.
    """

    def __init__(self, minValue = 0, maxValue = 100, totalWidth=80):
        self.progBar = "[]"   # This holds the progress bar string
        self.min = minValue
        self.max = maxValue
        self.span = maxValue - minValue
        self.width = totalWidth
        self.amount = 0       # When amount == max, we are 100% done
        self.updateAmount(0)  # Build progress bar string

    def updateAmount(self, newAmount = 0):
        """ Update the progress bar with the new amount (with min and max
            values set at initialization; if it is over or under, it takes the
            min or max value as a default. """
        if newAmount < self.min: newAmount = self.min
        if newAmount > self.max: newAmount = self.max
        self.amount = newAmount

        # Figure out the new percent done, round to an integer
        diffFromMin = float(self.amount - self.min)
        percentDone = (diffFromMin / float(self.span)) * 100.0
        percentDone = int(round(percentDone))

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2
        numHashes = (percentDone / 100.0) * allFull
        numHashes = int(round(numHashes))

        # Build a progress bar with an arrow of equal signs; special cases for
        # empty and full
        if numHashes == 0:
            self.progBar = "[>%s]" % (' '*(allFull-1))
        elif numHashes == allFull:
            self.progBar = "[%s]" % ('='*allFull)
        else:
            self.progBar = "[%s>%s]" % ('='*(numHashes-1),
                                        ' '*(allFull-numHashes))

        # figure out where to put the percentage, roughly centered
        percentPlace = (len(self.progBar) / 2) - len(str(percentDone))
        percentString = str(percentDone) + "%"

        # slice the percentage into the bar
        self.progBar = ''.join([self.progBar[0:percentPlace], percentString,
                                self.progBar[percentPlace+len(percentString):]
                                ])

    def __str__(self):
        return str(self.progBar)

    def __call__(self, value):
        """ Updates the amount, and writes to stdout. Prints a carriage return
            first, so it will overwrite the current line in stdout."""
        print '\r',
        self.updateAmount(value)
        sys.stdout.write(str(self))
        sys.stdout.flush()

def smooth(x,window_len=11,window='flat'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also: 

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string   
    """
    #t0 = time.time()
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=numpy.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    #print "Smooth of len " + str(window_len) + " in " + str(time.time() - t0)
    return y[window_len:-window_len+1]

def smooth2(x, window):
    vec = x
    for i in range(1, window):
        vec[:-i] += x[i:]
    return vec/float(window)

def smooth3(x, window):
    t0 = time.time()
    window = (window/2)*2 # force window to be a multiple of 2
    vec = numpy.zeros(len(x))
    prev = numpy.zeros(window)
    prev[0] = numpy.sum(x[0:window/2])
    fwindow = float(window)
    for i in range(window/2):
        vec[i] = prev[i%window]/float(i + window/2)
        prev[(i+1)%window] = prev[(i)%window] + x[i+window/2]
    for i in range(window/2, len(x)-window/2):
        vec[i] = prev[i%window]/fwindow
        prev[(i+1)%window] = prev[(i)%window] + x[i+window/2]
        prev[(i+1)%window] -= x[i-window/2]
    for i in range(len(x)-window/2, len(x)):
        vec[i] = prev[i%window]/(float(len(x) - i))
        prev[(i+1)%window] = prev[(i)%window] - x[i-window/2]
    print "Smooth of len " + str(window) + " in " + str(time.time() - t0)
    return vec

def smooth4(x, window):
    t0 = time.time()
    window = (window/2)*2 # force window to be a multiple of 2
    vec = numpy.zeros(len(x))
    prev = numpy.sum(x[0:window/2])
    fwindow = float(window)
    for i in range(window/2):
        vec[i] = prev/float(i + window/2)
        prev += x[i+window/2]
    for i in range(window/2, len(x)-window/2):
        vec[i] = prev/fwindow
        prev += x[i+window/2] - x[i-window/2]
    for i in range(len(x)-window/2, len(x)):
        vec[i] = prev/(float(len(x) - i))
        prev -= x[i-window/2]
    print "Smooth of len " + str(window) + " in " + str(time.time() - t0)
    return vec

def GCvsDepthPlot(GCcont, depth):
    yval = numpy.zeros(100)
    ynorm = numpy.zeros(100)
    for i in range(len(GCcont)):
        yval[math.floor(GCcont[i]*100.0)] += depth[i]
        ynorm[math.floor(GCcont[i]*100.0)] += 1.0
    for i in range(100):
        if ynorm[i] == 0.0:
            ynorm[i] += 1.0
    yval = yval/ynorm
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(yval, "*")
    #plt.xlabel("GC %")
    #plt.ylabel("Average depth at given percentage")
    #plt.title("Average Depth vs GC %")

    depthMax = numpy.max(yval)

    GCcorrection = numpy.ones(len(depth))

    for i in range(len(GCcont)):
        GCcorrection[i] = depthMax/yval[math.floor(GCcont[i]*100.0)]

    return GCcorrection

def GCcontentFunc(infile, window):
    handle = open(infile)
    parsed = SeqIO.parse(handle, "fasta")
    pos = 0
    GCcontent = []
    seqBuff = []
    GCtot = 0
    for seq_record in parsed:
        print len(seq_record.seq)
        for base in seq_record.seq:
            pos += 1
            if pos < window:
                seqBuff.append(base)
                if base == 'G' or base == 'g' or base == 'C' or base == 'c':
                    GCtot += 1
            else:
                GCcontent.append(float(GCtot)/float(window))
                seqBuff.append(base)
                if base == 'G' or base == 'g' or base == 'C' or base == 'c':
                    GCtot += 1
                popped = seqBuff.pop(0)
                if popped == 'G' or popped == 'g' or popped == 'C' or popped == 'c':
                    GCtot -= 1
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(smooth(numpy.asarray(GCcontent), 10000))
    #ax.set_title('Average GC content (window 10kbp)')
    #ax.set_xlabel('Position')
    #ax.set_ylabel('Avg GC content %')
    #plt.show()
    return GCcontent

if __name__ == '__main__':
    if len(sys.argv) not in (2,4):
        print __usage__
        sys.exit(0)
    c = ReadInInfo(sys.argv[1])
    if len(sys.argv) == 2:
        ScoreContigs(c)
        ALEcontigs, lens = c[0].thresholdMe()
        print len(ALEcontigs)
        lens.sort()
        lens.reverse()
        print lens
        print len(ALEcontigs)
        summer = 0
        for i in range(len(lens)):
            summer += lens[i]
            if summer > float(c[0].length)/2.0:
                print lens[i]
                break

    elif len(sys.argv) == 4:
        #ScoreContigs(c)
        ScoreContigsSection(c, int(sys.argv[2]), int(sys.argv[3]))
