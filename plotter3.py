#!/usr/bin/python

# (C) 2011 Scott Clark

import matplotlib.pylab as plt
import numpy
from Bio.Seq import Seq
from Bio import SeqIO
import scipy
import time, sys, math
import mpmath

class Contig():
  def __init__(self, length, name = "unnamed", probNorm = -10.0):
    self.name = name
    self.length = length
    self.depth = numpy.ones(length)
    self.depthProb = numpy.zeros(length)
    self.placementProb = numpy.zeros(length)
    self.kmerProb = numpy.zeros(length)
    self.totalProb = numpy.zeros(length)
    self.avgProbNorm = probNorm
    self.depthProbNormalizer = 1.0
    self.placementProbNormalizer = 1.0
    self.kmerProbNormalizer = 1.0

  def plot(self, start = 0, end = 0, plotType = "tdpk", saver = False, smoothLevelD = 10000, smoothLevelP = 1000, smoothLevelK = 1000):
        # correct for edge effects
        largestSmooth = max(smoothLevelD, smoothLevelP)
        largestSmooth = max(largestSmooth, smoothLevelK)

        if start == 0:
          start += largestSmooth
        if end < 1 or end < start or end > self.length:
            end = self.length
            end -= largestSmooth

        # make a new figure
        fig = plt.figure()
        ax2 = fig.add_subplot(111)
        legendList = []
        
        # smooth them out!
        depthProb = smooth(self.depthProb[max(0,start-largestSmooth):min(self.length-1,end+largestSmooth)], smoothLevelD)[start:end]
        placementProb = smooth(self.placementProb[max(0,start-largestSmooth):min(self.length-1,end+largestSmooth)], smoothLevelP)[start:end]
        kmerProb = smooth(self.kmerProb[max(0,start-largestSmooth):min(self.length-1,end+largestSmooth)], smoothLevelK)[start:end]

        # find totals
        totalProb = depthProb + placementProb + kmerProb
        totalSigma = numpy.std(totalProb)

        numOn = 0
        colorDict = {'t':'m', 'd':'r', 'p':'b', 'k':'g'}
        dataDict = {'t':totalProb, 'd':depthProb, 'p':placementProb, 'k':kmerProb}
        meanDict = {'t':numpy.mean(totalProb), 'd':numpy.mean(depthProb), 'p':numpy.mean(placementProb), 'k':numpy.mean(kmerProb)}
        colors = []

        alpha = findThreshold(totalProb, thresh = .99999, method = "median", plot = True)
        ax2.plot([0, end - start], [4 + alpha/totalSigma,4 + alpha/totalSigma], 'black')

        for typer in plotType:
          color = colorDict[typer]
          ax2.plot(getPlotFotmatted(numOn, totalSigma, dataDict[typer]), color)
          plotStdMarks(numOn, ax2, end - start, color)
          colors.append(color)
          numOn += 1

        # fix labels
        setLabels(ax2, numOn, colors, plotType, meanDict, totalSigma)

        ax2.set_title('Average Likelihoods')
        ax2.set_xlabel('Position')
        ax2.set_ylabel('Avg (log) Likelihood')
        ax2.set_xlim((0, end-start))
        ax2.set_ylim((0,numOn*7))

        plt.show()
        if saver == True:
            plt.savefig(self.name + '.png')

def setLabels(ax, numOn, colors, plotType, dataDict, sigma, number = 3):
  assert(numOn == len(colors))
  ticks = []
  labels = []
  i = -1
  for typer in plotType:
    i += 1
    for j in range(-number, number + 1):
      ticks.append(4 + 7*i + j)
      if j < 0:
        #labels.append(str(j) + '$\sigma$ = ' + str(dataDict[typer] - j*sigma)[0:5])
        labels.append(str(dataDict[typer] - j*sigma)[0:5])
      else:
        #labels.append('+' + str(j) + '$\sigma$ = ' + str(dataDict[typer] + j*sigma)[0:5])
        labels.append(str(dataDict[typer] + j*sigma)[0:5])

  i = 0
  specialLabels = {'t':'Total', 'd':'Depth', 'p':'Place', 'k':'K-mer'}
  for typer in plotType:
    labels[3 + 7*i] = specialLabels[typer] + ' ' + str(dataDict[typer])[0:5]
    i += 1

  ax.set_yticks(ticks)
  ax.set_yticklabels(labels)

  labels = ax.get_yticklabels()
  for i in range(numOn):
    for j in range(2*number + 1):
      labels[i*number*2 + i + j].set_color(colors[i])

def plotStdMarks(numOn, ax, length, color, number = 3):
  for i in range(-number, number + 1):
    ax.plot([0, length], [4 + 7*numOn + i, 4 + 7*numOn + i], color + '--', alpha = 0.25)

def getPlotFotmatted(numOn, sigma, data):
  # scale
  data = 1.0/sigma*data

  mu = numpy.mean(data)

  # translate
  # want to put a line every 7 (3 sig buffer), with a 1 sig buffer on top/bottom
  data = (4 + 7*numOn) - mu + data

  return data

def findAlpha(data, threshold):
  mu = numpy.mean(data)
  std = numpy.std(data)
  return numpy.sqrt(2*std**2) * mpmath.erfinv(2*threshold - 1) + mu

def findThreshold(data, method = "min", thresh = .99, bins = 1000, plot = False):
  aSet = numpy.zeros(bins)
  for i in range(bins):
    if method == "median":
      aSet[i] = numpy.median(data[len(data)*i/bins:(len(data)*(i+1))/bins])
    elif method == "mean":
      aSet[i] = numpy.mean(data[len(data)*i/bins:(len(data)*(i+1))/bins])
    elif method == "min":
      aSet[i] = numpy.min(data[len(data)*i/bins:(len(data)*(i+1))/bins])
  aSet.sort()
  #aSet = list(abs(aSet))
  aSet = list(aSet - numpy.mean(data))
  aSet.reverse()
  
  alpha0 = -findAlpha(aSet, thresh)
  alpha1 = -findAlpha(aSet[50:-50], thresh)
  print alpha0, alpha1

  if plot:
    ySet = numpy.zeros(bins)
    for i in range(bins):
      ySet[i] = float(i + 1)/float(bins)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(aSet, ySet)
    ax.plot([aSet[0], aSet[-1]], [thresh, thresh], 'r') # threshold
    ax.plot([0,0], [0,1], 'black') # 0-line
    ax.plot([alpha0,alpha0], [0,1], 'm') # a0-line
    ax.plot([alpha1,alpha1], [0,1], 'g') # a1-line
    ax.set_xlim((aSet[0], aSet[-1]))
    ax.set_ylim(0,1)

  return alpha1

def ReadInInfo(placementFile):
    fPl = open(placementFile, 'r')

    contigs = []

    for line in fPl:
        if line[0] == '#':
            print ""
            try:
              tName = line.split(' ')[2]
            except IndexError:
              tName = ""
            if tName != "position":
              tLen = int(line.split(' ')[3])
              contigs.append(Contig(tLen, name = tName))
              place = 0
              print "Reading in contig: " + tName + " len " + str(tLen)
              bar = progressBar(0, tLen, 42)
        else:
            data = line.split(' ')
            contigs[-1].depth[place] = numpy.double(data[2])
            for i in range(1,5):
                if "-nan"==data[i] or "nan"==data[i] or "inf"==data[i] or "-inf"==data[i] or numpy.double(data[i]) != numpy.double(data[i]):
                    data[i] = -60.0
            contigs[-1].depthProb[place] = numpy.double(data[3])
            contigs[-1].placementProb[place] = numpy.double(data[4])
            contigs[-1].kmerProb[place] = numpy.double(data[5])
            contigs[-1].totalProb[place] = numpy.double(data[6])
            place += 1
            if (place)%(tLen/40) == 0:
                bar(place)
    print "\nYou now have a list of contigs, try contig[0].plot()"
    return contigs

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


def main():
  print "Executed Sucessfully"

if __name__ == '__main__':
  main()

