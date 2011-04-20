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
    def thresholdMe(self, threshold = -60, minLen = 1000):
        ALEcontigs = []
        lens = []
        passed = self.totalProb > threshold
        print passed[20:50]
        start = 0
        for i in range(len(passed)):
            if passed[i] == True and i > 0 and passed[i-1] == False:
                start = i
            elif passed[i] == False and i > 0 and passed[i-1] == True:
                end = i
                if end - start >= minLen:
                    ALEcontigs.append([start,end])
                    lens.append(end-start)
                start = end
        return ALEcontigs, lens
    def plotMe(self, start = 0, end = 0):
        if end < 1 or end < start or end > self.length:
            end = self.length
        fig = plt.figure()
        ax2 = fig.add_subplot(111)
        ax2.plot(self.depthProb[start:end], 'r')
        ax2.plot(self.placementProb[start:end], 'b')
        ax2.plot(self.kmerProb[start:end], 'g')
        ax2.plot(self.totalProb[start:end], 'm')
        ax2.legend(['Average Depth Likelihood', 'Average Placement Likelihood', 'Average k-mer Likelihood', 'total'], 'lower left')
        ax2.set_title('Average Likelihoods')
        ax2.set_xlabel('Position')
        ax2.set_ylabel('Avg (log) Likelihood')
        newXtick = []
        partitions = len(ax2.get_xticklabels())
        for i in range(partitions):
            newXtick.append(start + i*(end-start)/(partitions-1))
        ax2.set_xticklabels(newXtick)
        #plt.plot(numpy.log(contigs[i].placementProb/(contigs[i].inDepth + contigs[i].outDepth)))
        plt.show()
        plt.savefig(self.name + '.png')
    def scoreMe(self, start = 0, end = 0):
        if end < 1 or end < start or end > self.length:
            end = self.length
        scoreVec = numpy.zeros(end-start)
        print "Scoring contig " + self.name + "..."
        bar = progressBar(0, end-start, 42)
        for i in range(start,end):
            if (i - start)%(end/40) == 0:
                bar(i - start)
            if not math.isnan(self.depthProb[i]):
                scoreVec[i-start] += self.depthProb[i]
            else:
                scoreVec[i-start] += -20.0
            if not math.isnan(self.placementProb[i]):
                scoreVec[i-start] += self.placementProb[i]
            else:
                scoreVec[i-start] += -20.0
            if not math.isnan(self.kmerProb[i]):
                scoreVec[i-start] += self.kmerProb[i]
            else:
                scoreVec[i-start] += -20.0
            if math.isnan(scoreVec[i-start]) or scoreVec[i-start] < -60:
                scoreVec[i-start] = -60
        score = (1.0/(float(end - start)))*numpy.sum(scoreVec)
        bar(end-start)
        print ""
        return score
        
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
        if line[0] == '>':
            print ""
            tName = line.split('>')[1].split(' ')[0]
            tLen = int(line.split(' ')[1])
            contigs.append(Contig(tLen, name = tName))
            place = 0
            print "Reading in contig: " + tName
            bar = progressBar(0, tLen, 42)
        else:
            data = line.split(' ')
            contigs[-1].depth[place] = numpy.double(data[0])
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
        
if __name__ == '__main__':
    if len(sys.argv) not in (2,4):
        print __usage__
        sys.exit(0)
    c = ReadInInfo(sys.argv[1])
    if len(sys.argv) == 2:
        # ScoreContigs(c)
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
