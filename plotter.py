import matplotlib.pylab as plt
import numpy
from Bio.Seq import Seq
from Bio import SeqIO
import scipy
import time, sys, math

__usage__ = "(C) Scott Clark 2010.\nUsage:\npython dataPlotter.py [paramFile] [placementFile]\nWhere these are files inputed and outputed by ALE\n"

class Library():
    def __init__(self):
        self.inwardInsertLen = {0:0}
        self.outwardInsertLen = {0:0}
        self.placementsPerRead = {0:0,1:0}
    def printMe(self, start = 0, end = 0):
        maxValIn = max(libraries[i].inwardInsertLen.keys())
        maxValOut = max(libraries[i].outwardInsertLen.keys())
        if end == 0:
            end = int(numpy.max(maxValIn, maxValOut))
        distIn = numpy.zeros(maxValIn+1)
        distOut = numpy.zeros(maxValOut+1)
        for j in libraries[i].inwardInsertLen.keys():
            distIn[j] = libraries[i].inwardInsertLen[j]
        for j in libraries[i].outwardInsertLen.keys():
            distOut[j] = libraries[i].outwardInsertLen[j]
        fig = plt.figure()
        plt.plot(distIn[start:end])
        plt.plot(distOut[start:end])

class Contig():
    def __init__(self, length):
        self.inDepth = numpy.ones(length)
        self.outDepth = numpy.ones(length)
        self.poissonProb = numpy.zeros(length)
        self.placementProb = numpy.zeros(length)
        self.placementNormalizer = numpy.zeros(length)
        self.kmerProb = numpy.zeros(length)
        self.kmerVec = numpy.zeros(257)
        self.seq = []
    def printMe(self, start = 0, end = 0):
        if end == 0:
            end = len(self.seq)
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax1.plot(self.inDepth[start:end], 'g')
        ax1.plot(self.outDepth[start:end], 'b')
        ax1.legend(['Inward Depth', 'Outward Depth'], 'upper right')
        ax1.set_title('Contig Depths')
        #ax1.set_xlabel('Position')
        ax1.set_ylabel('Depth')
        newXtick = []
        partitions = len(ax1.get_xticklabels())
        for i in range(partitions):
            newXtick.append(start + i*(end-start)/(partitions-1))
        ax1.set_xticklabels(newXtick)
        ax2 = fig.add_subplot(212)
        ax2.plot(self.poissonProb[start:end], 'r')
        ax2.plot(numpy.log(numpy.exp(-20) + self.placementProb/self.placementNormalizer)[start:end], 'b')
        ax2.plot(numpy.log(numpy.exp(-20) + self.kmerProb[start:end]), 'g')
        #ax2.legend(['Average Depth Likelihood', 'Average Placement Likelihood', 'Average k-mer Likelihood'], 'lower left')
        ax2.set_title('Average Likelihoods')
        ax2.set_xlabel('Position')
        ax2.set_ylabel('Avg (log) Likelihood')
        newXtick = []
        partitions = len(ax2.get_xticklabels())
        for i in range(partitions):
            newXtick.append(start + i*(end-start)/(partitions-1))
        ax2.set_xticklabels(newXtick)
        #plt.plot(numpy.log(contigs[i].placementProb/(contigs[i].inDepth + contigs[i].outDepth)))
        plt.savefig('ntest.png')
    def scoreMe(self, start = 0, end = 0):
        if end == 0:
            end = len(self.seq)
        scoreVec = numpy.zeros(end-start)
        print "Scoring contig..."
        bar = progressBar(0, end-start, 42)
        for i in range(start,end):
            if (i - start)%(end/40) == 0:
                bar(i - start)
            if not math.isnan(self.poissonProb[i]):
                scoreVec[i-start] += self.poissonProb[i]
            else:
                scoreVec[i-start] += -20.0
            if not math.isnan(numpy.log(self.placementProb[i]/self.placementNormalizer[i])):
                scoreVec[i-start] += numpy.log(self.placementProb[i]/self.placementNormalizer[i])
            else:
                scoreVec[i-start] += -20.0
            if not math.isnan(numpy.log(self.kmerProb[i])):
                scoreVec[i-start] += numpy.log(self.kmerProb[i])
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
        c[i].printMe(start=start, end=end)
        
def ScoreContigs(c):
    print "Scoring all contigs:"
    for i in range(len(c)):
        print "Contig " + str(i+1) + " score:" + str(c[i].scoreMe())

def ReadInInfo(paramFile, placementFile):
    fPa = open(paramFile, 'r')
    fPl = open(placementFile, 'r')
    fPl.readline() # explaination line
    libraries = []
    contigs = []
    contigFile = fPa.readline()[1:-1]
    numReadFiles = int(fPa.readline())
    inserts = []
    for i in range(numReadFiles):
        fPa.readline()
        infos = fPa.readline().split(',')
        inserts.append(float(infos[2]))
        inserts.append(float(infos[5]))
    maxInsert = max(numpy.asarray(inserts))
    
    print "Reading in contigs..."
    
    handle = open(contigFile)
    parsed = SeqIO.parse(handle, "fasta")
    for seq_record in parsed:
        contigs.append(Contig(len(seq_record.seq)))
        for i in range(len(seq_record.seq)):
            contigs[-1].seq.append(seq_record.seq[i])
            
    print "Done reading in contigs."
    print "Reading in libraries..."
    
    # read in the data
    it = 0
    normalizer = 0.0
    infoer = []
    for line in fPl:
        it+= 1
        if it%10000 == 0:
            print "Read in " + str(it) + " lines..."
        
        if line[0] == 'L':
            print "Reading in library " + str(len(libraries) + 1) + " ..."
            libraries.append(Library())
        else:
            infos = line.split(',') # Read Number, Placement number, total number of placements, likelihood, area1 start, area1 end, area 2 start, area 2 end, map type (0,1,4,5,8,9,12,13), assembly part
            for i in [0,1,2,4,5,6,7,8,9]:
                infos[i] = int(infos[i])
            infos[3] = float(infos[3])
            if infos[2] == 1: # a single placement
                libraries[-1].placementsPerRead[1] += 1
                # get the contig depths and insert dist
                if infos[8] == 13 or infos[8] == 9 or infos[8] == 12 or infos[8] == 8: # inward
                    contigs[infos[9]].inDepth[infos[4]:infos[5]] += 1
                    contigs[infos[9]].inDepth[infos[6]:infos[7]] += 1
                    if infos[7] - infos[4] not in libraries[-1].inwardInsertLen.keys():
                        libraries[-1].inwardInsertLen[infos[7] - infos[4]] = 1
                    else:
                        libraries[-1].inwardInsertLen[infos[7] - infos[4]] += 1
                else:
                    contigs[infos[9]].outDepth[infos[4]:infos[5]] += 1
                    contigs[infos[9]].outDepth[infos[6]:infos[7]] += 1
                    if infos[7] - infos[4] not in libraries[-1].outwardInsertLen.keys():
                        libraries[-1].outwardInsertLen[infos[7] - infos[4]] = 1
                    else:
                        libraries[-1].outwardInsertLen[infos[7] - infos[4]] += 1
                contigs[infos[9]].placementProb[infos[4]:infos[5]] += infos[3]
                contigs[infos[9]].placementProb[infos[6]:infos[7]] += infos[3]
                contigs[infos[9]].placementNormalizer[infos[4]:infos[5]] += 1.0
                contigs[infos[9]].placementNormalizer[infos[6]:infos[7]] += 1.0
            else:
                normalizer += infos[3]
                infoer.append(infos)
                if infos[2] - 1 == infos[1] and normalizer > 0.0:
                    for infos in infoer:
                        if infos[8] == 13 or infos[8] == 9 or infos[8] == 12 or infos[8] == 8: # inward
                            contigs[infos[9]].inDepth[infos[4]:infos[5]] += infos[3]/normalizer
                            contigs[infos[9]].inDepth[infos[6]:infos[7]] += infos[3]/normalizer
                            if infos[7] - infos[4] not in libraries[-1].inwardInsertLen.keys():
                                libraries[-1].inwardInsertLen[infos[7] - infos[4]] = infos[3]/normalizer
                            else:
                                libraries[-1].inwardInsertLen[infos[7] - infos[4]] += infos[3]/normalizer
                        else:
                            contigs[infos[9]].outDepth[infos[4]:infos[5]] += infos[3]/normalizer
                            contigs[infos[9]].outDepth[infos[6]:infos[7]] += infos[3]/normalizer
                            if infos[7] - infos[4] not in libraries[-1].outwardInsertLen.keys():
                                libraries[-1].outwardInsertLen[infos[7] - infos[4]] = infos[3]/normalizer
                            else:
                                libraries[-1].outwardInsertLen[infos[7] - infos[4]] += infos[3]/normalizer
                        contigs[infos[9]].placementProb[infos[4]:infos[5]] += infos[3]*infos[3]/normalizer
                        contigs[infos[9]].placementProb[infos[6]:infos[7]] += infos[3]*infos[3]/normalizer
                        contigs[infos[9]].placementNormalizer[infos[4]:infos[5]] += infos[3]/normalizer
                        contigs[infos[9]].placementNormalizer[infos[6]:infos[7]] += infos[3]/normalizer
                    normalizer = 0.0
                    infoer = []
                
    print "Done reading in libraries."
    
    print "Making Depth Statistics..."
    
    it = 0
    for c in contigs:
        it += 1
        depthNorm = numpy.sum(c.inDepth[maxInsert:len(c.seq)-maxInsert] + c.outDepth[maxInsert:len(c.seq)-maxInsert])/float(len(c.seq)-2*maxInsert)
        #print depthNorm
        print "Contig " + str(it) + " with E[depth] = " + str(depthNorm)+ ":"
        theLength = len(c.seq)
        bar = progressBar(0, theLength-1, 42)
        for i in range(theLength):
            if i%(theLength/40) == 0:
                bar(i)
            #if i > maxInsert and i < len(c.seq)-maxInsert:
            c.poissonProb[i] = poissonPMF(1 + c.inDepth[i] + c.outDepth[i], depthNorm)
        bar(theLength-1)
    
    print "\nMaking k-mer Statistics..."
    
    it = 0
    for c in contigs:
        it += 1
        print "Building vector for contig " + str(it)
        theLength = len(c.seq)
        bar = progressBar(0, theLength-1, 42)
        for i in range(theLength - 4):
            if i%(theLength/40) == 0:
                bar(i)
            c.kmerVec[kmerTranslate(c.seq[i:i+4])] += 1
        bar(theLength-1)
        print "\nCalculating statistics for contig " + str(it)
        # start up
        if kmerTranslate(c.seq[0:4]) != 256:
            c.kmerProb[0] += c.kmerVec[kmerTranslate(c.seq[0:4])]/float(theLength - 4)
            c.kmerProb[1] += 0.5*c.kmerVec[kmerTranslate(c.seq[0:4])]/float(theLength - 4)
            c.kmerProb[2] += (1/3.0)*c.kmerVec[kmerTranslate(c.seq[0:4])]/float(theLength - 4)
        if kmerTranslate(c.seq[1:5]) != 256:
            c.kmerProb[1] += 0.5*c.kmerVec[kmerTranslate(c.seq[1:5])]/float(theLength - 4)
            c.kmerProb[2] += (1/3.0)*c.kmerVec[kmerTranslate(c.seq[1:5])]/float(theLength - 4)
        if kmerTranslate(c.seq[2:6]) != 256:
            c.kmerProb[2] += (1/3.0)*c.kmerVec[kmerTranslate(c.seq[2:6])]/float(theLength - 4)
        bar = progressBar(0, theLength-1, 42)
        for i in range(3,theLength-4):
            if i%(theLength/40) == 0:
                bar(i)
            for j in range(4):
                if kmerTranslate(c.seq[i-j:i-j+4]) != 256:
                    c.kmerProb[i] += 0.25*c.kmerVec[kmerTranslate(c.seq[i-j:i-j+4])]/float(theLength - 4)
        bar(theLength-1)
    
    print "\nDone making k-mer Statistics."
                
                
    return contigs, libraries
    
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
    
def kmerTranslate(ins):
    # A = 00
    # T = 10
    # C = 01
    # G = 11
    num = 0
    for i in range(4):
        if ins[i] == 'T':
            num += pow(2,2*i)
        elif ins[i] == 'C':
            num += pow(2,2*i+1)
        elif ins[i] == 'G':
            num += pow(2,2*i) + pow(2,2*i+1)
        elif ins[i] == 'N':
            return 256
    return num
    
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
    if len(sys.argv) < 3:
        print __usage__
        sys.exit(0)
    c, l = ReadInInfo(sys.argv[1], sys.argv[2])
    if len(sys.argv) == 3:
        ScoreContigs(c)
    else:
        ScoreContigs(c)
        ScoreContigsSection(c, int(sys.argv[3]), int(sys.argv[4]))
