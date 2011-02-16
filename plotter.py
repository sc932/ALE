import matplotlib.pylab as plt
import numpy
from Bio.Seq import Seq
from Bio import SeqIO
import scipy
import time, sys, math

__usage__ = "(C) Scott Clark 2010,2011.\nUsage:\npython dataPlotter.py [paramFile] [placementFile]\nWhere these are files inputed and outputed by ALE\n"

class Contig():
    def __init__(self, length):
        self.depth = numpy.ones(length)
        self.depthProb = numpy.zeros(length)
        self.placementProb = numpy.zeros(length)
        self.kmerProb = numpy.zeros(length)
    def printMe(self, start = 0, end = 0):
        if end == 0:
            end = len(self.seq)
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax1.plot(self.depth[start:end], 'g')
        ax1.set_title('Contig Depths')
        #ax1.set_xlabel('Position')
        ax1.set_ylabel('Depth')
        newXtick = []
        partitions = len(ax1.get_xticklabels())
        for i in range(partitions):
            newXtick.append(start + i*(end-start)/(partitions-1))
        ax1.set_xticklabels(newXtick)
        ax2 = fig.add_subplot(212)
        ax2.plot(self.depthProb[start:end], 'r')
        ax2.plot(numpy.log(numpy.exp(-30) + self.placementProb)[start:end], 'b')
        ax2.plot(numpy.log(numpy.exp(-30) + self.kmerProb[start:end]), 'g')
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
        plt.savefig('plotter_output.png')
    def scoreMe(self, start = 0, end = 0):
        if end == 0:
            end = len(self.seq)
        scoreVec = numpy.zeros(end-start)
        print "Scoring contig..."
        bar = progressBar(0, end-start, 42)
        for i in range(start,end):
            if (i - start)%(end/40) == 0:
                bar(i - start)
            if not math.isnan(self.depthProb[i]):
                scoreVec[i-start] += self.depthProb[i]
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
        print "score", score
        return score