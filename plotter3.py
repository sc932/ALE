#!/usr/bin/python

# (C) 2011 Scott Clark

"""ALE_plotter - a plotting package for ALE scoring output

    __copyright__

    __full_usage__

    Tested on 64-bit linux (ubuntu)

    Depends upon:
        matplotlib
            http://matplotlib.sourceforge.net/
        numpy
            http://numpy.scipy.org/
        mpmath
            http://code.google.com/p/mpmath/
        python
            2.6+
            http://www.python.org/

    Author:
        __author__

    Version:
        __version__ of __version_date__
"""

__version__ = "0.1"
__version_date__ = "28 September 2011"
__usage__ = """Usage: ./ALE_plotter.py [-options] <inputfile>

where basic options are:
-h : show brief help on version and full usage
"""
__author__ = "Scott Clark <sc932 at cornell dot edu>"
__copyright__ = """
                                   ALE

                             COPYRIGHT NOTICE
                               (insert GPL)

Scott Clark Copyright 2011
"""

import matplotlib.pylab as plt
import numpy
import sys
import mpmath
import logging

class Contig():
    """A contig from an ALE assembly

    more info

    Attributes:
        name: The name of the contig (defaults to unnamed)
        length: The length of the contig
        depth: A numpy vector of depths for each position
        depth_prob: A numpy vector of depth probabilities for each position
        placement_prob: A numpy vector of placement probabilities for each position
        kmer_prob: A numpy vector of kmer probabilities for each position
        total_prob: A numpy vector of total probabilities for each position
        total_prob = depth_prob + placement_prob + kmer_prob
    """
    def __init__(self, length=0, name = "unnamed"):
        """Inits Contig with name and length (and prob vectors)
        
        Kwargs:
            length: The length of the contig (>=0)
            
        Raises:
            ValueError: length must be >= 0
        """

        if length < 0:
          raise ValueError("length must be >= 0")

        self.name = name
        self.length = length
        self.depth = numpy.ones(length)
        self.depth_prob = numpy.zeros(length)
        self.placement_prob = numpy.zeros(length)
        self.kmer_prob = numpy.zeros(length)
        self.total_prob = numpy.zeros(length)

    def plot(self, start = 0, end = 0, plot_type = "tdpk", save_figure = False, depth_smoothing_width = 10000, placement_smoothing_width = 1000, kmer_smoothing_width = 1000, thresh = 0.99999):
        """Plots the contig

        Kwargs:
            start: The start of the plot (position) (>0, <end)
            end: Then end of the plot (position) (<=length)
            plot_type: Type (d)epth, (t)otal, (p)lacement, (k)mer in form "dpkt" or similar
            save_figure: Whether to save the figure as a .png image
            depth_smoothing_width: Width of window for averaging of depth scores
            placement_smoothing_width: Width of window for averaging of depth scores
            kmer_smoothing_width: Width of window for averaging of depth scores
            thresh: Threshold for error line using assumed normal distribution of data

        Raises:
            ValueError: plot_type must be some combination of 't','d','p','k', found: %s.
            ValueError: start must be less than end and greater than 0.
            ValueError: end must be less than length of contig.
            ValueError: placement_smoothing_width must be >= 0
            ValueError: kmer_smoothing_width must be >= 0
            ValueError: depth_smoothing_width must be >= 0
        """
        # sanitize input

        # check plot_type
        for letter in plot_type:
            if letter not in "tdpk":
                raise ValueError(("plot_type must be some combination of 't','d','p','k', found: %s." % letter))

        # check start, end
        if end < start or start < 0:
            raise ValueError("start must be less than end and greater than 0.")
          
        if end > self.length:
            raise ValueError("end must be less than length of contig.")

        # check smoothing
        if placement_smoothing_width < 0:
            raise ValueError("placement_smoothing_width must be >= 0")

        if kmer_smoothing_width < 0:
            raise ValueError("kmer_smoothing_width must be >= 0")

        if depth_smoothing_width < 0:
            raise ValueError("depth_smoothing_width must be >= 0")

        # sub methods

        def set_labels(ax, current_subplot, colors, plot_type, data_dict, sigma, number=3):
            """Sets labels on y-axis for each subplot"""
            assert(current_subplot == len(colors))
            ticks = []
            labels = []
            i = -1
            for typer in plot_type:
                i += 1
                for j in range(-number, number + 1):
                    ticks.append(4 + 7*i + j)
                    if j < 0:
                        #labels.append(str(j) + '$\sigma$ = ' + str(data_dict[typer] - j*sigma)[0:5])
                        labels.append(str(data_dict[typer] - j*sigma)[0:5])
                    else:
                        #labels.append('+' + str(j) + '$\sigma$ = ' + str(data_dict[typer] + j*sigma)[0:5])
                        labels.append(str(data_dict[typer] + j*sigma)[0:5])

            i = 0
            special_labels = {'t':'Total', 'd':'Depth', 'p':'Place', 'k':'K-mer'}
            for typer in plot_type:
                labels[3 + 7*i] = special_labels[typer] + ' ' + str(data_dict[typer])[0:5]
                i += 1

            ax.set_yticks(ticks)
            ax.set_yticklabels(labels)

            labels = ax.get_yticklabels()
            for i in range(current_subplot):
                for j in range(2*number + 1):
                    labels[i*number*2 + i + j].set_color(colors[i])

        def plot_std_marks(current_subplot, ax, length, color, number=3):
            """Plots std marks every sigma from the mean for each subplot"""
            for i in range(-number, number + 1):
                ax.plot([0, length], [4 + 7*current_subplot + i, 4 + 7*current_subplot + i], color + '--', alpha = 0.25)

        def format_data_for_plot(current_subplot, sigma, data):
            """Formats data to be plotted on the predefined grid"""
            # scale
            data = 1.0/sigma*data

            mu = numpy.mean(data)

            # translate
            # want to put a line every 7 (3 sig buffer), with a 1 sig buffer on top/bottom
            data = (4 + 7*current_subplot) - mu + data

            return data

        def find_alpha(data, threshold):
            """Finds the alpha threshold by fitting data to a normal distribution

            Args:
                data: The data to be fit to a normal
                threshold: The value we solve the inverse empirical cdf for

            Returns:
                The x value that corresponds to threshold in
                the empirical normal distribution defined by the data
            """
            mu = numpy.mean(data)
            std = numpy.std(data)
            # inv normal cdf
            return numpy.sqrt(2*std**2) * mpmath.erfinv(2*threshold - 1) + mu

        def find_threshold(data, method="min", thresh=.99, bins=1000, plot_figure=False):
            """Finds the two thresholds for errors given the data

            Args:
                data: The data to fit

            Kwargs:
                method: Whether to us [min,median,mean] of data in each bin
                thresh: Threshold for find_alpha
                bins: Number of pieces of the data we look at
                plot: Whether to plot the cdf and the two alpha cutoffs

            Returns:
                A soft threshold (alpha0) and A strong threshold (alpha1)

            Raises:
                ValueError: method needs to be in [\"min\",\"median\",\"mean\"], found %s
            """

            if method not in ["min","median","mean"]:
                raise ValueError("method needs to be in [\"min\",\"median\",\"mean\"], found %s" % method)

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
            
            alpha0 = -find_alpha(aSet, thresh)
            alpha1 = -find_alpha(aSet[50:-50], thresh)
            print alpha0, alpha1

            if plot_figure:
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

        # Main plotting code

        # correct for edge effects
        largest_smooth = max(depth_smoothing_width, placement_smoothing_width)
        largest_smooth = max(largest_smooth, kmer_smoothing_width)

        if start == 0:
            start += largest_smooth
        if end == 0:
            end = self.length
            end -= largest_smooth

        # make a new figure
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # smooth them out!
        starting_smooth_point = max(0,start-largest_smooth)
        ending_smooth_point = min(self.length-1,end+largest_smooth)

        depth_prob = smooth(self.depth_prob[starting_smooth_point:ending_smooth_point], depth_smoothing_width)[start:end]

        placement_prob = smooth(self.placement_prob[starting_smooth_point:ending_smooth_point], placement_smoothing_width)[start:end]

        kmer_prob = smooth(self.kmer_prob[starting_smooth_point:ending_smooth_point], kmer_smoothing_width)[start:end]

        # find totals
        total_prob = depth_prob + placement_prob + kmer_prob
        totalSigma = numpy.std(total_prob)

        current_subplot = 0
        colorDict = {'t':'m', 'd':'r', 'p':'b', 'k':'g'}
        data_dict = {'t':total_prob, 'd':depth_prob, 'p':placement_prob, 'k':kmer_prob}
        meanDict = {'t':numpy.mean(total_prob), 'd':numpy.mean(depth_prob), 'p':numpy.mean(placement_prob), 'k':numpy.mean(kmer_prob)}
        colors = []

        alpha = find_threshold(total_prob, thresh = thresh, method = "median", plot_figure = True)
        ax.plot([0, end - start], [4 + alpha/totalSigma,4 + alpha/totalSigma], 'black')

        for typer in plot_type:
            color = colorDict[typer]
            ax.plot(format_data_for_plot(current_subplot, totalSigma, data_dict[typer]), color)
            plot_std_marks(current_subplot, ax, end - start, color)
            colors.append(color)
            current_subplot += 1

        # fix labels
        set_labels(ax, current_subplot, colors, plot_type, meanDict, totalSigma)

        ax.set_title('Average Likelihoods')
        ax.set_xlabel('Position')
        ax.set_ylabel('Avg (log) Likelihood')
        ax.set_xlim((0, end-start))
        ax.set_ylim((0,current_subplot*7))

        plt.show()
        if save_figure:
            plt.savefig(self.name + '.png')



def read_in_info(placement_file):
    """Reads in an ALE placement file, returns a list of Contigs.

    Args:
        placement_file: An ALE placement file (*.ale)
            must be in the following format:

            # Reference: gi|170079663|ref|NC_010473.1| 350000
            # contig position depth ln(depthLike) ln(placeLike) ln(kmerLike) ln(totalLike)
            0 0 1.000000 -60.000000 0.194888 -5.760798 -65.565910
            0 1 3.000000 -60.000000 0.466271 -5.608334 -65.142063
            0 2 5.000000 -60.000000 0.010585 -5.541655 -65.531071
            0 3 12.000000 -60.000000 -0.057731 -5.380759 -65.438491

            Specific lines:
            0. The length of the contig is int(line[0].split(' ')[3]) == 350000
               name CANNOT be "position"
               The name of the contig is line.split(' ')[2]
            1. This line is ignored and lists what is in the columns of following lines
            2+. See 1

    Returns:
        A list of Contigs (see class Contig)

    Raises:
        IOError: An error occured accessing the placement file.
        FormattingError: The placement file was not formatted correctly.
    """

    MINIMUM_VALUE = -60.0 # minimum value we allow a position to have

    ale_placement_file = open(placement_file, 'r')

    contigs = []

    for line in ale_placement_file:
        if line[0] == '#':
            tName = line.split(' ')[2]
            if tName != "position":
                tLen = int(line.split(' ')[3])
                contigs.append(Contig(tLen, name = tName))
                place = 0
                print "Reading in contig: " + tName + " len " + str(tLen)
                print ""
                bar = progressBar(0, tLen, 42)
        else:
            data = line.split(' ')
            contigs[-1].depth[place] = numpy.double(data[2])
            for i in range(1,5):
                if "-nan"==data[i] or "nan"==data[i] or "inf"==data[i] or "-inf"==data[i] or numpy.double(data[i]) != numpy.double(data[i]):
                    data[i] = MINIMUM_VALUE # Predefined threshold
            contigs[-1].depth_prob[place] = numpy.double(data[3])
            contigs[-1].placement_prob[place] = numpy.double(data[4])
            contigs[-1].kmer_prob[place] = numpy.double(data[5])
            contigs[-1].total_prob[place] = numpy.double(data[6])
            place += 1
            if (place)%(int(tLen)/40) == 0:
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

class Error(Exception):
    pass

def main():
    if len(sys.argv) < 2:
        print __usage__
        sys.exit(0)
        
    if sys.argv[1] == '--help' or sys.argv[1] == '-h' or sys.argv[1] == '-help' or sys.argv[1] == '--h':
        print __fullUsage__
        sys.exit(0)
    print "Executed Sucessfully"

if __name__ == '__main__':
    main()

