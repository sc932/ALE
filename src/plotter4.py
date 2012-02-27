#!/usr/bin/python

# (C) 2011 Scott Clark

"""plotter3 - a plotting package for ALE scoring output

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
        pymix
            http://www.pymix.org/pymix/
"""

__version__ = "0.2"
__usage__ = """Usage: ./ALE_plotter.py [-options] <inputfile.ale>

where basic options are:
  -h      : show brief help on version and full usage
  -nosave : do not save the figure as a pdf (instead plot to screen)
"""
__full_usage__="""Usage: ./ALE_plotter.py [-options] <inputfile.ale>

where basic options are:
  -h      : show brief help on version and full usage
  -nosave : do not save the figure as a pdf (instead plot to screen)

parameter options accepting <f>loats and <i>ntegers and <s>trings (default):
  -s <i>   : the starting position to plot (for all contigs, ie a single insert length)
  -e <i>   : the ending position of the plot
  -rs      : recursive search (for errors)
  -pt <s>  : plot type 'i'nsert 'k'mer 'p'lacement 'd'epth (-pt dpkt)
  -dsw <i> : depth smoothing window, averaging over position (-dsw 10000)
  -psw <i> : placement smoothing window (-psw 1000)
  -ksw <i> : kmer smoothing window (-ksw 1000)
  -isw <i> : insert smoothing window (-ksw 1000)
  -t <f>   : threshold percentage, see paper (-t 0.99999)
  -pt <f>  : plot threshold, only plot if more than % of errors (-pt 0.0)
  -st <i>  : number of standard deviations to engage threshold (-st 5)
  -fn <s>  : figure name (default: contig name)
  -mps <i> : minimum plot size in bp (-mps 20000)
  -sc <s>  : plot only a specific contig (ie -sc contigName213)
  -pmo     : plot meta information only (off)
  -dpm     : don't plot meta information at all (off)
  -wo      : turn on score weighting (off)
  -pw <f>  : placement weighting (1.0)
  -kw <f>  : kmer weighting (1.0)
  -dw <f>  : depth weighting (1.0)
  -iw <f>  : insert weighting (1.0)
"""
__author__ = "Scott Clark <sc932 at cornell dot edu>"
__copyright__ = """
                                   ALE

                             COPYRIGHT NOTICE
                               (insert GPL)

Scott Clark Copyright 2011
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import sys
import mixture # http://www.pymix.org/pymix/
import Smooth
import ProgressBar

#import logging

class CommandLineParameters(object):
    """User params"""
    def __init__(self):
        # parameters
        self.parameter_list = []
        self.type_list = []
        self.lookup_by_name = {}
        self.lookup_by_command = {}

        # inputs
        self.input_list = []
        self.input_names = {}

    def read_sys_args(self):
        if len(sys.argv) < 2:
            print __usage__
            sys.exit(0)
            
        if sys.argv[1] == '--help' or sys.argv[1] == '-h' or sys.argv[1] == '-help' or sys.argv[1] == '--h':
            print __full_usage__
            sys.exit(0)

        if len(sys.argv) == 2:
            arg_on = 1
        else:
            # read in command line arguments
            arg_on = 1
            while(arg_on + 1 < len(sys.argv)):            
                if sys.argv[arg_on] in self.get_val_set():
                    # user defined value
                    self.set_by_command(sys.argv[arg_on], sys.argv[arg_on + 1])
                    arg_on += 2
                elif sys.argv[arg_on] in self.get_bool_set():
                    # flip the switch
                    self.set_by_command(sys.argv[arg_on], None)
                    arg_on += 1
                else:
                    print "Did not recognize command line argument %s." % sys.argv[arg_on]
                    print "Try -h for help."
                    exit(0)

        for i, ins in enumerate(sys.argv[arg_on:]):
            self.input_list[i] = ins        

    def add_input(self, name):
        number_on = len(self.input_list)
        self.input_list.append("")
        self.input_names[name] = number_on

    def get_input(self, name):
        if name not in self.input_names:
            raise ValueError("%s not in input_names" % name)
        return self.input_list[self.input_names[name]]

    def add_parameter(self, name=None, command=None, val_type=str, init_val=None):
        if not name:
            raise ValueError("add_parameter needs a name supplied")
        if not command:
            raise ValueError("add_parameter needs a command supplied")
        number_on = len(self.parameter_list)
        self.parameter_list.append(init_val)
        self.type_list.append(val_type)
        self.lookup_by_name[name] = number_on
        self.lookup_by_command[command] = number_on

    def get_val_set(self):
        val_set = []
        for command in self.lookup_by_command:
            if self.type_list[self.lookup_by_command[command]] in (str, int, float):
                val_set.append(command)
        return val_set

    def get_bool_set(self):
        bool_set = []
        for command in self.lookup_by_command:
            if self.type_list[command] not in (str, int, float):
                bool_set.append(command)
        return bool_set

    def set_by_lookup_val(self, lookup_val, value):
        if value != None:
            try:
                value = self.type_list[lookup_val](value) # try to cast to type
            except:
                raise ValueError("Error setting parameter of %s with %s", (str(self.type_list[lookup_val]), str(value)))
            self.parameter_list[lookup_val] = value
        else:
            self.parameter_list[lookup_val] = not self.parameter_list[lookup_val]

    def set_by_command(self, command=None, value=None):
        if command not in self.lookup_by_command:
            raise ValueError("%s command not in param table")
        else:
            lookup_val = self.lookup_by_command[command]

        self.set_by_lookup_val(lookup_val, value)

    def get_by_command(self, command=None):
        if command not in self.lookup_by_command:
            raise ValueError("%s command not in param table")
        return self.parameter_list[self.lookup_by_command[command]]

    def set(self, name=None, value=None):
        if name not in self.lookup_by_name:
            raise ValueError("%s name not in param table")
        else:
            lookup_val = self.lookup_by_name[name]

        self.set_by_lookup_val(lookup_val, value)

    def get(self, name=None):
        if name not in self.lookup_by_name:
            raise ValueError("%s name not in param table")
        return self.parameter_list[self.lookup_by_name[name]]

class ALEFigure(object):
    """an ALE figure"""

    def __init__(self, start, end):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax2 = self.ax.twinx()

        self.start = start
        self.end = end
        
        self.threshold_window_set = []

    def add_threshold_windows(self, windows):
        self.threshold_window_set.extend(windows)

    def threshold_windows_in_last_half(self):
        for thresh_window in self.threshold_window_set:
            if thresh_window.start < (self.end-self.start)/2:
                return False
        return True

    def threshold_windows_in_first_half(self):
        for thresh_window in self.threshold_window_set:
            if thresh_window.start > (self.end-self.start)/2:
                return False
        return True

    def plot_threshold_windows(self):        
        # plot the thresholds
        for thresh_window in self.threshold_window_set:
            print "thresh window [%d,%d] %s" % (thresh_window.start, thresh_window.end, thresh_window.type_of)
            self.ax.axvspan(thresh_window.start, thresh_window.end, facecolor='r', alpha=0.1)          

    def get_percent_thresholded(self):
        total_below_threshold = numpy.zeros(self.end - self.start)
        for thresh_window in self.threshold_window_set:
            total_below_threshold[thresh_window.start:thresh_window.end] += 1
        summer = 0
        for pos in total_below_threshold:
            if pos:
                summer += 1
        return float(summer)/float(self.end - self.start)

    def set_labels(self, ax, plot_type, prob_vecs, std_witdh=2, twin=False, std_width=2):
        """Sets labels on y-axis for each subplot"""
        number_marks=3
        ticks = []
        labels = []
        i = -1
        for typer in plot_type:
            i += 1
            for j in range(-number_marks, number_marks + 1):
                ticks.append(4 + 7*i + j)
                if not twin:
                    if j < 0:
                        labels.append(str(j*std_width) + '$\sigma$')
                        #labels.append(str(j) + '$\sigma$ = ' + str(data_dict[typer] - j*sigma)[0:5])
                        #labels.append(str(data_dict[typer] - j*sigma)[0:5])
                    else:
                        labels.append(' ') # leave positive labels blank
                        #labels.append('+' + str(j) + '$\sigma$')
                        #labels.append('+' + str(j) + '$\sigma$ = ' + str(data_dict[typer] + j*sigma)[0:5])
                        #labels.append(str(data_dict[typer] + j*sigma)[0:5])
                else:
                    # TODO clean up
                    if j == 0:
                        labels.append(str(prob_vecs[typer].thresh_main_mean)[:7])
                    elif j < 0:
                        labels.append('(' + str(-j*std_width*prob_vecs[typer].thresh_main_std)[:6] + ')')
                    else:
                        labels.append(' ')

        i = 0
        if not twin:
            special_labels = {'i':'Insert', 'd':'Depth', 'p':'Place', 'k':'K-mer'}
            for typer in plot_type:
                labels[3 + 7*i] = special_labels[typer] + ' ' + '$\mu$' #str(data_dict[typer])[0:5]
                i += 1

        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)

        labels = ax.get_yticklabels()
        for i, typer in enumerate(list(plot_type)):
            for j in range(2*number_marks + 1):
                labels[i*number_marks*2 + i + j].set_color(prob_vecs[typer].color)

    def plot_std_marks(self, current_subplot, length, color, number=3):
        """Plots std marks every sigma from the mean for each subplot"""
        for i in range(-number, number + 1):
            self.ax.plot([0, length], [4 + 7*current_subplot + i, 4 + 7*current_subplot + i], color + '--', alpha = 0.1)

    def format_data_for_plot(self, current_subplot, mean, sigma, data, start, end, std_width=2.0):
        """Formats data to be plotted on the predefined grid"""
        # scale
        data = 1.0/(sigma*std_width)*data[start:end]

        #mu = numpy.mean(data)

        # translate
        # want to put a line every 7 (3 sig buffer), with a 1 sig buffer on top/bottom
        data = (4 + 7*current_subplot) - 1.0/(sigma*std_width)*mean + data

        for point in data:
            if point < 0.0:
                point = 0.0

        return data

    def plot_cleanup(self, start, end, num_subplots):
        self.ax.set_ylabel('Avg (log) Likelihood')
        self.ax.set_xlim((0, end-start))
        self.ax.set_ylim((0,num_subplots*7))

class ThresholdViolation(object):
    """A threshold of an ALE score vector"""
    def __init__(self, type_of, start, end, score):
        self.type_of = type_of
        self.start = start
        self.end = end
        self.score = score

class LikelihoodVector(object):
    """A likelihood vector of some type in a contig"""
    def __init__(self, length=None, color=None):
        if not length or length < 0:
            raise ValueError("length must be >= 0")
        if not color:
            raise ValueError("need to specify color")
        self.length = length
        self.color = color
        self.prob = numpy.zeros(length)
        self.prob_smoothed = numpy.zeros(length)
        self.smoothing_width = 1

        self.thresh_main_mean = None
        self.thresh_main_std = None

    def smooth_prob(self, smoothing_width=None):
        """Smooth a specific prob vector"""
        if smoothing_width and smoothing_width < 0:
            raise ValueError("smoothing_width must be >= 0")

        if not smoothing_width:
            self.smoothing_width = int(numpy.max((1, self.length/500.0)))
        else:
            self.smoothing_width = smoothing_width

        self.prob_smoothed = Smooth.smooth(self.prob, self.smoothing_width)[self.smoothing_width:-self.smoothing_width]

    def get_threshold_windows(self, data, user_params, typer='?', thresh_mult=-5.0):
        """Returns the start and end points of the windows that cross the threshold with some constraints"""
        # TODO make more pythonic
        # get the starts
        total_sigma = self.thresh_main_std
        data_mean = self.thresh_main_mean
        # TODO base off of user_params
        threshold = thresh_mult*total_sigma
        cross_thresh = user_params.get("threshold_percent")
        len_thresh = user_params.get("threshold_width")

        starts = []
        ends = []
        started = False
        end_started = False
        end_point = len(data)
        start_point = 0
        for position, point in enumerate(data):
            # starts
            if point/total_sigma < threshold/total_sigma + data_mean/total_sigma:
                if not started:
                    started = True
                    start_point = position
                    window_len = 1
                    cross_total = 1
                else:
                    window_len += 1
                    cross_total += 1
            else:
                if not started:
                    pass
                else:
                    window_len += 1
                    if float(cross_total)/float(window_len) < cross_thresh:
                        started = False
                        if window_len - 1 > len_thresh:
                            starts.append(start_point)
            # ends
        if started:
            starts.append(start_point)

        for position, point in enumerate(data):
            reverse_position = len(data) - 1 - position
            if reverse_position in starts and end_started:
                end_started = False
                ends.append(end_point)
            else:
                if data[reverse_position]/total_sigma < threshold/total_sigma + data_mean/total_sigma:
                    if not end_started:
                        end_started = True
                        end_point = reverse_position
                        end_window_len = 1
                        end_cross_total = 1
                    else:
                        end_window_len += 1
                        end_cross_total += 1
                else:
                    if not end_started:
                        pass
                    else:
                        end_window_len += 1
                        if float(end_cross_total)/float(end_window_len) < cross_thresh:
                            end_started = False
                            if end_window_len - 1 > len_thresh:
                                ends.append(end_point)

        if end_started:
            ends.append(end_point)
        ends.reverse()

        print "starts and ends of windowing threshold"
        print starts, ends

        threshold_windows = []

        for i in range(len(starts)):
            score = numpy.sum(data[starts[i]:ends[i]])/float(ends[i]-starts[i])
            threshold_windows.append(ThresholdViolation(typer, starts[i], ends[i], score))

        return threshold_windows

    def find_threshold(self, user_params):
        """Finds the thresholds for errors given the data using Gaussian Mixture Model

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
            
        """

        max_gauss_mixtures = user_params.get("max_gauss_mixtures")
        data = self.prob_smoothed

        print data

        # http://www.pymix.org/pymix/index.php?n=PyMix.Tutorial

        # make two gaussains
        gaussian_one = mixture.NormalDistribution(numpy.mean(data),numpy.std(data))
        gaussian_two = mixture.NormalDistribution(2.0*numpy.mean(data),numpy.std(data))

        mixture_model = mixture.MixtureModel(2, [0.5,0.5], [gaussian_one, gaussian_two])

        # print mixture_model

        EM_tuned = False
        while not EM_tuned:
            try:
                # make mix_data from a random 10% of the original data
                index_array = numpy.arange(data.size)
                numpy.random.shuffle(index_array)
                mix_data = mixture.DataSet()
                mix_data.fromArray(data[index_array[:int(numpy.floor(data.size/10.0))]])

                mixture_model.randMaxEM(mix_data, max_gauss_mixtures, 40, 0.001, silent=True)

                EM_tuned = True
            except AssertionError:
                # pymix likes to throw assertion errors when it has small machine precision errors...
                print "Caught an assertion error in pymix, randomizing input and trying again"

        #print mixture_model

        # hacky, no good api access to the model components
        gauss_one_mean = float(str(mixture_model.components[0][0]).split('[')[1].split(',')[0])
        gauss_one_std = float(str(mixture_model.components[0][0]).split(', ')[1].split(']')[0])

        gauss_two_mean = float(str(mixture_model.components[1][0]).split('[')[1].split(',')[0])
        gauss_two_std = float(str(mixture_model.components[1][0]).split(', ')[1].split(']')[0])

        #print "Gauss1: mu: %f, std: %f" % (gauss_one_mean, gauss_one_std)
        #print "Gauss2: mu: %f, std: %f" % (gauss_two_mean, gauss_two_std)

        #print "Using threshold %f" % threshold

        # inv normal cdf
        if gauss_one_mean > gauss_two_mean or mixture_model.pi[1] < 0.60:
            self.thresh_main_mean = gauss_one_mean
            self.thresh_main_std = gauss_one_std
        else:
            self.thresh_main_mean = gauss_two_mean
            self.thresh_main_std = gauss_two_std

        if self.thresh_main_std == 0.1:
            self.thresh_main_std = numpy.std(data)

class Contig():
    """A contig from an ALE assembly

    more info

    Attributes:
        name (str): The name of the contig (defaults to unnamed)

        length (int): The length of the contig

        depth (numpy.array): A numpy vector of depths for each position

        depth_prob (numpy.array): A numpy vector of depth probabilities for each position

        placement_prob (numpy.array): A numpy vector of placement probabilities for each position

        insert_prob (numpy.array): A numpy vector of insert probabilities for each position

        kmer_prob (numpy.array): A numpy vector of kmer probabilities for each position

        total_prob (numpy.array): A numpy vector of total probabilities for each position

    """
    def __init__(self, length=0, name="unnamed"):
        """Inits Contig with name and length (and prob vectors)
        
        Kwargs:
            length (int): The length of the contig (>=0)

            name (str): The name of the contig
            
        Raises:
            ValueError: length must be >= 0
        """

        if length < 0:
          raise ValueError("length must be >= 0")

        self.name = name
        self.length = length

        self.depth = numpy.zeros(length)

        self.prob_vecs={'d':LikelihoodVector(length=length, color='r'),
                   'p':LikelihoodVector(length=length, color='b'),
                   'i':LikelihoodVector(length=length, color='m'),
                   'k':LikelihoodVector(length=length, color='g')}

        self.main_figure = ALEFigure(0, length)

        self.pre_plot_run = False

    def pre_plot(self, user_params):
        """smooth the prob for the main figure and find the thresholds"""
        if self.pre_plot_run:
            print "Only need to run pre_plot() once! Ignoring..."
            return 0
        else:
            self.pre_plot_run = True
        # load parameters
        start = user_params.get("start")
        self.start = start
        end = user_params.get("end")
        self.end = end
        if end == 0:
            end = self.length

        # SMOOTHING
        # if there are no smoothing widths set them to a fraction of the window
        smoothing_dict = {'d':user_params.get("depth_smoothing_width"), 'p':user_params.get("placement_smoothing_width"), 'i':user_params.get("insert_smoothing_width"), 'k':user_params.get("kmer_smoothing_width")}

        largest_smooth = 0
        for typer in self.prob_vecs:
            if typer in user_params.get("plot_type"):
                print "smoothing %s" % typer
                self.prob_vecs[typer].smooth_prob(smoothing_width=smoothing_dict[typer])
                if self.prob_vecs[typer].smoothing_width > largest_smooth:
                    largest_smooth = self.prob_vecs[typer].smoothing_width

        if start < largest_smooth:
            self.start = largest_smooth
        if end > self.length - largest_smooth:
            self.end = self.length - largest_smooth

        self.largest_smooth = largest_smooth

        # THRESHOLDING
        for typer in user_params.get("plot_type"):
            print "Thresholding for %s" % typer
            self.prob_vecs[typer].find_threshold(user_params)
            print "mean, std", self.prob_vecs[typer].thresh_main_mean, self.prob_vecs[typer].thresh_main_std

    def sub_plots(self, user_params, start, end, pdf_stream = None):
        # make a sub plot
        increment = user_params.get("subplot_length")
        exists_plot = False
        plotted_previous = False
        if end == 0:
            end = self.length
        for s in numpy.arange(start, end, increment):

            sub_plot_start = s
            sub_plot_end = numpy.min((s+increment, end))
            if exists_plot and not self.main_figure.threshold_windows_in_last_half():
                if self.main_figure.threshold_windows_in_first_half() and plotted_previous:
                    plotted_previous = False
                else:
                    self.save_plot(user_params, pdf_stream=pdf_stream)
                    plotted_previous = True
            else:
                plotted_previous = False
            print "making figure %d-%d of %d" % (sub_plot_start, sub_plot_end, end)
            self.plot(user_params, sub_plot_start, sub_plot_end, is_subplot=True)
            exists_plot = True

            sub_plot_start = numpy.min((s+increment/2, end))
            sub_plot_end = numpy.min((s+3*increment/2, end))
            if exists_plot and not self.main_figure.threshold_windows_in_last_half():
                if self.main_figure.threshold_windows_in_first_half() and plotted_previous:
                    plotted_previous = False
                else:
                    self.save_plot(user_params, pdf_stream=pdf_stream)
                    plotted_previous = True
            else:
                plotted_previous = False
            print "making figure %d-%d of %d" % (sub_plot_start, sub_plot_end, end)
            self.plot(user_params, numpy.min((s+increment/2, end)), numpy.min((s+3*increment/2, end)), is_subplot=True)
            exists_plot = True

    def save_plot(self, user_params, pdf_stream=None):
        save_figure = user_params.get("save_figure")
        plot_threshold = user_params.get("plot_threshold")
        percent_thresholded = self.main_figure.get_percent_thresholded()

        if percent_thresholded > plot_threshold:
            if save_figure:
                pdf_stream.savefig()
            else:
                plt.show()
    
    def plot(self, user_params, start, end, is_subplot=False):
        if not self.pre_plot_run:
            print "Running pre_plot..."
            self.pre_plot(user_params)

        # load parameters
        if end == 0:
            end = self.length

        # make a new figure
        self.main_figure = ALEFigure(start, end)

        plot_type = user_params.get("plot_type")
        min_plot_size = user_params.get("min_plot_size")

        if end-start < min_plot_size:
            print "lower than min plot size"
            return 0, []

        if is_subplot:
            thresh_mult = user_params.get("sub_threshold_depth")
            std_width = int(-2.0*thresh_mult/5.0)
        else:
            thresh_mult = user_params.get("threshold_depth")
            std_width = int(-2.0*thresh_mult/5.0)
        
        # Main plotting code

        # THRESHOLD WINDOWS
        plot_type = list(plot_type)
        plot_type.sort() # make sure t comes last, if at all
        for typer in plot_type:
            if not is_subplot:
                threshold_windows = self.prob_vecs[typer].get_threshold_windows(self.prob_vecs[typer].prob_smoothed[start:end], user_params, typer=typer, thresh_mult=thresh_mult)
            else:
                threshold_windows = self.prob_vecs[typer].get_threshold_windows(self.prob_vecs[typer].prob[start:end], user_params, typer=typer, thresh_mult=thresh_mult)
            self.main_figure.add_threshold_windows(threshold_windows)
        self.main_figure.plot_threshold_windows()
        percent_thresholded = self.main_figure.get_percent_thresholded()

        # DATA LINES
        current_subplot = 0
        for typer in plot_type:
            if not is_subplot:
                self.main_figure.ax.plot(self.main_figure.format_data_for_plot(current_subplot, self.prob_vecs[typer].thresh_main_mean, self.prob_vecs[typer].thresh_main_std, self.prob_vecs[typer].prob_smoothed, start, end), self.prob_vecs[typer].color)
            else:
                self.main_figure.ax.plot(self.main_figure.format_data_for_plot(current_subplot, self.prob_vecs[typer].thresh_main_mean, self.prob_vecs[typer].thresh_main_std, self.prob_vecs[typer].prob, start, end, std_width=std_width), self.prob_vecs[typer].color)
            self.main_figure.ax.plot([0, end - start], [4 + 7*current_subplot - 2.5,4 + 7*current_subplot - 2.5], 'black')
            self.main_figure.plot_std_marks(current_subplot, end - start, self.prob_vecs[typer].color)
            current_subplot += 1

        # LABELS, TITLES, ETC
        self.main_figure.set_labels(self.main_figure.ax2, plot_type, self.prob_vecs, twin=True, std_width=std_width)
        self.main_figure.set_labels(self.main_figure.ax, plot_type, self.prob_vecs, std_width=std_width)

        self.main_figure.ax.set_title('Average Likelihoods (' + str(percent_thresholded)[:6] + ' below threshold)')
        self.main_figure.ax.set_xlabel('Position (bp) (+' + str(numpy.max([start, self.largest_smooth])) + ')')
        self.main_figure.plot_cleanup(start, end, current_subplot)

        return percent_thresholded, self.main_figure.threshold_window_set

def plot_histogram(input_data, save_figure=False, pdf_stream=None):
    max_val = numpy.max(input_data)
    min_val = numpy.min(input_data)
    bin_size = (max_val - min_val)/100.0
    if bin_size < 1e-6: bin_size = 0.01
    histogram = numpy.zeros(100)
    for value in input_data:
        histogram[min((max((0,int(numpy.floor((value - min_val)/bin_size)))),99))] += 1

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.hist(histogram, 100, normed=1, facecolor='green', alpha=0.75)

    ax.set_xlabel('threshold')
    ax.set_ylabel('distribution')

    if save_figure:
        pdf_stream.savefig()
    else:
        plt.show()

def read_in_info(placement_file):
    """Reads in an ALE placement file, returns a list of Contigs.

    Args:
        placement_file: An ALE placement file (*.ale)
            must be in the following format::

                # comments/metadata
                # can have multiple lines, all starting with #
                # Reference: gi|170079663|ref|NC_010473.1| 350000
                # contig position depth ln(depthLike) ln(placeLike) ln(insertLike) ln(kmerLike) 
                0 0 1.000000 -60.000000 0.194888 -5.760798 -65.565910
                0 1 3.000000 -60.000000 0.466271 -5.608334 -65.142063
                0 2 5.000000 -60.000000 0.010585 -5.541655 -65.531071
                0 3 12.000000 -60.000000 -0.057731 -5.380759 -65.438491

            Specific lines (using the above as an example):
                0. Any number of comment lines starting with #
                1. The length of the contig is::

                       length = int(line.split(' ')[3]) == 350000

                   The name of the contig is::

                       name = line.split(' ')[2] == gi|170079663|ref|NC_010473.1|

                   name **cannot** be 'position'

                2. The following line is ignored and lists what is in the columns of following lines

                3. The data corresponding to the column headers for each position in the contig

                4. See 2.

    Returns:
        A list of Contigs (see class :py:mod:`plotter3.Contig`)

    Raises:
        IOError: An error occured accessing the placement file.

        FormattingError: The placement file was not formatted correctly.
    """

    MINIMUM_VALUE = -60.0 # minimum value we allow a position to have

    ale_placement_file = open(placement_file, 'r')

    contigs = []
    previous_line_one = ""
    previous_line_two = ""

    for line in ale_placement_file:
        if line[0] == '#':
            if previous_line_one == "":
                previous_line_one = line
            else:
                previous_line_two = previous_line_one
                previous_line_one = line           
        else:
            if previous_line_two != "":
                tName = previous_line_two.split(' ')[2]              
                tLen = int(previous_line_two.split(' ')[-1])
                contigs.append(Contig(tLen, name = tName))
                place = 0
                print "Reading in contig: " + tName + " len " + str(tLen)
                print ""
                bar = ProgressBar.progressBar(0, tLen, 42)
                previous_line_two = ""
            data = line.split(' ')
            contigs[-1].depth[place] = numpy.double(data[2])
            for i in range(1,7):
                if "-nan"==data[i] or "nan"==data[i] or "inf"==data[i] or "-inf"==data[i] or numpy.double(data[i]) != numpy.double(data[i]):
                    data[i] = MINIMUM_VALUE # Predefined threshold
            contigs[-1].prob_vecs['d'].prob[place] = numpy.double(data[3])
            contigs[-1].prob_vecs['p'].prob[place] = numpy.double(data[4])
            contigs[-1].prob_vecs['i'].prob[place] = numpy.double(data[5])
            contigs[-1].prob_vecs['k'].prob[place] = numpy.double(data[6])           
            place += 1
            if tLen > 40:
                if (place)%(int(tLen)/40) == 0:
                    bar(place)
    print "\nYou now have a list of contigs, try contig[0].plot()"
    return contigs

def main():
    # default parameter values
    user_params = CommandLineParameters()
    user_params.add_parameter("start", "-s", int, 0)
    user_params.add_parameter("end", "-e", int, 0)
    user_params.add_parameter("plot_type", "-pt", str, "idpk")
    user_params.add_parameter("save_figure", "-nosave", None, True)
    user_params.add_parameter("depth_smoothing_width", "-dsw", int, 0)
    user_params.add_parameter("placement_smoothing_width", "-psw", int, 0)
    user_params.add_parameter("insert_smoothing_width", "-isw", int, 0)
    user_params.add_parameter("kmer_smoothing_width", "-ksw", int, 0)
    user_params.add_parameter("threshold_depth", "-td", float, -5.0)
    user_params.add_parameter("threshold_percent", "-tp", float, 0.2)
    user_params.add_parameter("threshold_width", "-tw", int, 1)
    user_params.add_parameter("sub_threshold_depth", "-std", float, -30.0)
    user_params.add_parameter("subplot_length", "-sl", int, 5000)
    user_params.add_parameter("plot_threshold", "-plt", float, 0.0)
    user_params.add_parameter("figure_name", "-fn", str, "")
    user_params.add_parameter("plot_meta", "-dpm", None, True)
    user_params.add_parameter("plot_meta_only", "-pmo", None, False)
    user_params.add_parameter("specific_contig", "-sc", str, None)
    user_params.add_parameter("min_plot_size", "-mps", int, 100)
    user_params.add_parameter("max_gauss_mixtures", "-mgm", int, 2)
    user_params.add_parameter("N_worst_positions", "-nwp", int, 10)
    user_params.add_input("ale_file")

    # read in command line arguments
    user_params.read_sys_args()

    # read in contigs
    contigs = read_in_info(user_params.get_input("ale_file"))

    save_figure = user_params.get("save_figure")
    # open up a pdf_stream and file for output
    if save_figure:
        figure_name = user_params.get("figure_name")
        if figure_name == "":
            figure_name = sys.argv[-1] + '.pdf'
        pdf_stream = PdfPages(figure_name)
    else:
        pdf_stream = None

    print "Generating figures..."

    meta_percent = []
    meta_number = []

    fout = open(sys.argv[-1] + ".thresholds", 'w')

    # generate the data, given the user params
    if not user_params.get("plot_meta_only"):
        for contig in contigs:
            if contig.length >= user_params.get("min_plot_size"):
                if not user_params.get("specific_contig") or user_params.get("specific_contig") == contig.name:
                    contig.pre_plot(user_params)
                    percent_thresholded, threshold_windows = contig.plot(user_params, user_params.get("start"), user_params.get("end"))
                    contig.save_plot(user_params, pdf_stream=pdf_stream)
                    contig.sub_plots(user_params, user_params.get("start"), user_params.get("end"), pdf_stream=pdf_stream)
                    print "%s had (%f) thresholded." % (contig.name, percent_thresholded)
                    meta_percent.append(percent_thresholded)
                    for window in threshold_windows:
                        fout.write("%s:%i-%i\t%s\n" % (contig.name, window.start, window.end, window.type_of))

    fout.close()

    
    # plot the meta data
    if len(meta_percent) > 1:
        if user_params.get("plot_meta") or user_params.get("plot_meta_only"):
            plot_histogram(meta_percent, save_figure=save_figure, pdf_stream=pdf_stream)
            plot_histogram(meta_number, save_figure=save_figure, pdf_stream=pdf_stream)

    # save and close the output
    if save_figure:
        print "saved file %s" % figure_name
        pdf_stream.close()

    print "Executed Sucessfully"

if __name__ == '__main__':
    main()
