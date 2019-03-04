"""
#nanmean_vs_mean_time.py

Test for difference in computation time, if any, between numpy.mean() and numpy.nanmean() in the plot_time_series() function from filterbank module

"""

import blimpy as bp
import pylab as plt
import numpy as np
from blimpy.utils import db, lin, rebin, closest, unpack_2to8
import time


def plot_time_series_nanmean(filename, f_start=None, f_stop=None, if_id=0, logged=True, orientation='h', MJD_time=False, **kwargs):
    """ Plot the time series.

     Args:
        filename (string): name of file to run test on
        f_start (float): start frequency, in MHz
        f_stop (float): stop frequency, in MHz
        logged (bool): Plot in linear (False) or dB units (True),
        kwargs: keyword args to be passed to matplotlib imshow()
    """
    start = time.time()
    fil = bp.Waterfall(filename, f_start=f_start, f_stop=f_stop)

    ax = plt.gca()
    plot_f, plot_data = fil.grab_data(f_start, f_stop, if_id)

    if logged and fil.header[b'nbits'] >= 8:
        plot_data = db(plot_data)

    #Since the data has been squeezed, the axis for time goes away if only one bin, causing a bug with axis=1
    if len(plot_data.shape) > 1:
        plot_data = np.nanmean(plot_data, axis=1)
    else:
        plot_data = np.nanmean(plot_data)

    #Make proper time axis for plotting (but only for plotting!). Note that this makes the values inclusive.
    extent = fil._calc_extent(plot_f=plot_f,plot_t=fil.timestamps,MJD_time=MJD_time)
    plot_t = np.linspace(extent[2],extent[3],len(fil.timestamps))

    if MJD_time:
        tlabel = "Time [MJD]"
    else:
        tlabel = "Time [s]"

    if logged:
        plabel = "Power [dB]"
    else:
        plabel = "Power [counts]"

    # Reverse oder if vertical orientation.
    if 'v' in orientation:
        plt.plot(plot_data, plot_t, **kwargs)
        plt.xlabel(plabel)

    else:
        plt.plot(plot_t, plot_data, **kwargs)
        plt.xlabel(tlabel)
        plt.ylabel(plabel)

    ax.autoscale(axis='both',tight=True)
    end = time.time()
    print("Time for np.nanmean: %.2f" % (end - start))

def plot_time_series_mean(filename, f_start=None, f_stop=None, if_id=0, logged=True, orientation='h', MJD_time=False, **kwargs):
    """ Plot the time series.

     Args:
        filename (string): name of file to run test on
        f_start (float): start frequency, in MHz
        f_stop (float): stop frequency, in MHz
        logged (bool): Plot in linear (False) or dB units (True),
        kwargs: keyword args to be passed to matplotlib imshow()
    """

    start = time.time()
    fil = bp.Waterfall(filename, f_start=f_start, f_stop=f_stop)

    ax = plt.gca()
    plot_f, plot_data = fil.grab_data(f_start, f_stop, if_id)

    if logged and fil.header[b'nbits'] >= 8:
        plot_data = db(plot_data)

    #Since the data has been squeezed, the axis for time goes away if only one bin, causing a bug with axis=1
    if len(plot_data.shape) > 1:
        plot_data = np.mean(plot_data, axis=1)
    else:
        plot_data = np.mean(plot_data)

    #Make proper time axis for plotting (but only for plotting!). Note that this makes the values inclusive.
    extent = fil._calc_extent(plot_f=plot_f,plot_t=fil.timestamps,MJD_time=MJD_time)
    plot_t = np.linspace(extent[2],extent[3],len(fil.timestamps))

    if MJD_time:
        tlabel = "Time [MJD]"
    else:
        tlabel = "Time [s]"

    if logged:
        plabel = "Power [dB]"
    else:
        plabel = "Power [counts]"

    # Reverse oder if vertical orientation.
    if 'v' in orientation:
        plt.plot(plot_data, plot_t, **kwargs)
        plt.xlabel(plabel)

    else:
        plt.plot(plot_t, plot_data, **kwargs)
        plt.xlabel(tlabel)
        plt.ylabel(plabel)

    ax.autoscale(axis='both',tight=True)
    end = time.time()
    print("Time for np.mean: %.2f" % (end - start))


plot_time_series_nanmean("Voyager1.single_coarse.fine_res.fil")
plot_time_series_mean("Voyager1.single_coarse.fine_res.fil")
