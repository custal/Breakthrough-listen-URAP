import numpy as np
import blimpy as bp
import pylab as plt
import copy
from RFI_search_via_SNR import *

MAX_IMSHOW_POINTS   = (8192, 4096)

def plot_waterfall_clean_six(fil, f_start=None, f_stop=None, if_id=0, logged=True, MJD_time=False, **kwargs):
    """ Plot waterfall of data

    Args:
        f_start (float): start frequency, in MHz
        f_stop (float): stop frequency, in MHz
        logged (bool): Plot in linear (False) or dB units (True),
        cb (bool): for plotting the colorbar
        kwargs: keyword args to be passed to matplotlib imshow()
    """

    plot_f, plot_data = fil.grab_data(f_start, f_stop, if_id)

    #Using accending frequency for all plots.
    if fil.header[b'foff'] < 0:
        plot_data = plot_data[..., ::-1] # Reverse data
        plot_f = plot_f[::-1]

    if logged:
        plot_data = db(plot_data)

    # Make sure waterfall plot is under 4k*4k
    dec_fac_x, dec_fac_y = 1, 1
    if plot_data.shape[0] > MAX_IMSHOW_POINTS[0]:
        dec_fac_x = int(plot_data.shape[0] / MAX_IMSHOW_POINTS[0])

    if plot_data.shape[1] > MAX_IMSHOW_POINTS[1]:
        dec_fac_y =  int(plot_data.shape[1] /  MAX_IMSHOW_POINTS[1])

    plot_data = rebin(plot_data, dec_fac_x, dec_fac_y)

    extent = fil._calc_extent(plot_f=plot_f,plot_t=fil.timestamps,MJD_time=MJD_time)

    return plot_data, extent


raw_files = [bp.Waterfall("spliced_blc0001020304050607_guppi_57660_63677_HIP58684_0039.gpuspec.0002.h5"),
bp.Waterfall("spliced_blc0001020304050607_guppi_57660_64028_HIP57601_0040.gpuspec.0002.h5"),
bp.Waterfall("spliced_blc0001020304050607_guppi_57660_64379_HIP58684_0041.gpuspec.0002.h5"),
bp.Waterfall("spliced_blc0001020304050607_guppi_57660_64723_HIP57780_0042.gpuspec.0002.h5"),
bp.Waterfall("spliced_blc0001020304050607_guppi_57660_65067_HIP58684_0043.gpuspec.0002.h5"),
bp.Waterfall("spliced_blc0001020304050607_guppi_57660_65408_HIP57845_0044.gpuspec.0002.h5")]


for fil in raw_files:
    plt.figure()
    f, axarr = plt.subplots(2, sharex=True)
    f.suptitle('Comparison of raw data to cleaned data')
    plot_data, extent = plot_waterfall_clean_six(fil)
    axarr[0].imshow(plot_data,
        aspect='auto',
        origin='lower',
        rasterized=True,
        interpolation='nearest',
        extent=extent,
        cmap='viridis'
    )

    marked_indices = mark_RFI(fil, 100, 0.3, 0.7, 10, coarse_chan_num=10)
    marked_indices = add_freq_range(fil, marked_indices, 2175, 2205)
    marked_indices = add_freq_range(fil, marked_indices, 2307, 2364)
    fil = frequency_cut(fil, marked_indices, fill_mask=True)
    print("hello world")

    plot_data, extent = plot_waterfall_clean_six(fil)
    axarr[1].imshow(plot_data,
        aspect='auto',
        origin='lower',
        rasterized=True,
        interpolation='nearest',
        extent=extent,
        cmap='viridis'
    )


plt.show()
