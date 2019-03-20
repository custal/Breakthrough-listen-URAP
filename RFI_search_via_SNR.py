import numpy as np
import blimpy as bp
import pylab as plt
from blimpy.utils import db, lin, rebin, closest, unpack_2to8
from blimpy.sigproc import *

def get_spikes(arr, period):
    """Finds the indices of the DC spikes of a GBT spectrum based on a priori
    knowledge of the periodicity of the power spectrum, assuming spikes occur
    in the center of each coarse channel. e.g., if the period is 1024 points,
    then the half-period is 512, and spikes occur at indices i={511, 511+1024,
    511+2*(1024), ...}.

    Args:
    -----
    arr (1-D iterable): data
    period (int): number of bins in a coarse channel
            (e.g. 1024 for mid-resolution product)

    Returns:
    --------
    spikes (np.array): List of indices at which the spikes occur.
    """

    last_index = len(arr)-1
    #Spikes are assumed to start half a period from beginning
    #of data and repeat with the given periodicity.
    spikes = np.arange(0, last_index, period)
    return spikes.astype(int)


def remove_RFI(file, t=0):
    """Removes potential RFI signals from data by removing data points with a
    high standard deviation above the mean signal to noise ratio.

    Args:
    -----
    file (waterfall): Waterfall object to have data removed from.
    t (int): integration of data to take

    Returns:
    --------
    waterfall object with RFI removed
    """

    #find spike points defining coarse channels
    chans_per_coarse = int(file.file_header[b'nchans'])/64
    spikes = get_spikes(file.data[0][0], chans_per_coarse)

    freqs = file.populate_freqs()
    if file.header[b'foff'] < 0:
        freqs = freqs[::-1]

    coarse_channels = []
    #Make array of coarse channels (sorts only first integration at the moment)
    for i in range(len(spikes)-1):
         package = []
         coarse_chan_f, coarse_chan_data = file.grab_data(freqs[spikes[i]], freqs[spikes[i+1]])
         package.append(coarse_chan_f)
         package.append(db(coarse_chan_data[t]))
         coarse_channels.append(package)


    channel = 60
    #Remove RFI (29)
    #for coarse in coarse_channels:
    coarse = coarse_channels[channel]
    coarse_chan_data = coarse[1]
    coarse_chan_data_sorted = np.sort(coarse_chan_data)

    #remove top 10% and bottom 20% of array, this cuts RFI and FFT spikes
    coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[int(chans_per_coarse*0.9):int(chans_per_coarse+1)])
    coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[0:int(chans_per_coarse*0.2)])

    std = np.std(coarse_chan_data_sorted)
    print(std)
    #mask values that maybe RFI above a specified standard deviation
    median = np.median(coarse_chan_data)
    print(median)
    coarse = np.ma.array(coarse_chan_data)
    coarse = np.ma.masked_where((coarse_chan_data - median)/std > 5 , coarse_chan_data)

    median_array = [median] * len(coarse_channels[channel][0])
    plt.plot(coarse_channels[channel][0], median_array)
    plt.plot(coarse_channels[channel][0], coarse)


#used to test to view spikes
# x = []
# y = []
# a = 60
# for i in range(2048):
#     print(plot_f[int(spikes[a]+i)], plot_data[int(spikes[a]+i)])
#     if i != 0:
#         x.append(plot_f[int(spikes[a]+i)])
#         y.append(plot_data[int(spikes[a]+i)])
#     else:
#         print("^")
#         x.append(plot_f[int(spikes[a]+i)])
#         y.append(plot_data[int(spikes[a]+i)])
#         x0 = plot_f[int(spikes[a]+i)]
#         y0 = plot_data[int(spikes[a]+i)]
#
#
# plt.plot(x, y)
# plt.scatter(x0, y0, s=50)


fil = bp.Waterfall("blc00_guppi_57872_20242_DIAG_2MASS_1502+2250_0024.gpuspec.0002.fil")
remove_RFI(fil)
plt.show()
