"""
#RFI_search_via_SNR.py

Find RFI by excluding any points with a large signal to noise ratio.

TODO: Once a point has been removed mark its frequency value.

"""

import numpy as np
import blimpy as bp
import pylab as plt
from blimpy.utils import db, lin, rebin, closest, unpack_2to8
from blimpy.sigproc import *



#np.set_printoptions(threshold=np.inf)

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

    last_index = len(arr)+1
    #Spikes are assumed to start half a period from beginning
    #of data and repeat with the given periodicity.
    spikes = np.arange(0, last_index, period)
    return spikes.astype(int)

def grab_data_index(array, start_index, end_index, direction):
    """Grabs data between specified indexes.

    Args:
    -----
    array (array of floats): Array to grab data from.
    start_index (int): Grab data starting from this index.
    end_index (int): Grab data ending at this index.
    direction (float): direction of data (file.header[b'foff'])

    Returns:
    --------
    Data between specified indexes
    """

    data = []
    for i in range(start_index, end_index+1):
        data.append(array[i])

    if direction < 0:
        data[::-1]

    return data


def remove_RFI(file, t=0):
    """Removes potential RFI signals from data by removing data points with a
    high standard deviation above the mean signal to noise ratio.

    Args:
    -----
    file (waterfall): Waterfall object to have data removed from.
    t (int): time integration of data to take

    Returns:
    --------
    waterfall object with RFI removed

    TODO: make std threshold and number of coarse channels variables (and %
    data cut above and below)
    TODO: number of coarse channels, 64, is hard coded in, any way to discern
    this from file header?
    """

    #find spike points defining coarse channels
    chans_per_coarse = int(file.file_header[b'nchans'])/64
    spikes = get_spikes(file.data[t][0], chans_per_coarse)

    n_coarse_chan = file.calc_n_coarse_chan()
    file.blank_dc(n_coarse_chan)

    coarse_channels = []
    #Make array of coarse channels (sorts only first integration at the moment)
    for i in range(len(spikes)-1):
         coarse_chan_data = grab_data_index(file.data[t][0], spikes[i], spikes[i+1]-1, file.header[b'foff'])
         coarse_channels.append(db(coarse_chan_data))

    coarse_channels_masked = []
    for coarse in coarse_channels:
        coarse_chan_data_sorted = np.sort(coarse)

        #remove top 10% and bottom 20% of array, this cuts RFI and FFT spikes
        coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[int(chans_per_coarse*0.9):int(chans_per_coarse+1)])
        coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[0:int(chans_per_coarse*0.2)])

        std = np.std(coarse_chan_data_sorted)
        median = np.median(coarse_chan_data_sorted)

        #mask values that maybe RFI above a specified standard deviation
        coarse_masked = np.ma.array(coarse)
        coarse_masked = np.ma.masked_where((coarse - median)/std > 10 , coarse)
        coarse_channels_masked.append(coarse_masked)


    all_channels = np.ma.concatenate(coarse_channels_masked)
    file.data[t][0] = all_channels

    freqs = file.populate_freqs()
    if file.header[b'foff'] < 0:
        freqs[::-1]

    plt.plot(freqs, all_channels)
    return
    # test_channels = []
    # for i in range(24,25):
    #     test_channels.append(coarse_channels_masked[i])
    #
    # test_channels = np.ma.concatenate(test_channels)
    # x = range(1024)
    # plt.plot(x, test_channels)
    #
    # plt.subplot(2, 1, 1)
    # plt.plot(freqs, db(file.data[t][0]))
    # plt.title('RFI search via SNR')
    #
    # plt.subplot(2, 1, 2)
    # plt.plot(freqs, all_channels)

def mark_RFI(file, t=0):
    """Similar to above function but marks index of RFI instead of removing it.
    Additional feature: finds the RFI from the median spectrum.

    Args:
    -----
    file (waterfall): Waterfall object to have data removed from.
    t (int): time integration of data to take

    Returns:
    --------
    array of indexes which are marked as being at RFI frequencies.

    TODO: make std threshold and number of coarse channels variables (and %
    data cut above and below)
    TODO: number of coarse channels, 64, is hard coded in, any way to discern
    this from file header?
    """

    #find spike points defining coarse channels
    chans_per_coarse = int(file.file_header[b'nchans'])/64
    spikes = get_spikes(file.data[t][0], chans_per_coarse)

    n_coarse_chan = file.calc_n_coarse_chan()
    file.blank_dc(n_coarse_chan)
    data = np.median(file.data, axis=0)[0]
    #data = file.data[0][0]

    coarse_channels = []
    #Make array of coarse channels (sorts only first integration at the moment)
    for i in range(len(spikes)-1):
         coarse_chan_data = grab_data_index(data, spikes[i], spikes[i+1]-1, file.header[b'foff'])
         coarse_channels.append(db(coarse_chan_data))

    marked_indexes = []
    for i, coarse in enumerate(coarse_channels):
        coarse_chan_data_sorted = np.sort(coarse)

        #remove top 10% and bottom 20% of array, this cuts RFI and FFT spikes
        coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[int(chans_per_coarse*0.9):int(chans_per_coarse+1)])
        coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[0:int(chans_per_coarse*0.2)])

        std = np.std(coarse_chan_data_sorted)
        median = np.median(coarse_chan_data_sorted)

        #return indexes that maybe RFI above a specified standard deviation
        coarse_normalized = (coarse - median)/std
        coarse_marked = np.argwhere(coarse_normalized > 10)

        #create list of marked indexes
        if len(coarse_marked) != 0:
            coarse_marked = coarse_marked.flatten()
            coarse_marked = coarse_marked + spikes[i]
            marked_indexes.append(coarse_marked)

    marked_indexes = np.concatenate(marked_indexes)
    return marked_indexes


def mark_RFI_aggregated(file, aggregate=1, t=0):
    """Similar to above function but marks index of RFI instead of removing it.
    Additional feature: finds the RFI from the median spectrum but aggregated
    into blocks.

    Args:
    -----
    file (waterfall): Waterfall object to have data removed from.
    aggregate (int): Number to divide data in the time axis. Larger number makes
    RFI search more refined and more likely to find transient RFI.
    t (int): time integration of data to take.

    Returns:
    --------
    array of indexes which are marked as being at RFI frequencies.

    TODO: make std threshold and number of coarse channels variables (and %
    data cut above and below)
    TODO: number of coarse channels, 64, is hard coded in, any way to discern
    this from file header?
    """

    #find spike points defining coarse channels
    chans_per_coarse = int(file.file_header[b'nchans'])/64
    spikes = get_spikes(file.data[t][0], chans_per_coarse)

    n_coarse_chan = file.calc_n_coarse_chan()
    file.blank_dc(n_coarse_chan)

    #aggregate data to hopefully detect more RFI
    data_split = np.array_split(file.data, aggregate, axis=0)

    marked_indexes = []
    for data in data_split:
        data = np.median(data, axis=0)[0]

        coarse_channels = []
        #Make array of coarse channels (sorts only first integration at the moment)
        for i in range(len(spikes)-1):
             coarse_chan_data = grab_data_index(data, spikes[i], spikes[i+1]-1, file.header[b'foff'])
             coarse_channels.append(db(coarse_chan_data))

        for i, coarse in enumerate(coarse_channels):
            coarse_chan_data_sorted = np.sort(coarse)

            #remove top 10% and bottom 20% of array, this cuts RFI and FFT spikes
            coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[int(chans_per_coarse*0.9):int(chans_per_coarse+1)])
            coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[0:int(chans_per_coarse*0.2)])

            std = np.std(coarse_chan_data_sorted)
            median = np.median(coarse_chan_data_sorted)

            #return indexes that maybe RFI above a specified standard deviation
            coarse_normalized = (coarse - median)/std
            coarse_marked = np.argwhere(coarse_normalized > 10)

            #create list of marked indexes
            if len(coarse_marked) != 0:
                coarse_marked = coarse_marked.flatten()
                coarse_marked = coarse_marked + spikes[i]
                marked_indexes.append(coarse_marked)

    marked_indexes = np.concatenate(marked_indexes)
    marked_indexes = np.unique(marked_indexes)
    return marked_indexes


def frequency_cut(file, f_index):
    """ Cut entire frequency column for a given index.

    Args:
        file (waterfall): waterfall object to cut data from.
        f_index (arrray of ints): Array of frequency indexes to be cut.


    Returns:
        (file): instance of Waterfall with cut data
    """

    #make mask
    mask = np.zeros(file.data.shape, dtype=bool)
    for index in f_index:
        for i in range(len(file.data)):
            mask[i][0][index] = True

    file.data = np.ma.array(file.data, mask=mask)
    return file


fil = bp.Waterfall("blc00_guppi_57872_20242_DIAG_2MASS_1502+2250_0024.gpuspec.0002.fil")
marked_indexes = mark_RFI_aggregated(fil, 100)
fil = frequency_cut(fil, marked_indexes)
fil.plot_waterfall()

plt.show()
