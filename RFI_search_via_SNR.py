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
#from spec_fit import _smart_split



#np.set_printoptions(threshold=np.inf)

class Marked_index(object):
    """Class for holding data on indices to later be cut and the mask to replace
    them with. Zeroe by default, set to median value of coarse channel by
    mark_RFI or assigned when calling add_indices or add_freq_range  """

    def __init__(self, indices, mask=0):
        self.indices = indices
        self.mask = mask



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
        coarse_masked = np.ma.masked_where((coarse - median)/std > 5 , coarse)
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


def mark_RFI(file, aggregate=1, central_lower=0.2, central_upper=0.8, sigma=10, t=0):
    """Marks potential RFI signals from data by removing data points with a high
    standard deviation above the mean signal to noise ratio.

    Args:
    -----
    file (waterfall): Waterfall object to have data removed from.
    aggregate (int): Number to divide data in the time axis. Larger number makes
    RFI search more refined and more likely to find transient RFI. Must be less
    than number of time indices.
    central_lower (flaot): Lower bound of data (as a fraction of total data) to
    keep when calculating the central median of a coarse channel. eg. 0.3
    removes bottom 30% of data before median calculation.
    central_upper (float): Upper bound of data (as a fraction of total data) to
    keep when calculating the central median of a coarse channel. eg. 0.7
    removes top 30% of data before median calculation.
    sigma (float): Remove any data points with a value greater than sigma standard
    deviations away from the central median
    t (int): time integration of data to take.

    Returns:
    --------
    array of indexes which are marked as being at RFI frequencies.

    TODO: number of coarse channels, 64, is hard coded in, any way to discern
    this from file header? Maybe use Umran's _smart_split function
    """

    #check to make sure aggregate is less than number of time indices in the file.
    a = []
    time_indices = len(file.data)
    if aggregate > time_indices:
        print("Aggregate greater than number of time indices in file, %i, please choose a smaller or equal number." % (time_indices))
        return a
    #make sure central_lower < central_upper
    if central_lower > central_upper:
        print("Lower_central is greater than upper_central, please use new values")
        return a

    #find spike points defining coarse channels
    chans_per_coarse = int(file.file_header[b'nchans'])/64
    spikes = get_spikes(file.data[t][0], chans_per_coarse)

    n_coarse_chan = file.calc_n_coarse_chan()
    file.blank_dc(n_coarse_chan)

    #aggregate data to hopefully detect more RFI
    data_split = np.array_split(file.data, aggregate, axis=0)

    marked_indices = []
    for data in data_split:
        data = np.median(data, axis=0)[0]

        coarse_channels = []
        #Make array of coarse channels
        for i in range(len(spikes)-1):
             coarse_chan_data = grab_data_index(data, spikes[i], spikes[i+1]-1, file.header[b'foff'])
             coarse_channels.append(db(coarse_chan_data))

        for i, coarse in enumerate(coarse_channels):
            coarse_chan_data_sorted = np.sort(coarse)

            #remove top % and bottom % of array, this cuts RFI and FFT spikes
            coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[int(chans_per_coarse*central_upper):int(chans_per_coarse+1)])
            coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[0:int(chans_per_coarse*central_lower)])

            std = np.std(coarse_chan_data_sorted)
            median = np.median(coarse_chan_data_sorted)

            #return indexes that maybe RFI above a specified standard deviation
            coarse_normalized = (coarse - median)/std
            coarse_marked = np.argwhere(coarse_normalized > sigma)

            #create list of marked indexes
            if len(coarse_marked) != 0:
                coarse_marked = coarse_marked.flatten()
                coarse_marked = coarse_marked + spikes[i]
                marked_indices.append(coarse_marked)


    if len(marked_indices) != 0:
       marked_indices = np.concatenate(marked_indices)
       marked_indices = np.unique(marked_indices)

    marked_indices = marked_indices.tolist()
    return marked_indices


def frequency_cut(file, f_index, fill_median=False):
    """ Cut entire frequency column for a given index.

    Args:
        file (waterfall): waterfall object to cut data from.
        f_index (arrray of ints): Array of frequency indexes to be cut.
        fill_median (bool): If True then masked RFI will be replaced by the
        median value of the coarse channel without RFI or by whatever value is
        assigned to the index to be cut by the add_indices and add_freq_range
        functions.

    Returns:
        (file): instance of Waterfall with cut data
    """

    #make mask
    mask = np.zeros(file.data.shape, dtype=bool)
    for index in f_index:
        for i in range(len(file.data)):
            mask[i][0][index] = True

    file.data = np.ma.array(file.data, mask=mask)

    if fill_median == True:
        median = file.data
        file.data = np.ma.set_fill_value()
    return file

def add_indices(marked_indices, *frequencies):
    """ Add more indices to be cut to the marked_indices array.

    Args:
        marked_indices (list of ints): Array of marked indices .
        *frequencies (ints or list of ints): Frequency indices to be cut. Takes list or single values

    Returns:
        (array of ints): Marked indices with new frequencies added
    """

    for freq in frequencies:
        flat_list = []
        if type(freq) is list:
            for item in freq:
                 marked_indices.append(item)
        else:
            marked_indices.append(freq)

    marked_indices = np.unique(marked_indices)

    return marked_indices

def frequency_to_index(frequency, f_increment, f0):
    return int(np.round((frequency - f0) / f_increment))

def add_freq_range(fil, marked_indices, frequency_lower, frequency_upper):
    """ Add an entire frequency range of indices to marked_indices.

    Args:
        fil (waterfall): File working on.
        marked_indices (list of ints): Array of marked indices.
        frequency_lower (float): Lower bound of frequency range to be cut.
        frequency_upper (float): Upper bound of frequency range to be cut.

    Returns:
        (array of ints): Marked indices with new frequencies added
    """
    f_increment = fil.file_header[b'foff']
    f0 = fil.container.f_stop
    if fil.file_header[b'foff'] < 0:
        a = frequency_lower
        frequency_lower = frequency_upper
        frequency_upper = a

    frequency_lower = frequency_to_index(frequency_lower, f_increment, f0)
    frequency_upper = frequency_to_index(frequency_upper, f_increment, f0)
    frequency_range = range(frequency_lower, frequency_upper)

    for item in frequency_range:
        marked_indices.append(item)

    marked_indices = np.unique(marked_indices)
    return marked_indices

#blc00_guppi_57872_20242_DIAG_2MASS_1502+2250_0024.gpuspec.0002.fil
#spliced_blc0001020304050607_guppi_57936_37003_HIP116719_0057.gpuspec.0002.fil Presence of time spike due to RFI filling entire coarse channel
fil = bp.Waterfall("spliced_blc0001020304050607_guppi_57936_37003_HIP116719_0057.gpuspec.0002.fil")
marked_indices = mark_RFI(fil, 100, 0.3, 0.7, 10)
marked_indices = add_freq_range(fil, marked_indices, 2178, 2202)
fil = frequency_cut(fil, marked_indices)
fil.plot_waterfall(logged=True)
# remove_RFI(fil)

plt.show()
