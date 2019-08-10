"""
#RFI_search_via_SNR.py

Find RFI by excluding any points with a large signal to noise ratio.


"""

import numpy as np
import blimpy as bp
import bisect
from collections import OrderedDict
from blimpy.utils import db, lin, rebin, closest, unpack_2to8
from blimpy.sigproc import *
import time


class RFI_cleaner(object):
    """Class which contains all indices to be marked for removal/masking, the
    masks to use and all methods to acquire these indices and masks
    """

    def __init__(self, waterfall_object: bp.Waterfall, transient_mask=False, mask_time_location: int =0):
        """
        Args:
        -----
        indices (list[int]): list of indices which are marked as RFI
        mask_dictionary (list(dict{range(int,int):int})): list of dictionaries where the keys are
        the range of which to use the mask contained in the values. Each dictionary corresponds to
        a different time aggregate.
        time_axis_split (list(range(int))): list of the ranges for which the data is split along the
        time axis when mark_RFI is called.
        transient_mask (bool): If True then blimpy will construct a seperate mask for each marked index
        depending on it's location along the time axis and the value of aggregate in the mark_RFI function.
        This will produce a better fitting mask, but on marginally compared to the increase in computation time.
        mask_time_location (int): If transient_mask is False then this umber will be used to decide where the mask
        shall be calculated from. Eg. =32 will use the mask calculated for t=32 as the mask for all time.
        Can't have t=len(time_axis)
        """

        self.waterfall = waterfall_object
        self.indices = np.array([], dtype=int)
        self.masks = []
        self.time_axis_split = []
        self.transient_mask = transient_mask
        self.mask_time_location = mask_time_location


    def get_spikes(self, data_length, chans_per_coarse_chan_num):
        """Finds the indices of the DC spikes of a GBT spectrum based on a priori
        knowledge of the periodicity of the power spectrum, assuming spikes occur
        in the center of each coarse channel. e.g., if the period is 1024 points,
        then the half-period is 512, and spikes occur at indices i={511, 511+1024,
        511+2*(1024), ...}.

        Args:
        -----
        data_length (int): Length of data in frequency axis.
        chans_per_coarse_chan_num (int): number of channels per coarse channel
        selection size. eg. for a 0002 file with 5 coarse chans selected
        chans_per_coarse_chan_num = 1024*5.

        Returns:
        --------
        spikes (np.array): List of indices at which the spikes occur.
        """

        #last index calculated to be the highest number of coarse channels that fit into data.
        last_index = data_length+1
        spikes = np.arange(0, last_index, chans_per_coarse_chan_num)

        return spikes.astype(int)


    def grab_data_index(self, array, start_index, end_index, direction):
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


    def mark_RFI(self, aggregate: int =1, central_lower: float =0.2, central_upper: float =0.8, sigma: float =10, coarse_chan_num: int =1) -> None:
        """Marks potential RFI signals from data by removing data points with a high
        standard deviation above the mean signal to noise ratio.

        Args:
        -----
        aggregate (int): Number to divide data in the time axis. Larger number makes
        RFI search more refined and more likely to find transient RFI. Must be less
        than number of time indices.
        coarse_chan_num (int): Number of coarse channels to use when calculating SNR
        in order to remove RFI.
        central_lower (float): Lower bound of data (as a fraction of total data) to
        keep when calculating the central median of a coarse channel. eg. 0.3
        removes bottom 30% of data before median calculation.
        central_upper (float): Upper bound of data (as a fraction of total data) to
        keep when calculating the central median of a coarse channel. eg. 0.7
        removes top 30% of data before median calculation.
        sigma (float): Remove any data points with a value greater than sigma standard
        deviations away from the central median.
        """

        #check to make sure aggregate is less than number of time indices in the file.
        a = []
        time_indices = len(self.waterfall.data)
        if aggregate > time_indices:
            print("Aggregate greater than number of time indices in file, %i, please choose a smaller or equal number." % (time_indices))
            return a

        if isinstance(aggregate, int) == False:
            print("aggregate is not an interger")
            return a

        #make sure central_lower < central_upper
        if central_lower > central_upper:
            print("Lower_central is greater than upper_central")
            return a

        n_coarse_chan = self.waterfall.calc_n_coarse_chan()
        if coarse_chan_num > n_coarse_chan:
            print("coarse_chan_num greater than number of coarse channels")
            return a

        if isinstance(coarse_chan_num, int) == False:
            print("coarse_chan_num is not an interger")
            return a

        fine_chans_per_coarse_chan = int(self.waterfall.header[b'nchans']/n_coarse_chan)
        chans_per_coarse_chan_num = fine_chans_per_coarse_chan*coarse_chan_num
        data_length = len(self.waterfall.data[0][0])

        #find spike points defining coarse channels, can grab more than one. For 0002 files this is 1024
        spikes = self.get_spikes(data_length, chans_per_coarse_chan_num)

        #aggregate along time axis data to hopefully detect more RFI
        data_split = np.array_split(self.waterfall.data, aggregate, axis=0)

        #Array of indices to be removed. The elements of this array are memebers of the Marked_index class.
        marked_indices = []


        # Keeps track of how far up the time axis we are for calculating self.time_axis_split
        time = 0
        if not self.transient_mask:
            self.time_axis_split.append(range(len(self.waterfall.data)))

        for data in data_split:
            if self.transient_mask:
                self.time_axis_split.append(range(time, time+len(data)))
            time += len(data)
            data = np.median(data, axis=0)[0]
            mask_dictionary = OrderedDict()
            coarse_channels = []
            indices = np.array([], dtype=int)
            #Make array of coarse channels
            for i in range(len(spikes)-1):
                 coarse_chan_data = self.grab_data_index(data, spikes[i], spikes[i+1]-1, self.waterfall.header[b'foff'])
                 coarse_channels.append(db(coarse_chan_data))

            #Append the end data points in the chans_per_coarse_chan_num does not divide into the length of the data
            if data_length%chans_per_coarse_chan_num != 0:
                coarse_channels.append(db(self.grab_data_index(data, spikes[-1], data_length-1, self.waterfall.header[b'foff'])))

            for i, coarse in enumerate(coarse_channels):
                coarse_chan_data_sorted = np.sort(coarse)

                #remove top % and bottom % of array, this cuts RFI and FFT spikes
                coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[int(len(coarse)*central_upper):int(len(coarse)+1)])
                coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[0:int(len(coarse)*central_lower)])


                std = np.std(coarse_chan_data_sorted)
                median = np.median(coarse_chan_data_sorted)

                #return indexes that maybe RFI above a specified standard deviation
                coarse_normalized = (coarse - median)/std
                coarse_marked = np.argwhere(coarse_normalized > sigma)


                #create list of marked indexes
                if len(coarse_marked) != 0:
                    coarse_marked = coarse_marked.flatten()
                    coarse_marked = coarse_marked + spikes[i]


                indices = np.append(indices, coarse_marked)
                if i < len(spikes)-1:
                    range_end = spikes[i+1]
                # catch the last coarse channel.
                else:
                    range_end = data_length-1
                mask_dictionary[range(spikes[i], range_end)] = lin(median)

            self.indices = np.append(self.indices, indices)
            if self.transient_mask or (time-len(data)) <= self.mask_time_location < time:
                self.masks.append(mask_dictionary)


    def frequency_cut(self, fill_mask=True):
        """ Cut entire frequency column for a given index. This must be called after mark_RFI.

        Args:
            fill_mask (bool): If True then masked RFI will be replaced by the
            median value of the coarse channel without RFI or by whatever value is
            assigned to the index to be cut by the add_indices and add_freq_range
            functions.

        Returns:
            (file): instance of Waterfall with cut data
        """

        #make mask
        if not fill_mask:
            self.indices.flatten()
            self.indices = np.unique(self.indices)
            self.waterfall.whichmask = True
            mask = np.zeros(self.waterfall.data.shape, dtype=bool)
            for index in self.indices:
                for i in range(len(self.waterfall.data)):
                    mask[i][0][index] = True
            self.waterfall.data = np.ma.array(self.waterfall.data, mask=mask)
        else:
            #Remove duplicate indices if not wanting transient mask
            if not self.transient_mask:
                self.indices.flatten()
                self.indices = [np.unique(self.indices)]
            # Take advantage of the fact that the data is sorted
            for i, range in enumerate(self.time_axis_split):
                indices = self.indices[i]
                for coarse_channel, mask in self.masks[i].items():
                    while len(indices) != 0:
                        index = indices[0]
                        if index in coarse_channel:
                            for j in range:
                                    self.waterfall.data[j][0][index] = mask
                            indices = np.delete(indices, 0)
                        else:
                            break
                    if len(indices) == 0:
                        break

        return self.waterfall


    # def add_indices(marked_indices, *frequencies):
    #     """ Add more indices to be cut to the marked_indices array.
    #     TODO: Needs updating to include Marked_index class
    #
    #     Args:
    #         marked_indices (list of ints): Array of marked indices .
    #         *frequencies (ints or list of ints): Frequency indices to be cut. Takes list or single values
    #
    #     Returns:
    #         (array of ints): Marked indices with new frequencies added
    #     """
    #
    #     for freq in frequencies:
    #         flat_list = []
    #         if type(freq) is list:
    #             for item in freq:
    #                  marked_indices.append(item)
    #         else:
    #             marked_indices.append(freq)
    #
    #     marked_indices = np.unique(marked_indices)
    #
    #     return marked_indices


    def frequency_to_index(self, frequency, f_increment, f0):
        return int(np.round((frequency - f0) / f_increment))


    def add_freq_range(self, frequency_lower, frequency_upper, mask=True):
        """ Add an entire frequency range of indices to marked_indices.

        Args:
            frequency_lower (float): Lower bound of frequency range to be cut (MHz).
            frequency_upper (float): Upper bound of frequency range to be cut (MHz).
            mask (bool or float): Determines mask for frequency range. If set to
            a bool value then the mask is calculated based on the the coarse
            channels enclosing the frequency range. Can give own value in units of dB

        Returns:
            (array of ints): Marked indices with new frequencies added

        TODO: Make this work for the transient_mask option as well
        """

        f_increment = self.waterfall.file_header[b'foff']
        f0 = self.waterfall.container.f_stop
        if f_increment < 0:
            a = frequency_lower
            frequency_lower = frequency_upper
            frequency_upper = a

        #Find indices corresponding to frequency range.
        frequency_lower = self.frequency_to_index(frequency_lower, f_increment, f0)
        frequency_upper = self.frequency_to_index(frequency_upper, f_increment, f0)
        frequency_range = np.arange(frequency_lower, frequency_upper)

        #Create mask based on enclosing coarse channels
        if mask == True or mask == False:
            data = np.median(self.waterfall.data, axis=0)[0]

            n_coarse_chan = self.waterfall.calc_n_coarse_chan()
            chans_per_coarse = int(self.waterfall.header[b'nchans']/n_coarse_chan)
            spikes = self.get_spikes(len(self.waterfall.data[0][0]), chans_per_coarse)
            spike_lower = spikes[bisect.bisect_left(spikes, (frequency_lower+0.1))]
            spike_upper = spikes[bisect.bisect_right(spikes, (frequency_upper-0.1))]

            spike_lower = self.grab_data_index(data, spike_lower - chans_per_coarse, spike_lower, f_increment)
            spike_upper = self.grab_data_index(data, spike_upper, spike_upper + chans_per_coarse, f_increment)

            median_lower = np.median(spike_lower)
            median_upper = np.median(spike_upper)

            #Calculate average of median values to make mask
            mask = (median_lower + median_upper)/2
        #want all masks in units of counts but given in dB so need to convert
        else:
            mask = lin(mask)

        #Create cut
        self.indices = np.append(self.indices, frequency_range)

        #Add mask to mask dictionary
