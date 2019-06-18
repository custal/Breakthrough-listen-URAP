import pylab as plt


# def remove_RFI(file, t=0):
#     """Removes potential RFI signals from data by removing data points with a
#     high standard deviation above the mean signal to noise ratio.
#
#     Args:
#     -----
#     file (waterfall): Waterfall object to have data removed from.
#     t (int): time integration of data to take
#
#     Returns:
#     --------
#     waterfall object with RFI removed
#
#     TODO: make std threshold and number of coarse channels variables (and %
#     data cut above and below)
#     TODO: number of coarse channels, 64, is hard coded in, any way to discern
#     this from file header?
#     """
#
#     #find spike points defining coarse channels
#     chans_per_coarse = int(file.file_header[b'nchans'])/64
#     spikes = get_spikes(file.data[t][0], chans_per_coarse)
#
#     n_coarse_chan = file.calc_n_coarse_chan()
#     file.blank_dc(n_coarse_chan)
#
#     coarse_channels = []
#     #Make array of coarse channels (sorts only first integration at the moment)
#     for i in range(len(spikes)-1):
#          coarse_chan_data = grab_data_index(file.data[t][0], spikes[i], spikes[i+1]-1, file.header[b'foff'])
#          coarse_channels.append(db(coarse_chan_data))
#
#     coarse_channels_masked = []
#     for coarse in coarse_channels:
#         coarse_chan_data_sorted = np.sort(coarse)
#
#         #remove top 10% and bottom 20% of array, this cuts RFI and FFT spikes
#         coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[int(chans_per_coarse*0.9):int(chans_per_coarse+1)])
#         coarse_chan_data_sorted = np.delete(coarse_chan_data_sorted, np.s_[0:int(chans_per_coarse*0.2)])
#
#         std = np.std(coarse_chan_data_sorted)
#         median = np.median(coarse_chan_data_sorted)
#
#         #mask values that maybe RFI above a specified standard deviation
#         coarse_masked = np.ma.array(coarse)
#         coarse_masked = np.ma.masked_where((coarse - median)/std > 5 , coarse)
#         coarse_channels_masked.append(coarse_masked)
#
#
#     all_channels = np.ma.concatenate(coarse_channels_masked)
#     file.data[t][0] = all_channels
#
#     freqs = file.populate_freqs()
#     if file.header[b'foff'] < 0:
#         freqs[::-1]
#
#     plt.plot(freqs, all_channels)
#     return
#     # test_channels = []
#     # for i in range(24,25):
#     #     test_channels.append(coarse_channels_masked[i])
#     #
#     # test_channels = np.ma.concatenate(test_channels)
#     # x = range(1024)
#     # plt.plot(x, test_channels)
#     #
#     # plt.subplot(2, 1, 1)
#     # plt.plot(freqs, db(file.data[t][0]))
#     # plt.title('RFI search via SNR')
#     #
#     # plt.subplot(2, 1, 2)
#     # plt.plot(freqs, all_channels)
