import numpy as np
import blimpy as bp
import pylab as plt
from RFI_cleaner import *


#blc00_guppi_57872_20242_DIAG_2MASS_1502+2250_0024.gpuspec.0002.fil
#spliced_blc0001020304050607_guppi_57936_37003_HIP116719_0057.gpuspec.0002.fil Presence of time spike due to RFI filling entire coarse channel
fil = bp.Waterfall("blc20_guppi_57991_48821_3C161_0006.gpuspec.0002.fil")
n_coarse_chan = fil.calc_n_coarse_chan()
fil.blank_dc(n_coarse_chan)
marked_indices = mark_RFI(fil, 56, 0.3, 0.7, 10, coarse_chan_num=1)
#marked_indices = add_freq_range(fil, marked_indices, 6613, 6617)
fil = frequency_cut(fil, marked_indices, fill_mask=True)


plt.figure()
fil.plot_spectrum(logged=True)
plt.figure()
fil.plot_waterfall(logged=True)
plt.show()

# chan_sizes = [1024, 1792, 3584, 5488, 10976, 21952]
# for i in chan_sizes:
#     marked_indices = mark_RFI(fil, 100, 0.3, 0.7, 10, chans_per_coarse=i)
#     fil = frequency_cut(fil, marked_indices, fill_mask=True)
#     plt.figure()
#     fil.plot_waterfall(logged=True)

# plt.show()
