import numpy as np
import blimpy as bp
import pylab as plt
from RFI_cleaner import *
import time


# Managed to reduce calculation time down to 20%!!!!
fil = bp.Waterfall("/home/custal/Documents/Code/Waterfall_data/spliced_blc0001020304050607_guppi_57936_37003_HIP116719_0057.gpuspec.0002.fil") # Create a Waterfall object
n_coarse_chan = fil.calc_n_coarse_chan() # Define the number of coarse channels. May be a bug where if a frequency range is defined when creating the Waterfall object this function does not work correctly, I have not investiagted this fully.
fil.blank_dc(n_coarse_chan) # Call this to remove all DC bins which are artefacts from the FFT done on the data. They appear as spikes at the centre of each coarse channel.
RFI = RFI_cleaner(fil)
RFI.mark_RFI(100, 0.3, 0.7, 10, 15) # This function does the bulk of the computation. The input paramters have been found to work best for this file.
RFI.add_freq_range(2307, 2364) # This line demonstartes the ability to manually remove a specific range of frequencies. In this case we are removing the missing data which the mark_RFI does not detect.
fil = RFI.frequency_cut() # After defining all the indices we want to remove in the data we call this function to mask/erase them. The Waterfall object now works as normal but with the RFI removed


plt.figure()
fil.plot_spectrum(logged=True)
plt.figure()
fil.plot_waterfall(logged=True)
plt.show()
