import numpy as np
import blimpy as bp
import pylab as plt
from RFI_search_via_SNR import *

files = [bp.Waterfall("spliced_blc0001020304050607_guppi_57660_63677_HIP58684_0039.gpuspec.0002.h5"),
bp.Waterfall("spliced_blc0001020304050607_guppi_57660_64028_HIP57601_0040.gpuspec.0002.h5"),
bp.Waterfall("spliced_blc0001020304050607_guppi_57660_64379_HIP58684_0041.gpuspec.0002.h5"),
bp.Waterfall("spliced_blc0001020304050607_guppi_57660_64723_HIP57780_0042.gpuspec.0002.h5"),
bp.Waterfall("spliced_blc0001020304050607_guppi_57660_65067_HIP58684_0043.gpuspec.0002.h5"),
bp.Waterfall("spliced_blc0001020304050607_guppi_57660_65408_HIP57845_0044.gpuspec.0002.h5")]

for fil in files:
    marked_indices = mark_RFI(fil, 100, 0.3, 0.7, 10, chans_per_coarse=5488)
    fil = frequency_cut(fil, marked_indices)
    plt.figure()
    fil.plot_waterfall(logged=True)

plt.show()
