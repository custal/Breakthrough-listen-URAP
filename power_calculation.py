"""
#power_calculation.py

Calculates the total intensity over an entire observation of a filterbank file
in order to compare to on off observations to hopefully detect flare stars which
may have greater power output.

TODO: remove RFI before calculation.

"""

import blimpy as bp
import pylab as plt
import numpy as np
from blimpy.utils import db, lin, rebin, closest, unpack_2to8


def power_integration(*fil_array):
    """ Calculate power integration across an observation.

    Args:

        fil (waterfall): Watefall objects to be integrated over in order to find power and compared

    Returns:

        power (int): Total power of observation
    """
    for fil in fil_array:
        power = np.sum(fil.data)
        print(power)

spliced_blc0001020304050607_guppi_57660_64028_HIP57601_0040.gpuspec.0002.h5
spliced_blc0001020304050607_guppi_57660_64379_HIP58684_0041.gpuspec.0002.h5
spliced_blc0001020304050607_guppi_57660_64723_HIP57780_0042.gpuspec.0002.h5
spliced_blc0001020304050607_guppi_57660_65067_HIP58684_0043.gpuspec.0002.h5
spliced_blc0001020304050607_guppi_57660_65408_HIP57845_0044.gpuspec.0002.h5

fil0 = bp.Waterfall("spliced_blc0001020304050607_guppi_57660_63677_HIP58684_0039.gpuspec.0002.h5")
fil1 = bp.Waterfall("spliced_blc0001020304050607_guppi_57660_64028_HIP57601_0040.gpuspec.0002.h5")
fil2 = bp.Waterfall("spliced_blc0001020304050607_guppi_57660_64379_HIP58684_0041.gpuspec.0002.h5")
fil3 = bp.Waterfall("spliced_blc0001020304050607_guppi_57660_64723_HIP57780_0042.gpuspec.0002.h5")
fil4 = bp.Waterfall("spliced_blc0001020304050607_guppi_57660_65067_HIP58684_0043.gpuspec.0002.h5")
fil5 = bp.Waterfall("spliced_blc0001020304050607_guppi_57660_65408_HIP57845_0044.gpuspec.0002.h5")


power_integration(fil0, fil1, fil2, fil3, fil4, fil5, fil6)
