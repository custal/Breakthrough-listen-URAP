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

fil1 = bp.Waterfall("blc00_guppi_57872_20242_DIAG_2MASS_1502+2250_0024.gpuspec.0002.fil")
fil2 = bp.Waterfall("blc05_guppi_57872_20242_DIAG_2MASS_1502+2250_0024.gpuspec.0002.fil")
power_integration(fil1, fil2)
