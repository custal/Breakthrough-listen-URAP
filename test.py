"""
#test.py

Just used to play with code.

"""

import blimpy as bp
import pylab as plt
import numpy as np
from blimpy.utils import db, lin, rebin, closest, unpack_2to8

# fil = bp.Waterfall("blc05_guppi_57872_20242_DIAG_2MASS_1502+2250_0024.gpuspec.0002.fil")
#
# fil.plot_spectrum(logged=True)
# fil.plot_spectrum(logged=True, t=4150)
#
# plt.show()
data = np.zeros(20)
print(data)
data_split = np.array_split(data, 30)
print(data_split)
