"""
#test.py

Just used to play with code.

"""

import blimpy as bp
import pylab as plt
import numpy as np
import bisect
from blimpy.utils import db, lin, rebin, closest, unpack_2to8

# fil = bp.Waterfall("spliced_blc0001020304050607_guppi_57936_37003_HIP116719_0057.gpuspec.0002.fil")
#
# fil.plot_spectrum(logged=True, t=0)
#
# plt.show()
# # data = np.zeros(10)
# data_split = np.array_split(data, 1)
# print(data)
# print(data_split[0])


# fil = bp.Waterfall("spliced_blc0001020304050607_guppi_57936_37003_HIP116719_0057.gpuspec.0002.fil")
# n_coarse_chan = fil.calc_n_coarse_chan()
# fil.blank_dc(n_coarse_chan)
# print(fil.header[b'nchans']/n_coarse_chan)
# spectrum = db(fil.data[0,0,:1024*7])
# plt.plot(spectrum)

# plt.show()

a = [0,5,10,15,20,25,30,35,40,45,50,55,60]
b = bisect.bisect_right(a, 14.9)
print(b)
