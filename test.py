"""
#test.py

Just used to play with code.

"""

import blimpy as bp
import pylab as plt
import numpy as np
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
def add_freqs(marked_indices, *frequencies):
    """ Add more frequencies to be cut to the marked_indices array.

    Args:
        marked_indices (list of ints): Array of marked indices .
        *frequencies: Frequency indices to be cut. Takes list or single values

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
a = range(5,6)
print(type(a))
