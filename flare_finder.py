"""
#flare_finder.py

Attempts to detect a flares by looking for large power spikes in time.

"""

import blimpy as bp
import pylab as plt
import numpy as np
from blimpy.utils import db, lin, rebin, closest, unpack_2to8
