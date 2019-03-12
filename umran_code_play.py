import spec_fit as sf
import numpy as np
import blimpy as bp
import pylab as plt



fil = bp.Waterfall("spliced_blc0001020304050607_guppi_58100_80058_OUMUAMUA_0015.gpuspec.0002.fil")
plot_f, plot_data = fil.plot_spectrum()

RFI_fit = sf.spec_fit_STL(plot_f, plot_data, 1024, chunk_size=65536)

plt.plot(plot_f, RFI_fit)
plt.show()
