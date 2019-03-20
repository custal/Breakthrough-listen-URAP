import spec_fit as sf
import numpy as np
import blimpy as bp
import pylab as plt



fil = bp.Waterfall("blc05_guppi_57872_20242_DIAG_2MASS_1502+2250_0024.gpuspec.0002.fil")
plot_f, plot_data = fil.plot_spectrum(logged=True)


RFI_fit = sf.spec_fit_STL(plot_f, plot_data, 1024, chunk_size=65536)
plt.plot(plot_f, RFI_fit)
plt.show()
