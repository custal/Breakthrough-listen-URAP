import spec_fit as sf
import numpy as np
import blimpy as bp
import pylab as plt



fil = bp.Waterfall("spliced_blc0001020304050607_guppi_57936_37003_HIP116719_0057.gpuspec.0002.fil")
plot_f, plot_data = fil.plot_spectrum(logged=True)


RFI_fit = sf.spec_fit_STL(plot_f, plot_data, 5488)
plt.plot(plot_f, RFI_fit)
plt.show()
