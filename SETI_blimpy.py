import blimpy as bp
import pylab as plt


obs = bp.Waterfall("Voyager1.single_coarse.fine_res.fil", 8419., 8420)
obs.plot_all()
plt.show()
