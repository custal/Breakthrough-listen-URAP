"""
#frequency_cut.py

function to remove undesired frequnecies from filterbank data.

Use similar to Waterfall()

   fil = frequency_cut()

TODO: Remove frequency to index function and replace with the attribute of waterfall object that Umran showed which I can't remeber the name of.

"""

import blimpy as bp
import pylab as plt
import numpy as np

def frequency_to_index(frequency, f_increment, f0):
    return int(np.round((frequency - f0) / f_increment))

def frequency_cut(filename, f_start=None, f_stop=None, *ft_cut):
    """ Cut data from specified frequency and time ranges.

    Args:
        filename (string): name of file to be cut
        f_start (float): start frequency in MHz
        f_stop (float): stop frequency in MHz
        ft_cut (array of floats): an nx4 array containing frequnecy and time interval to be cut ie. [[f0, f1, t0, t1], ..., [f2, f3, t2, t3]]

    Returns:
        (fil): instance of Waterfall with cut data
    """

    fil = bp.Waterfall(filename, f_start=f_start, f_stop=f_stop)

    #Cut requested frequency and time ranges
    if ft_cut:
        if f_stop == None:
            f_stop = fil.file_header[b'fch1']

        for ft_interval in ft_cut:

            #Flip frequency array intervals if frequency increment is negative
            if fil.file_header[b'foff'] < 0:
                a = ft_interval[0]
                ft_interval[0] = ft_interval[1]
                ft_interval[1] = a

            #If interval limit is not specified set it to upper/lower of data, if frequency limit is specified convert to index
            if ft_interval[0] == None:
                ft_interval[0] = 0
            else:
                ft_interval[0] = frequency_to_index(ft_interval[0], fil.file_header[b'foff'], f_stop)
            if ft_interval[1] == None:
                ft_interval[1] = len(fil.data[0, 0])
            else:
                ft_interval[1] = frequency_to_index(ft_interval[1], fil.file_header[b'foff'], f_stop)
            if ft_interval[2] == None:
                ft_interval[2] = 0
            if ft_interval[3] == None:
                ft_interval[3] = len(fil.data)

            for i in range(ft_interval[2], ft_interval[3]):
                for j in range(ft_interval[0], ft_interval[1]):
                    fil.data[i, 0, j] = 'nan'

    return fil

fil = frequency_cut("Voyager1.single_coarse.fine_res.fil", 8419, 8421, [8419.5, 8419.8, None, None], [None, None, 6, 9])

fil.plot_all()
plt.show()
