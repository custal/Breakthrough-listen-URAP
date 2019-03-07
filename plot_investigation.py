import numpy as np

a = np.array([5, 7, 8, 9, 16, 13, 56])
b = np.array([5, 7, 8, 9, 16, 13, 56, np.nan, np.nan, np.nan])

amean = a.mean()
bmean = b.mean()
ananmean = np.nanmean(a)
bnanmean = np.nanmean(b)

print(amean, bmean, ananmean, bnanmean)
