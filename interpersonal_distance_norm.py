#!/usr/bin/env python3

import numpy as np
from scipy.spatial.distance import cdist

# The text asserts that the average distance between 
# random points in a unit circle is 128 / 45 pi.
# This was solved analytically with Mathematica,
# but this can be checked quickly with numpy.

NMC = 50000

# Random on [0, 1], scaled by 2, - 1
#  ==> NMC x 2D uniform on [-1, 1]
mat = np.random.rand(NMC, 2) * 2 - 1

# Preserve only those within the unit circle.
mat = mat[(mat**2).sum(axis = 1) < 1,:]

# How many remain?
Ncirc = mat.shape[0]

# Calculate all distances
dist2 = cdist(mat, mat, metric = 'euclidean')

# Average distance between them.
print((dist2.sum() / Ncirc**2) / (128 / (45 * np.pi)))

