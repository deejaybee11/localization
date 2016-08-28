# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 17:37:19 2016

@author: dbro184
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt

radius = 10
length = 30
width = 5
x0 = 5
y0 = -15

x = np.linspace(-60, 60, 1024)
y = np.linspace(-60, 60, 1024)

rsq = (x[:, np.newaxis] - x0)**2 + (y - y0)**2

dumbell = np.ones((1024, 1024))

dumbell[rsq <= radius**2] = 0

rsq2 = (x[:, np.newaxis] - x0)**2 + (y - (y0 + length + 2*radius))**2

dumbell[rsq2 <= radius**2] = 0

dumbell[np.logical_and(np.logical_and(abs(x[:,np.newaxis] - x0) <= width / 2, ((y[np.newaxis,:] - y0) >= 0)), ((y[np.newaxis,:] - y0) <= (2*radius + length)))] = 0

plt.imshow(dumbell)

