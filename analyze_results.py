
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 09 13:01:42 2015

@author: Dylan Brown

The MIT License (MIT)

Copyright (c) 2015

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import os
import struct
from astropy.io import fits
from matplotlib import cm
import glob
import time
import sys
import re


def numerical_sort(val):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(val)
    parts[1::2] = map(int, parts[1::2])
    return parts

def analyze(cmd_args):

    
    plt.close("all")
    print("Analyzing folder " + str(cmd_args[0]))
    date = cmd_args[0]
    
    if not os.path.exists("pngs/" + date):
        os.makedirs("pngs/" + date)

    savedir = "pngs/" + date

    num_files = 0
    for i in sorted(glob.glob("fits/" + date + "/psi*.fit"), key = numerical_sort):
        num_files += 1
        fits_image = fits.open(i, mode = 'readonly')
        image = fits_image[0].data
        fits_image.close()

        plt.ioff
        fig, ax = plt.subplots()
        ax.imshow(image, cmap = cm.afmhot)

        if num_files <= 10:
            plt.savefig(savedir + "/psi00" + str(num_files) + ".png", dpi = 300)
        elif num_files >= 10:
            plt.savefig(savedir + "/psi0" + str(num_files) + ".png", dpi = 300)
        elif num_files >= 100:
            plt.savefig(savedir + "/psi" + str(num_files) + ".png", dpi = 300)

        plt.close('all')

if __name__ == "__main__":
    
    num_args = len(sys.argv)
    cmd_args = sys.argv[1:num_args]
    print(cmd_args)
    analyze(cmd_args)
