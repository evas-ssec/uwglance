#!/usr/bin/env python
# encoding: utf-8
"""
Plotting routines for difference values using matplotlib

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging
from pylab import *

LOG = logging.getLogger(__name__)



def spectral_diff_plot( mean_diff, std_diff, max_diff, min_diff, acceptable_diff=None, x=None ):
    """plot spectral difference in current axes, wrapped in std
    >>> x = arange(0,9,0.1)
    >>> y1 = sin(x)
    >>> y2 = sin(x+0.05)
    >>> d = y2-y1
    >>> s = std(d)
    >>> spectral_diff_plot( d, s, array([0.1] * len(d)), x=x )
    """
    cla()
    if x is None: x = range(len(mean_diff))
    if acceptable_diff is not None:
        plot(x, acceptable_diff, 'g.', hold=True, alpha=0.2)
        plot(x, -acceptable_diff, 'g.', hold=True, alpha=0.2)
    plot(x, mean_diff+std_diff, 'r', hold=True, alpha=0.5)
    plot(x,min_diff,'b')
    plot(x,max_diff,'b')
    plot(x, mean_diff-std_diff, 'r', hold=True, alpha=0.5)
    plot(x, mean_diff, 'k', hold=True)

def compare_spectra(actual, desired=None, acceptable=None, x=None):
    """ given an actual[spectrum][wnum], desired[spectrum][wnum], plot comparisons of differences
    """
    delta = actual-desired if (desired is not None) else actual
    d_mean = mean(delta,axis=0)
    d_max = delta.max(axis=0)
    d_min = delta.min(axis=0)
    d_std = std(delta,axis=0)
    des_mean = mean(desired,axis=0)
    subplot(211)
    cla()
    if x is None: x = range(len(des_mean))
    # plot(x,des_mean+d_max,'b')
    # plot(x,des_mean+d_min,'b')
    plot(x,des_mean,'k')
    grid()
    title("mean spectrum")
    subplot(212)
    spectral_diff_plot(d_mean, d_std, d_max, d_min, acceptable, x)
    grid()
    title("difference min-max (blue), mean (black), mean +/- std (red)")
    
    


if __name__=='__main__':
    import doctest
    doctest.testmod()
