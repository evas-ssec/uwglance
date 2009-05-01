#!/usr/bin/env python
# encoding: utf-8
"""
Routines to do assorted difference and comparison calculations and statistics

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging
from numpy import *

LOG = logging.getLogger(__name__)


def _missing(x,missing_value=None):
    if missing_value is not None:
        return isnan(x) | (x==missing_value)
    return isnan(x)

def diff(a, b, epsilon=0., (amissing,bmissing)=(None,None)):
    """
    take two arrays of similar size and composition
    return difference array filled with nans where differences aren't valid,
    good mask where values are finite in both a and b
    trouble mask where missing values or nans don't match or delta > epsilon
    (a-notfinite-mask, b-notfinite-mask)
    (a-missing-mask, b-missing-mask)
    """
    shape = a.shape
    assert(b.shape==shape)
    assert(a.dtype==b.dtype)
    anfin, bnfin = ~isfinite(a), ~isfinite(b)
    amis, bmis = zeros(shape,dtype=bool), zeros(shape,dtype=bool)
    if amissing is not None:
        amis[a==amissing] = True
    if bmissing is not None:
        bmis[b==bmissing] = True
    d = empty_like(a)
    mask = ~(anfin | bnfin | amis | bmis)
    d[~mask] = nan
    d[mask] = b[mask] - a[mask]
    # trouble areas - mismatched nans, mismatched missing-values, differences > epsilon
    trouble = (anfin ^ bnfin) | (amis ^ bmis) | (abs(d)>epsilon)
    return d, mask, trouble, (anfin, bnfin), (amis, bmis)

def stats(dif, mask, bad, *etc):
    rms = sum(abs(dif[mask] ** 2)) / dif.size    
    return {    'rms_diff': rms, 
                'std_diff': std(dif[mask]), 
                'mean_diff': mean(dif[mask]), 
                'median_diff': median(dif[mask]) 
                }


def summarize(a, b, epsilon=0., (amiss,bmiss)=(None,None)):
    """return dictionary of similarity statistics
    stats not including 'nan' in name exclude nans in either arrays
    """
    d, mask, trouble, (anfin, bnfin), (amis, bmis) = nfo = diff(a,b,epsilon,(amiss,bmiss))
    a_xor_b_finite = sum(anfin ^ bnfin)
    stadic = stats(*nfo)
    n_o_e = sum(trouble)
    a_nans = isnan(a)
    b_nans = isnan(b)
    n = a.size
        
    out = { 'a_xor_b_finite_count': a_xor_b_finite,
             'finite_count': sum(mask),
             'a_missing_count': sum(amis),
             'b_missing_count': sum(bmis),
             'outside_epsilon_count': n_o_e,
             'outside_epsilon_fraction': n_o_e / float(n),
             'max_count': n,
             'a_nan_count': sum(a_nans),
             'b_nan_count': sum(b_nans),
             'a_and_b_nan_count': sum(a_nans & b_nans),
             'shape': a.shape
             }
    out.update(stadic)
    return out

STATISTICS_DOC = {  'general': "Finite values are non-missing and finite (not NaN or +-Inf)",
                    'mean_percent_change': "Percent change from A to B for finite values, averaged",
                    'max_percent_change': "Percent change from A to B for finite values, maximum value",
                    'a_xor_b_finite_count': "number of values that changed finite-ness between A and B",
                    'mean_diff': "Mean difference of finite values",
                    'std_diff': "Stdev of difference of finite values",
                    'max_diff': "Maximum difference of finite values",
                    'median_diff': "Median difference of finite values",
                    'rms_diff': "RMS difference of finite values",
                    'finite_count': "number of finite values in common between A and B",
                    'outside_epsilon_count': "number of finite differences falling outside epsilon",
                    'outside_epsilon_fraction': "fraction of values falling outside epsilon (outside_epsilon_count/max_count)",
                    'max_count': "number of values (cumprod(shape))",
                    'a_missing_count': "number of values flagged missing in A",
                    'b_missing_count': "number of values flagged missing in B",
                    'a_nan_count': "number of NaNs in A",
                    'b_nan_count': "number of NaNs in B",
                    'a_and_b_nan_count': "number of NaNs in common between A and B",
                    'shape': "shape of A"
                    }
STATISTICS_DOC_STR = '\n'.join( '%s:\n    %s' % x for x in sorted(list(STATISTICS_DOC.items())) ) + '\n'

if __name__=='__main__':
    import doctest
    doctest.testmod()
