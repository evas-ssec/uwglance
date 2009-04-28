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


def statistics(a, b, epsilon=0., (amiss,bmiss)=(None,None)):
    """return dictionary of similarity statistics
    stats not including 'nan' in name exclude nans in either arrays
    """
    shape = a.shape
    aflat = a.flatten()
    bflat = b.flatten()
    n_a_missing = 0
    n_b_missing = 0
    if amiss is not None:
        mvm = (amiss==aflat)
        n_a_missing = mvm.sum()
        aflat[mvm] = nan
    if bmiss is not None:
        mvm = (bmiss==bflat)
        n_b_missing = mvm.sum()
        bflat[mvm] = nan
    dflat = bflat-aflat
    a_xor_b_finite = sum(isfinite(aflat) ^ isfinite(bflat))
    del aflat
    del bflat
    n = len(dflat)
    dflat = array(dflat, float64)
    a_nans = isnan(a.flatten())
    b_nans = isnan(b.flatten())
    fin = isfinite(dflat) 
    
    dflat = dflat[fin]
    perc = dflat/a.flatten()[fin]*100.
    rms = sqrt(sum(dflat*dflat)/float(len(dflat)))  # should n be len(dflat) or n?
    outside_epsilon = abs(dflat)>epsilon
    n_o_e = outside_epsilon.sum()
    
    return { 'mean_percent_change': perc.mean(),
             'max_percent_change': abs(perc).max(),
             'a_xor_b_finite_count': a_xor_b_finite,
             'mean_diff': dflat.mean(),
             'std_diff': std(dflat),
             'max_diff': abs(dflat).max(),
             'rms_diff': rms,
             'finite_count': len(dflat),
             'a_missing_count': n_a_missing,
             'b_missing_count': n_b_missing,
             'outside_epsilon_count': n_o_e,
             'outside_epsilon_fraction': n_o_e / float(n),
             'max_count': n,
             'a_nan_count': a_nans.sum(),
             'b_nan_count': b_nans.sum(),
             'a_and_b_nan_count': (a_nans & b_nans).sum(),
             'shape': shape
             }

STATISTICS_DOC = {  'general': "Finite values are non-missing and finite (not NaN or +-Inf)",
                    'mean_percent_change': "Percent change from A to B for finite values, averaged",
                    'max_percent_change': "Percent change from A to B for finite values, maximum value",
                    'a_xor_b_finite_count': "number of values that changed finite-ness between A and B",
                    'mean_diff': "Mean difference of finite values",
                    'std_diff': "Stdev of difference of finite values",
                    'max_diff': "Maximum difference of finite values",
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
