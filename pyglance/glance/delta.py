#!/usr/bin/env python
# encoding: utf-8
"""
Routines to do assorted difference and comparison calculations and statistics

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging
from numpy import *
from scipy.stats import pearsonr, spearmanr, pointbiserialr

compute_r = spearmanr

LOG = logging.getLogger(__name__)


def _missing(x,missing_value=None):
    if missing_value is not None:
        return isnan(x) | (x==missing_value)
    return isnan(x)

def diff(a, b, epsilon=0., (amissing,bmissing)=(None,None), ignoreMask=None):
    """
    take two arrays of similar size and composition
    if an ignoreMask is passed in values in the mask will not be analysed to
    form the various return masks and the corresponding spots in the
    "difference" return data array will contain nan values.
    return difference array filled with nans where differences aren't valid,
    good mask where values are finite in both a and b
    trouble mask where missing values or nans don't match or delta > epsilon
    (a-notfinite-mask, b-notfinite-mask)
    (a-missing-mask, b-missing-mask)
    """
    shape = a.shape
    assert(b.shape==shape)
    assert(a.dtype==b.dtype)
    
    # if the ignore mask does not exist, set it to include none of the data
    if (ignoreMask is None) :
        ignoreMask = zeros(shape,dtype=bool)
    
    # deal with the basic masks
    anfin, bnfin = ~isfinite(a) & ~ignoreMask, ~isfinite(b) & ~ignoreMask
    amis, bmis = zeros(shape,dtype=bool), zeros(shape,dtype=bool)
    if amissing is not None:
        amis[a==amissing] = True
        amis[ignoreMask] = False # don't analyse the ignored values
    if bmissing is not None:
        bmis[b==bmissing] = True
        bmis[ignoreMask] = False # don't analyse the ignored values
    
    # build the comparison data that includes the "good" values
    d = empty_like(a)
    mask = ~(anfin | bnfin | amis | bmis | ignoreMask)
    d[~mask] = nan
    d[mask] = b[mask] - a[mask]
    
    # trouble areas - mismatched nans, mismatched missing-values, differences > epsilon
    trouble = (anfin ^ bnfin) | (amis ^ bmis) | (abs(d)>epsilon)
    return d, mask, trouble, (anfin, bnfin), (amis, bmis)

def corr(x,y,mask):
    "compute correlation coefficient"
    gf = mask.flatten()
    xf = x.flatten()[gf]
    yf = y.flatten()[gf]
    assert(sum(~isfinite(yf))==0)
    assert(sum(~isfinite(xf))==0)
    # don't try to build a correlation if
    # masking left us with insufficient data
    # to do so
    if (xf.size < 2) or (yf.size < 2) :
        return nan
    return compute_r(xf,yf)[0]

def rms_corr_withnoise(truth, actual, noiz, epsilon=0., (amissing,bmissing)=(None,None), plot=None):
    """ compute RMS and R statistics for truth vs actual and truth+noiz vs actual
    """
    x=truth
    y=actual
    z=noiz
    d,good,bad,_,_ = diff(x,y,epsilon,(amissing,bmissing))
    # compute RMS error
    rmse = sqrt(sum(d[good]**2)) / d.size
    gf = good.flatten()
    # raise NotImplementedError # FIXME we're getting NaN for R
    xf = x.flatten()[gf]
    yf = y.flatten()[gf]
    assert(sum(~isfinite(yf))==0)
    assert(sum(~isfinite(xf))==0)
    # compute correlation coefficient
    r = compute_r(xf,yf)[0]
    # create xpn, x plus noise
    xpn = array(x)
    xpn[good] += z[good]
    xpnf = xpn.flatten()[gf]
    # compute RMS error versus noise
    dpn,good,bad,_,_ = diff(xpn,y,epsilon,(amissing,bmissing))
    rmsepn = sqrt(sum(dpn[good]**2)) / d.size
    assert(sum(~isfinite(xpnf))==0)
    rpn = compute_r(xpnf,yf)[0]
    if plot: plot(xf,xpnf,yf)
    return { 'rms_error': rmse,
             'correlation': r,
             'rms_error_with_noise': rmsepn,
             'correlation_with_noise': rpn,
             }

def stats(diffData, mask, *etc):
    rms = sum(abs(diffData[mask] ** 2)) / diffData.size    
    return {    'rms_diff': rms, 
                'std_diff': std(diffData[mask]), 
                'mean_diff': mean(diffData[mask]), 
                'median_diff': median(diffData[mask]),
                'max_diff': max(diffData[mask])
                }

def _get_num_perfect(a, b, ignoreMask=None):
    numPerfect = 0
    if not (ignoreMask is None) :
        numPerfect = sum(a[~ignoreMask] == b[~ignoreMask])
    else :
        numPerfect = sum(a == b)
    return numPerfect

def _get_nan_stats(a, b) :
    """
    Get a list of statistics about non-numerical values in data sets a and b,
    the return value will be a dictionary of statistics
    """
    # find the nan values in the data
    a_nans = isnan(a)
    num_a_nans = sum(a_nans)
    b_nans = isnan(b)
    num_b_nans = sum(b_nans)
    num_common_nans = sum(a_nans & b_nans)
    
    # make the assumption that a and b are the same size and only use the size of a
    total_num_values = a.size
    
    nan_stats = {'a_nan_count': num_a_nans,
                 'a_nan_fraction': (float(num_a_nans) / float(total_num_values)),
                 'b_nan_count': num_b_nans,
                 'b_nan_fraction': (float(num_b_nans) / float(total_num_values)),
                 'common_nan_count': num_common_nans,
                 'common_nan_fraction': (float(num_common_nans) / float(total_num_values))
                 }
    
    return nan_stats
    
def _get_missing_value_stats(a_missing_mask, b_missing_mask) :
    """
    Get a list of statistics about missing data values in data sets a and b,
    given masks describing the positions of the missing data (such masks
    can be gotten from the diff function), 
    the return value will be a dictionary of statistics
    """
    # calculate information about the missing values
    num_a_missing = sum(a_missing_mask)
    num_b_missing = sum(b_missing_mask)
    num_common_missing = sum(a_missing_mask & b_missing_mask)
    
    # make the assumption that a and b are the same size and only use the size of a's mask
    total_num_values = a_missing_mask.size
    
    missing_value_stats = {'a_missing_count': num_a_missing,
                           'a_missing_fraction': (float(num_a_missing) / float(total_num_values)),
                           'b_missing_count': num_b_missing,
                           'b_missing_fraction': (float(num_b_missing) / float(total_num_values)),
                           'common_missing_count': num_common_missing,
                           'common_missing_fraction': (float(num_common_missing) / float(total_num_values))
                           }
    
    return missing_value_stats
    
def _get_finite_data_stats(a_is_finite_mask, b_is_finite_mask) :
    """
    Get a list of statistics about finite data values in data sets a and b,
    given masks describing the positions of the finite data in each file (such masks
    can be gotten from the diff function),
    the return value will be a dictionary of statistics
    """
    # calculate information about the finite data
    num_a_finite = sum(a_is_finite_mask)
    num_b_finite = sum(b_is_finite_mask)
    num_common_finite = sum(a_is_finite_mask & b_is_finite_mask)
    num_finite_in_only_one = sum(a_is_finite_mask ^ b_is_finite_mask) # use an exclusive OR
    
    # make the assumption that a and b are the same size and only use the size of a's mask
    total_num_values = a_is_finite_mask.size
    
    finite_value_stats = {'a_finite_count': num_a_finite,
                          'a_finite_fraction': (float(num_a_finite) / float(total_num_values)),
                          'b_finite_count': num_b_finite,
                          'b_finite_fraction': (float(num_b_finite) / float(total_num_values)),
                          'common_finite_count': num_common_finite,
                          'common_finite_fraction': (float(num_common_finite) / float(total_num_values)),
                          'finite_in_only_one_count': num_finite_in_only_one,
                          'finite_in_only_one_fraction': (float(num_finite_in_only_one) / float(total_num_values)),
                          }
    
    return finite_value_stats

def _get_general_data_stats(a_missing_value, b_missing_value, epsilon, trouble_mask, ignore_in_a_mask, ignore_in_b_mask) :
    """
    Get a list of general statistics about a and b, given a and b and some other information
    about them.
    the return value will be a dictionary of statistics
    """
    # figure out how much trouble we had
    num_trouble = sum(trouble_mask)
    num_ignored_in_a = sum(ignore_in_a_mask)
    num_ignored_in_b = sum(ignore_in_b_mask)
    
    # make the assumption that a and b are the same size/shape as their trouble mask
    total_num_values = trouble_mask.size
    
    general_stats = {'a_missing_value': a_missing_value,
                     'b_missing_value': b_missing_value,
                     'epsilon': epsilon,
                     'num_data_points': total_num_values,
                     'shape': trouble_mask.shape,
                     'spatially_invalid_pts_ignored_in_a': num_ignored_in_a,
                     'spatially_invalid_pts_ignored_in_b': num_ignored_in_b,
                     'trouble_points_count': num_trouble,
                     'trouble_points_fraction': num_trouble/ float(total_num_values)
                     }
    
    return general_stats

def _get_numerical_data_stats(a, b, diff_data, data_is_finite_mask, epsilon, additional_statistics={}) :
    """
    Get a list of numerical comparison related statistics about a and b,
    given a and b and some other information about them.
    the return value will be a dictionary of statistics
    """
    # calculate our various statistics
    num_finite_values_too_different = sum(abs(diff_data[data_is_finite_mask]) > epsilon)
    num_perfect = _get_num_perfect(a, b, ~data_is_finite_mask)
    r_corr = corr(a, b, data_is_finite_mask)
    
    # we actually want the total number of _finite_ values rather than all the data
    total_num_finite_values = sum(data_is_finite_mask)
    
    comparison = {  'diff_outside_epsilon_count': num_finite_values_too_different,
                    'diff_outside_epsilon_fraction': num_finite_values_too_different / float(total_num_finite_values),
                    'perfect_match_count': num_perfect,
                    'perfect_match_fraction': num_perfect / float(total_num_finite_values),
                    'correlation': r_corr
                    }
    comparison.update(additional_statistics)
    
    return comparison

def summarize(a, b, epsilon=0., (a_missing_value, b_missing_value)=(None,None), ignoreInAMask=None, ignoreInBMask=None):
    """return dictionary of statistics dictionaries
    stats not including 'nan' in name exclude nans in either arrays
    """
    #print('a type: ' + str(a.dtype))
    #print('b type: ' + str(b.dtype))
    
    # select/build our ignore masks
    # if the user didn't send us any, don't ignore anything
    if (ignoreInAMask is None) :
        ignoreInAMask = zeros(a.shape, dtype=bool)
    if (ignoreInBMask is None) :
        ignoreInBMask = zeros(b.shape, dtype=bool)
    ignoreMask = ignoreInAMask | ignoreInBMask
    
    d, mask, trouble, (anfin, bnfin), (amis, bmis) = nfo = diff(a,b,epsilon,(a_missing_value, b_missing_value),ignoreMask)
    
    # build some other finite data masks that we'll need
    finite_a_mask = ~(anfin | amis)
    finite_b_mask = ~(bnfin | bmis)
    finite_mask = finite_a_mask & finite_b_mask
    if not (ignoreInAMask is None) :
        finite_a_mask = finite_a_mask & (~ ignoreInAMask)
    if not (ignoreInBMask is None) :
        finite_b_mask = finite_b_mask & (~ ignoreInBMask)
    if not (ignoreMask is None) :
        finite_mask = finite_mask & (~ ignoreMask)
    
    general_stats = _get_general_data_stats(a_missing_value, b_missing_value, epsilon, trouble, ignoreInAMask, ignoreInBMask) 
    additional_statistics = stats(*nfo) # grab some additional comparison statistics
    comparison_stats = _get_numerical_data_stats(a, b, d, finite_mask, epsilon, additional_statistics) 
    nan_stats = _get_nan_stats(a, b)
    missing_stats = _get_missing_value_stats(amis, bmis)
    finite_stats = _get_finite_data_stats(finite_a_mask, finite_b_mask) 
    
    out = {}
    out['NaN Statistics'] = nan_stats
    out['Missing Value Statistics'] = missing_stats
    out['Finite Data Statistics'] = finite_stats
    out['Numerical Comparison Statistics'] = comparison_stats
    out['General Statistics'] = general_stats
    
    return out

STATISTICS_DOC = {  'general': "Finite values are non-missing and finite (not NaN or +-Inf); fractions are out of all data, " +
                               "both finite and not, unless otherwise specified",
                    
                    # general statistics
                    'a_missing_value': 'the value that is considered \"missing\" data when it is found in A',
                    'b_missing_value': 'the value that is considered \"missing\" data when it is found in B',
                    'epsilon': 'amount of difference between matching data points in A and B that is considered acceptable',
                    'num_data_points': "number of data values in A",
                    'shape': "shape of A",
                    'spatially_invalid_pts_ignored_in_a': 'number of points with invalid latitude/longitude information in A that were' +
                                                            ' ignored for the purposes for data analysis and presentation',
                    'spatially_invalid_pts_ignored_in_b': 'number of points with invalid latitude/longitude information in B that were' +
                                                            ' ignored for the purposes for data analysis and presentation',
                    'trouble_points_count': 'number of points that are nonfinite or missing in either input data set (A or B),' +
                                            ' or are unacceptable when compared (according to the current epsilon value)',
                    'trouble_points_fraction': 'fraction of points that are nonfinite or missing in either input data set (A or B),' +
                                               ' or are unacceptable when compared (according to the current epsilon value)',
                    
                    # finite data stats descriptions
                    'a_finite_count': "number of finite values in A",
                    'a_finite_fraction': "fraction of finite values in A (out of all data points in A)",
                    'b_finite_count': "number of finite values in B",
                    'b_finite_fraction': "fraction of finite values in B (out of all data points in B)",
                    'common_finite_count': "number of finite values in common between A and B",
                    'common_finite_fraction': "fraction of finite values in common between A and B",
                    'finite_in_only_one_count': "number of values that changed finite-ness between A and B",
                    'finite_in_only_one_fraction': "fraction of values that changed finite-ness between A and B",
                    
                    # missing data value statistics
                    'a_missing_count': "number of values flagged missing in A",
                    'a_missing_fraction': "fraction of values flagged missing in A",
                    'b_missing_count': "number of values flagged missing in B",
                    'b_missing_fraction': "fraction of values flagged missing in B",
                    'common_missing_count': "number of missing values in common between A and B",
                    'common_missing_fraction': "fraction of missing values in common between A and B",
                    
                    # NaN related statistics
                    'a_nan_count': "number of NaNs in A",
                    'a_nan_fraction': "fraction of NaNs in A",
                    'b_nan_count': "number of NaNs in B",
                    'b_nan_fraction': "fraction of NaNs in B",
                    'common_nan_count': "number of NaNs in common between A and B",
                    'common_nan_fraction': "fraction of NaNs in common between A and B",
                    
                    # Numerical comparison statistics
                    'correlation': "Pearson correlation r-coefficient (0.0-1.0) for finite values of A and B",
                    'diff_outside_epsilon_count': "number of finite differences falling outside epsilon",
                    'diff_outside_epsilon_fraction': "fraction of finite differences falling outside epsilon (out of common_finite_count)",
                    'max_diff': "Maximum difference of finite values",
                    'mean_diff': "mean difference of finite values",
                    'median_diff': "median difference of finite values",
                    'perfect_match_count': "number of perfectly matched finite data points between A and B",
                    'perfect_match_fraction': "fraction of finite values perfectly matching between A and B (out of common_finite_count)",
                    'rms_diff': "root mean square (RMS) difference of finite values",
                    'std_diff': "standard deviation of difference of finite values",
                    
                    # note: the statistics described below may no longer be generated?
                    'mean_percent_change': "Percent change from A to B for finite values, averaged",
                    'max_percent_change': "Percent change from A to B for finite values, maximum value"
                    
                    }
STATISTICS_DOC_STR = '\n'.join( '%s:\n    %s' % x for x in sorted(list(STATISTICS_DOC.items())) ) + '\n'

if __name__=='__main__':
    import doctest
    doctest.testmod()
