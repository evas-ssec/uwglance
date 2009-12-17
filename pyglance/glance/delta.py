#!/usr/bin/env python
# encoding: utf-8
"""
Routines to do assorted difference and comparison calculations and statistics

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging
import numpy as np
from numpy import *
from scipy.stats import pearsonr, spearmanr, pointbiserialr

compute_r = pearsonr #spearmanr

LOG = logging.getLogger(__name__)


def _missing(x,missing_value=None):
    if missing_value is not None:
        return isnan(x) | (x==missing_value)
    return isnan(x)

def diff(aData, bData, epsilon=0.,
         (a_missing_value, b_missing_value)=(None, None),
         (ignore_mask_a, ignore_mask_b)=(None, None)):
    """
    take two arrays of similar size and composition
    if an ignoreMask is passed in values in the mask will not be analysed to
    form the various return masks and the corresponding spots in the
    "difference" return data array will contain fill values (selected
    based on data type).
    
    return difference array filled with fill data where differences aren't valid,
    good mask where values are finite in both a and b
    trouble mask where missing values or nans don't match or delta > epsilon
    (a-notfinite-mask, b-notfinite-mask)
    (a-missing-mask, b-missing-mask)
    """
    shape = aData.shape
    assert(bData.shape==shape)
    assert(can_cast(aData.dtype, bData.dtype) or can_cast(bData.dtype, aData.dtype))
    
    # if the ignore masks do not exist, set them to include none of the data
    if (ignore_mask_a is None) :
        ignore_mask_a = zeros(shape,dtype=bool)
    if (ignore_mask_b is None) :
        ignore_mask_b = zeros(shape,dtype=bool)
    
    # deal with the basic masks
    a_not_finite_mask, b_not_finite_mask = ~isfinite(aData) & ~ignore_mask_a, ~isfinite(bData) & ~ignore_mask_b
    a_missing_mask, b_missing_mask = zeros(shape,dtype=bool), zeros(shape,dtype=bool)
    # if we were given missing values, mark where they are in the data
    if a_missing_value is not None:
        a_missing_mask[aData == a_missing_value] = True
        a_missing_mask[ignore_mask_a] = False # don't analyse the ignored values
    if b_missing_value is not None:
        b_missing_mask[bData == b_missing_value] = True
        b_missing_mask[ignore_mask_b] = False # don't analyse the ignored values
    
    # build the comparison data that includes the "good" values
    valid_in_a_mask = ~(a_not_finite_mask | a_missing_mask | ignore_mask_a)
    valid_in_b_mask = ~(b_not_finite_mask | b_missing_mask | ignore_mask_b)
    valid_in_both = valid_in_a_mask & valid_in_b_mask
    
    # figure out our shared data type
    sharedType = aData.dtype
    if (aData.dtype is not bData.dtype) :
        sharedType = common_type(aData, bData)
    LOG.debug('Shared data type that will be used for diff comparison: ' + str(sharedType))
    
    # construct our diff'ed array
    raw_diff = zeros(shape, dtype=sharedType) #empty_like(aData)
    
    fill_data_value = select_fill_data(sharedType)
    
    LOG.debug('current fill data value: ' + str(fill_data_value))
    
    raw_diff[~valid_in_both] = fill_data_value # throw away invalid data
    raw_diff[valid_in_both] = bData[valid_in_both] - aData[valid_in_both]
    
    # the valid data which is too different between the two sets according to the given epsilon
    outside_epsilon_mask = (abs(raw_diff) > epsilon) & valid_in_both
    # trouble points = mismatched nans, mismatched missing-values, differences that are too large 
    trouble_pt_mask = (a_not_finite_mask ^ b_not_finite_mask) | (a_missing_mask ^ b_missing_mask) | outside_epsilon_mask
    
    return raw_diff, valid_in_both, (valid_in_a_mask, valid_in_b_mask), trouble_pt_mask, outside_epsilon_mask,  \
           (a_not_finite_mask, b_not_finite_mask), (a_missing_mask, b_missing_mask), (ignore_mask_a, ignore_mask_b)

def select_fill_data(dTypeValue) :
    """
    select a fill data value based on the type of data that is being
    inspected/changed
    """
    
    fill_value_to_return = None
    
    if issubdtype(dTypeValue, np.float) or issubdtype(dTypeValue, np.complex) :
        fill_value_to_return = nan
    elif issubdtype(dTypeValue, np.int) :
        fill_value_to_return = np.iinfo(dTypeValue).min
    elif issubdtype(dTypeValue, np.bool) :
        fill_value_to_return = True
    elif ((dTypeValue is np.uint8)  or
          (dTypeValue is np.uint16) or
          (dTypeValue is np.uint32) or
          (dTypeValue is np.uint64)) :
        fill_value_to_return = np.iinfo(dTypeValue).max
    
    return fill_value_to_return

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

# TODO, should this ultimately be removed?
def rms_corr_withnoise(truth, actual, noiz, epsilon=0., (amissing,bmissing)=(None,None), plot=None):
    """ compute RMS and R statistics for truth vs actual and truth+noiz vs actual
    """
    x=truth
    y=actual
    z=noiz
    d, good, _, bad, _, _, _, _ = diff(x,y,epsilon,(amissing,bmissing))
    #d,good,bad,_,_,_ = diff(x,y,epsilon,(amissing,bmissing))
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
    dpn, good, _, bad, _, _, _, _ = diff(xpn,y,epsilon,(amissing,bmissing))
    #dpn,good,bad,_,_,_ = diff(xpn,y,epsilon,(amissing,bmissing))
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
    
    # if there are no values after the mask,
    # we can't do any of these forms of
    # statistical analysis
    if sum(mask) <= 0 :
        return { }
    
    absDiffData = abs(diffData)
    rms = calculate_root_mean_square(diffData, mask)
    return {    'rms_diff': rms, 
                'std_diff': std(absDiffData[mask]), 
                'mean_diff': mean(absDiffData[mask]), 
                'median_diff': median(absDiffData[mask]),
                'max_diff': max(absDiffData[mask])
                }

def convert_mag_dir_to_U_V_vector(magnitude_data, direction_data, invalidMask=None):
    """
    This method is intended to convert magnitude and direction data into (U, V) vector data.
    An invalid mask may be given if some of the points in the set should be masked out.
    
    TODO, this method is not fully tested
    """
    
    if invalidMask is None :
        invalidMask = zeros(magnitude_data.shape, dtype=bool)
    
    uData = zeros(magnitude_data.shape, dtype=float)
    uData[invalidMask]  = nan
    uData[~invalidMask] = magnitude_data[~invalidMask] * cos (direction_data[~invalidMask])
    
    vData = zeros(magnitude_data.shape, dtype=float)
    vData[invalidMask]  = nan
    vData[~invalidMask] = magnitude_data[~invalidMask] * sin (direction_data[~invalidMask])
    
    return uData, vData

def calculate_root_mean_square (data, goodMask=None) :
    """
    calculate the root mean square of the data,
    possibly selecting only the points in the given
    goodMask, if no mask is given, all points will
    be used
    """
    if goodMask is None:
        goodMask = np.ones(data.shape, dtype=bool)
    
    rootMeanSquare = sqrt( sum( data[goodMask] ** 2 ) / sum( goodMask ) )
    
    return rootMeanSquare

def organize_dimensions_and_produce_rms_info(dataA, dataB, binDimensionIndexNum, tupleDimensionIndexNum,
                                             epsilon=0., (a_missing_value, b_missing_value)=(None, None),
                                             (ignore_mask_a, ignore_mask_b)=(None, None)) :
    """
    reorganize the dimensions in the given data sets in order to put the bin dimension first
    and the tuple dimension last
    collapse all other dimensions into a single "case" dimension, but preserve information on their
    previous shape so that the coresponding index of a given case can be recovered
    
    once the data is reshaped to the desired ordering with the appopriate set of cases, diff the
    data and calculate the rms on a per case basis
    
    return reshaped original and diffrence data sets, masks as corresponding to the output of diff and
    matching the shape of the reshaped data, and information on the shape of the collapsed case dimensions
    
    TODO, this method has not yet been tested or integrated with the rest of glance
    """
    
    # make sure we meet the minimal requirement for compatable shapes/index sizes
    assert(dataA.shape is dataB.shape)
    assert(binDimensionIndexNum   < len(dataA.shape))
    assert(tupleDimensionIndexNum < len(dataA.shape))
    assert(tupleDimensionIndexNum is not binDimensionIndexNum)
    
    # figure out the order to put the index we want first
    new_index_order = range(len(dataA.shape))
    smaller = tupleDimensionIndexNum
    larger  = binDimensionIndexNum
    if (tupleDimensionIndexNum > binDimensionIndexNum) :
        larger  = tupleDimensionIndexNum
        smaller = binDimensionIndexNum
    del(new_index_order[larger])
    del(new_index_order[smaller])
    new_index_order = [binDimensionIndexNum] + new_index_order + [tupleDimensionIndexNum]
    
    # copy our data, then put the indexes into the new order
    new_a_data = dataA.copy()
    new_a_data = new_a_data.transpose(new_index_order)
    new_b_data = dataB.copy()
    new_b_data = new_b_data.transpose(new_index_order)
    # if our masks need to be reordered, do that
    new_a_missing_mask = a_missing_mask
    if new_a_missing_mask is not None :
        new_a_missing_mask = new_a_missing_mask.copy()
        new_a_missing_mask = new_a_missing_mask.transpose(new_index_order)
    new_b_missing_mask = b_missing_mask
    if new_b_missing_mask is not None :
        new_b_missing_mask = new_b_missing_mask.copy()
        new_b_missing_mask = new_b_missing_mask.transpose(new_index_order)
    new_a_ignore_mask = ignore_mask_a
    if new_a_ignore_mask is not None :
        new_a_ignore_mask = new_a_ignore_mask.copy()
        new_a_ignore_mask = new_a_ignore_mask.transpose(new_index_order)
    new_b_ignore_mask = ignore_mask_b
    if new_b_ignore_mask is not None :
        new_b_ignore_mask = new_b_ignore_mask.copy()
        new_b_ignore_mask = new_b_ignore_mask.transpose(new_index_order)
    
    # get the shape information on the internal dimensions we're going to combine
    case_dimension_original_shape = new_b_data.shape[1:-1]
    
    # combine the internal dimensions, to figure out what shape things
    # will be with the flattened cases
    numDimensionsToFlatten = len(case_dimension_original_shape)
    sizeAfterFlattened = np.multiply.accumulate(case_dimension_original_shape)[-1]
    newShape = (new_a_data.shape[0], sizeAfterFlattened, new_a_data.shape[-1])
    
    # flatten the case dimensions
    new_a_data = new_a_data.reshape(newShape)
    new_b_data = new_b_data.reshape(newShape)
    if new_a_missing_mask is not None :
        new_a_missing_mask = new_a_missing_mask.reshape(newShape)
    if new_b_missing_mask is not None :
        new_b_missing_mask = new_b_missing_mask.reshape(newShape)
    if new_a_ignore_mask is not None :
        new_a_ignore_mask = new_a_ignore_mask.reshape(newShape)
    if new_b_ignore_mask is not None :
        new_b_ignore_mask = new_b_ignore_mask.reshape(newShape)
    
    # diff our data
    rawDiffData, goodInBothMask, (goodInAMask, goodInBMask), troubleMask, outsideEpsilonMask, \
        (aNotFiniteMask, bNotFiniteMask), (aMissingMask, bMissingMask), (finalAIgnoreMask, finalBIgnoreMask) = diffOutput = \
            diff(new_a_data, new_b_data, epsilon, (new_a_missing_mask, new_b_missing_mask), (new_a_ignore_mask, new_b_ignore_mask))
    
    # calculate the rms diffs for all the cases
    caseRMSInfo = np.zeros(rawDiffData.shape[:-1], dtype=np.float64)
    # TODO, how to do this without a loop?
    for caseNum in range(sizeAfterFlattened) :
        caseRMSInfo[caseNum] = calculate_root_mean_square(rawDiffData[caseNum], goodMask=goodInBothMask[caseNum])
    
    # return the reshaped original data,
    # information on the original shape of the case dimension that was flattened,
    # the rms diff info,
    # and then the full set of info that came out of the diff between the reshaped a and b data
    return (new_a_data, new_b_data), case_dimension_original_shape, caseRMSInfo, diffOutput
    """
    diffOutput is in the form:
    
        rawDiffData, goodInBothMask, (goodInAMask, goodInBMask), troubleMask, outsideEpsilonMask, \
        (aNotFiniteMask, bNotFiniteMask), (aMissingMask, bMissingMask), (finalAIgnoreMask, finalBIgnoreMask)
    """

def _get_num_perfect(a, b, ignoreMask=None):
    numPerfect = 0
    if not (ignoreMask is None) :
        numPerfect = sum(a[~ignoreMask] == b[~ignoreMask])
    else :
        numPerfect = sum(a == b)
    return numPerfect

def _get_nan_stats(a_nan_mask, b_nan_mask) :
    """
    Get a list of statistics about non-numerical values in data sets a and b,
    the return value will be a dictionary of statistics
    """
    # find the nan values in the data
    num_a_nans = sum(a_nan_mask)
    num_b_nans = sum(b_nan_mask)
    num_common_nans = sum(a_nan_mask & b_nan_mask)
    
    # make the assumption that a and b are the same size and only use the size of a
    total_num_values = a_nan_mask.size
    
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
    
def _get_finite_data_stats(a_is_finite_mask, b_is_finite_mask, common_ignore_mask) :
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
    num_finite_in_only_one = sum((a_is_finite_mask ^ b_is_finite_mask) & (~ common_ignore_mask)) # use an exclusive OR
    
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

def _get_general_data_stats(a, b,
                            a_missing_value, b_missing_value,
                            epsilon, 
                            spatial_ignore_in_a_mask, spatial_ignore_in_b_mask,
                            bad_in_a, bad_in_b
                            ) :
    """
    Get a list of general statistics about a and b, given a and b and some other information
    about them.
    the return value will be a dictionary of statistics
    """
    # figure out how much spatial trouble we had
    num_ignored_in_a = sum(spatial_ignore_in_a_mask)
    num_ignored_in_b = sum(spatial_ignore_in_b_mask)
    
    # get the number of data points
    total_num_values = a.size
    
    general_stats = {'a_missing_value': a_missing_value,
                     'b_missing_value': b_missing_value,
                     'epsilon': epsilon,
                     'max_a': max_with_mask(a, bad_in_a),
                     'max_b': max_with_mask(b, bad_in_b),
                     'min_a': min_with_mask(a, bad_in_a),
                     'min_b': min_with_mask(b, bad_in_b),
                     'num_data_points': total_num_values,
                     'shape': a.shape,
                     'spatially_invalid_pts_ignored_in_a': num_ignored_in_a,
                     'spatially_invalid_pts_ignored_in_b': num_ignored_in_b
                     }
    
    return general_stats

def _get_numerical_data_stats(a, b, diff_data,  data_is_finite_mask,
                              outside_epsilon_mask, trouble_mask,
                              additional_statistics={}) : 
    """
    Get a list of numerical comparison related statistics about a and b,
    given a and b and some other information about them.
    the return value will be a dictionary of statistics
    """
    # calculate our various statistics
    num_finite_values_too_different = sum(outside_epsilon_mask)
    num_perfect = _get_num_perfect(a, b, ~data_is_finite_mask)
    r_corr = corr(a, b, data_is_finite_mask)
    num_trouble = sum(trouble_mask)
    
    # we actually want the total number of _finite_ values rather than all the data
    total_num_finite_values = sum(data_is_finite_mask)
    
    # no dividing by 0!
    fraction_too_different = 0.0
    fraction_perfect = 0.0
    if total_num_finite_values > 0 :
        fraction_too_different = num_finite_values_too_different / float(total_num_finite_values)
        fraction_perfect = num_perfect / float(total_num_finite_values)
    
    comparison = {  'correlation': r_corr,
                    'diff_outside_epsilon_count': num_finite_values_too_different,
                    'diff_outside_epsilon_fraction': fraction_too_different,
                    'perfect_match_count': num_perfect,
                    'perfect_match_fraction': fraction_perfect,
                     'trouble_points_count': num_trouble, 
                     'trouble_points_fraction': float(num_trouble) / float(a.size)
                    }
    comparison.update(additional_statistics)
    
    return comparison

# get the min, ignoring the stuff in mask
def min_with_mask(data, mask) :
    temp = data[~mask]
    toReturn = None
    if len(temp) > 0 :
        toReturn = temp[temp.argmin()]
    return toReturn

# get the max, ignoring the stuff in mask
def max_with_mask(data, mask) :
    temp = data[~mask]
    toReturn = None
    if len(temp) > 0 :
        toReturn = temp[temp.argmax()]
    return toReturn

def summarize(a, b, epsilon=0., (a_missing_value, b_missing_value)=(None,None), ignoreInAMask=None, ignoreInBMask=None):
    """return dictionary of statistics dictionaries
    stats not including 'nan' in name exclude nans in either arrays
    """
    
    diffData, finite_mask, (finite_a_mask, finite_b_mask), \
    trouble, outside_epsilon, (anfin, bnfin), \
    (amis, bmis), (ignoreInAMask, ignoreInBMask) = nfo = diff(a, b, epsilon,
                                                              (a_missing_value, b_missing_value),
                                                              (ignoreInAMask, ignoreInBMask))
    '''
    d, valid_mask, trouble, (anfin, bnfin), (amis, bmis), outside_epsilon = nfo = diff(a,b,
                                                                                       epsilon,
                                                                                       (a_missing_value, b_missing_value),
                                                                                       (ignoreInAMask, ignoreInBMask))
                                                                                       '''
    
    general_stats = _get_general_data_stats(a, b, a_missing_value, b_missing_value, epsilon, 
                                            ignoreInAMask, ignoreInBMask, ~finite_a_mask, ~finite_b_mask) 
    additional_statistics = stats(*nfo) # grab some additional comparison statistics 
    comparison_stats = _get_numerical_data_stats(a, b, diffData, finite_mask, outside_epsilon, trouble, additional_statistics) 
    nan_stats = _get_nan_stats(anfin, bnfin)
    missing_stats = _get_missing_value_stats(amis, bmis)
    finite_stats = _get_finite_data_stats(finite_a_mask, finite_b_mask, (ignoreInAMask | ignoreInBMask)) 
    
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
                    'max_a': 'the maximum finite, non-missing value found in A',
                    'max_b': 'the maximum finite, non-missing value found in B',
                    'min_a': 'the minimum finite, non-missing value found in A',
                    'min_b': 'the minimum finite, non-missing value found in B',
                    'num_data_points': "number of data values in A",
                    'shape': "shape of A",
                    'spatially_invalid_pts_ignored_in_a': 'number of points with invalid latitude/longitude information in A that were' +
                                                            ' ignored for the purposes of data analysis and presentation',
                    'spatially_invalid_pts_ignored_in_b': 'number of points with invalid latitude/longitude information in B that were' +
                                                            ' ignored for the purposes of data analysis and presentation',
                    
                    # finite data stats descriptions
                    'a_finite_count': "number of finite values in A",
                    'a_finite_fraction': "fraction of finite values in A (out of all data points in A)",
                    'b_finite_count': "number of finite values in B",
                    'b_finite_fraction': "fraction of finite values in B (out of all data points in B)",
                    'common_finite_count': "number of finite values in common between A and B",
                    'common_finite_fraction': "fraction of finite values in common between A and B",
                    'finite_in_only_one_count': "number of values that changed finite-ness between A and B; " +
                                                "only the common spatially valid area is considerd for this statistic",
                    'finite_in_only_one_fraction': "fraction of values that changed finite-ness between A and B; " +
                                                "only the common spatially valid area is considerd for this statistic",
                    
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
                    'trouble_points_count': 'number of points that differ in finite/missing status between the input data sets A and B,' +
                                            ' or are unacceptable when compared according to the current epsilon value',
                    'trouble_points_fraction': 'fraction of points that differ in finite/missing status between the input data sets A and B,' +
                                            ' or are unacceptable when compared according to the current epsilon value',
                    
                    # note: the statistics described below may no longer be generated?
                    'mean_percent_change': "Percent change from A to B for finite values, averaged",
                    'max_percent_change': "Percent change from A to B for finite values, maximum value"
                    
                    }
STATISTICS_DOC_STR = '\n'.join( '%s:\n    %s' % x for x in sorted(list(STATISTICS_DOC.items())) ) + '\n'

if __name__=='__main__':
    import doctest
    doctest.testmod()
