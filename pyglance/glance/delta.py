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
    
    new_direction_data = direction_data[:] + 180
    
    print ("direction data: " + str(new_direction_data[~invalidMask]))
    
    uData = zeros(magnitude_data.shape, dtype=float)
    uData[invalidMask]  = nan
    uData[~invalidMask] = magnitude_data[~invalidMask] * np.sin (deg2rad(new_direction_data[~invalidMask]))
    
    vData = zeros(magnitude_data.shape, dtype=float)
    vData[invalidMask]  = nan
    vData[~invalidMask] = magnitude_data[~invalidMask] * np.cos (deg2rad(new_direction_data[~invalidMask]))
    
    return uData, vData

# a method to make a list of index numbers for reordering a multi-dimensional array
def _make_new_index_list(numberOfIndexes, firstIndexNumber=0, lastIndexNumber=None) :
    """
    the first and last index numbers represent the dimensions you want to be first and last (respectively)
    when the list is reordered; any other indexes will retain their relative ordering
    
    newIndexList = _make_new_index_list(numIndexes, binIndex, tupleIndex)
    """
    
    if lastIndexNumber is None:
        lastIndexNumber = numberOfIndexes - 1
    
    newIndexList = range(numberOfIndexes)
    maxSpecial   = max(firstIndexNumber, lastIndexNumber)
    minSpecial   = min(firstIndexNumber, lastIndexNumber)
    del(newIndexList[maxSpecial])
    del(newIndexList[minSpecial])
    newIndexList = [firstIndexNumber] + newIndexList + [lastIndexNumber]
    
    return newIndexList

def reorder_for_bin_tuple (data, binIndexNumber, tupleIndexNumber) :
    """
    reorder the data given so that the bin index is first, the tuple index is last,
    and any additional dimensions are flattened into a middle "case" index
    
    the reordered data and the shape of flattened case indexes will be returned
    (note if the original data was only 2 dimensional, None will be returned for the
    shape of the flattened case indexes, since there were no other dimensions to flatten)
    """
    
    # put the bin and tuple dimensions in the correct places
    newIndexList = _make_new_index_list(len(data.shape), binIndexNumber, tupleIndexNumber)
    newData = data.transpose(newIndexList)
    
    # get the shape information on the internal dimensions we're going to combine
    caseOriginalShape = newData.shape[1:-1]
    
    # combine the internal dimensions, to figure out what shape things
    # will be with the flattened cases
    sizeAfterFlattened = np.multiply.accumulate(caseOriginalShape)[-1]
    newShape = (newData.shape[0], sizeAfterFlattened, newData.shape[-1])
    
    # flatten the case dimensions
    newData = newData.reshape(newShape)
    
    # TODO, remove once this is tested
    #print ('original data shape: ' + str(data.shape))
    #print ('original case shape: ' + str(caseOriginalShape))
    #print ('new data shape:      ' + str(newData.shape))
    
    return newData, caseOriginalShape

def determine_case_indecies (flatIndex, originalCaseShape) :
    """
    determine the original indexes of the case
    given the flat index number and the original shape
    
    Note: this method is very memory inefficent
    TODO, find a better way of doing this? does numpy guarantee reshaping strategy?
    """
    
    # create a long flat array with the contents being the index number
    numCases = np.multiply.accumulate(originalCaseShape)[-1]
    temp = np.array(range(numCases))
    
    # reshape the flat array back to the original shape
    # then figure out where our index went
    temp = temp.reshape(originalCaseShape)
    positionOfIndex = np.where(temp == flatIndex)
    
    del temp
    
    return positionOfIndex

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
def min_with_mask(data, mask=None) :
    
    if (mask is None) and (type(data) is numpy.ma) :
        mask = ~data.mask
    if mask is None :
        mask = np.zeros(data.shape, dtype=bool)
    
    temp = data[~mask]
    toReturn = None
    if len(temp) > 0 :
        toReturn = temp[temp.argmin()]
    return toReturn

# get the max, ignoring the stuff in mask
def max_with_mask(data, mask=None) :
    
    if (mask is None) and (type(data) is numpy.ma) :
        mask = ~data.mask
    if mask is None :
        mask = np.zeros(data.shape, dtype=bool)
    
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

def create_colocation_mapping_within_epsilon((alongitude, alatitude),
                                             (blongitude, blatitude),
                                             lonlatEpsilon,
                                             invalidAMask=None, invalidBMask=None):
    """
    match points together based on their longitude and latitude values
    to match points must be within lonlatEpsilon degrees in both longitude and latitude
    
    if the longitude and latitude variables contain invalid data the invalidAMask and
    invalidBMask should be passed with the appropriate masking to remove the invalid values
    
    the return will be in the form of two dictionaries of points, one from a and one from b,
    indexed on the index number in the A or B data where they can be found. Each entry will
    consist of a list of:
        [longitudeValue, latitudeValue, indexNumber, [list of matching indexes in the other set]]
    
    Note: the return will include all pairs of points that match,
    this means an individual a or b point may be repeated if it matches
    multiple points within the lonlatEpsilon provided
    
    Warning: This algorithm will fail to find all matching points if the lonlatEpsilon is set to a
    value greater than or equal to 1.0 degrees. This is related to the bin size used for searching
    thoretically the bin size could be corrected to scale with the lonlatEpsilon in the future. TODO
    """
    assert(alongitude.shape == alatitude.shape)
    assert(blongitude.shape == blatitude.shape)
    assert(lonlatEpsilon >= 0.0)
    
    LOG.debug("Preparing to colocate longitude and latitude points (acceptable epsilon: " + str(lonlatEpsilon) + " degrees)")
    LOG.debug("size of A: " + str(alongitude.shape))
    LOG.debug("size of B: " + str(blongitude.shape))
    
    # make blank invalid masks if none were passed in
    if invalidAMask is None :
        invalidAMask = np.zeros(alongitude.shape, dtype=bool)
    if invalidBMask is None :
        invalidBMask = np.zeros(blongitude.shape, dtype=bool)
    
    # make flat versions of our longitude and latitude
    # so that our index correlations will be simple
    flatALatitude  =  alatitude.ravel()
    flatALongitude = alongitude.ravel()
    flatBLatitude  =  blatitude.ravel()
    flatBLongitude = blongitude.ravel()
    
    # find the ranges of the longitude and latitude
    minLatitude  = min(min(flatALatitude),  min(flatBLatitude))
    maxLatitude  = max(max(flatALatitude),  max(flatBLatitude))
    minLongitude = min(min(flatALongitude), min(flatBLongitude))
    maxLongitude = max(max(flatALongitude), max(flatBLongitude))
    
    # make the bins for the data in longitude/latitude space
    aBins = { }
    bBins = { }
    allAPts = { }
    allBPts = { }
    
    # loop to put all the aData in the bins
    for index in range(flatALatitude.size) :
        filingLat = int( flatALatitude[index])
        filingLon = int(flatALongitude[index])
        
        if (filingLat, filingLon) not in aBins :
            aBins[(filingLat, filingLon)] = [ ]
        
        # create the simple list holding that point in the form:
        # the lon/lat values (for ease of comparison), the index number in A, and the list of matches
        aPoint = [flatALatitude[index], flatALongitude[index], index, [ ]]
        
        # put the point in the list and bin
        allAPts[index] = aPoint
        aBins[(filingLat, filingLon)].append(aPoint)
    
    # loop to put all the bData in the bins
    for index in range(flatBLatitude.size) :
        filingLat = int( flatBLatitude[index])
        filingLon = int(flatBLongitude[index])
        
        if (filingLat, filingLon) not in bBins :
            bBins[(filingLat, filingLon)] = [ ]
        
        # create the simple list holding that point in the form:
        # the lon/lat values (for ease of comparison), the index number in A, and the list of matches
        bPoint = [flatBLatitude[index], flatBLongitude[index], index, [ ]]
        
        # put the point in the list and bin
        allBPts[index] = bPoint
        bBins[(filingLat, filingLon)].append(bPoint)
    
    # for debugging purposes
    totalMatches = 0
    
    # look in all the aData bins and select point pairs that match within epsilon
    for binLatitude, binLongitude in aBins.keys() :
        
        # figure out the lat/lon of the 9 bins that are "near" this one
        toSearch = [ ]
        for latValue in range(binLatitude - 1, binLatitude + 1) :
            for lonValue in range(binLongitude -1, binLongitude + 1) :
                toSearch.append((latValue, lonValue))
        
        # for each A pt in this bin
        for aLat, aLon, aIndex, aMatches in aBins[(binLatitude, binLongitude)] :
            
            # look through my nearby B bins and find any
            # "matching" points that fall within epsilon
            for latValue, lonValue in toSearch :
                
                # if there's anything in that B bin
                if ((latValue, lonValue) in bBins) and (bBins[(latValue, lonValue)] is not None) :
                    
                    # for each data point in the B bin, check if it matches our current A point
                    for bLat, bLon, bIndex, bMatches in bBins[(latValue, lonValue)] :
                        
                        if (abs(bLat - aLat) < lonlatEpsilon) and (abs(aLon - bLon) < lonlatEpsilon) :
                            totalMatches = totalMatches + 1
                            
                            # put the point on our matched lists
                            aMatches.append(bIndex)
                            bMatches.append(aIndex)
    
    LOG.debug('Found ' + str(totalMatches) + ' matched pairs.')
    
    return allAPts, allBPts, totalMatches

def create_colocated_lonlat_with_lon_lat_colocation(listOfColocatedALonLat, listOfColocatedBLonLat,
                                                    totalMatches,
                                                    aLongitude, aLatitude,
                                                    bLongitude, bLatitude) :
    """
    given a pre colocated list of A and B lon/lat info from create_colocation_mapping_within_epsilon,
    match up the longitude and latitude and return the colocated sets
    """
    
    # some general statistics
    multipleMatchesInA      = 0
    multipleMatchesInB      = 0
    totalValidMatchedPairs  = 0
    
    # our final data sets
    matchedLongitude   = np.zeros(totalMatches, dtype=aLongitude.dtype) 
    matchedLatitide    = np.zeros(totalMatches, dtype=aLatitude.dtype) 
    
    # we don't know how many unmatched points we may have
    unmatchedALongitude = [ ]
    unmatchedALatitude  = [ ]
    unmatchedBLongitude = [ ]
    unmatchedBLatitude  = [ ]
    
    # go through the A list and collect all the matches
    currentIndex = 0
    # look through all the A points
    for aIndex in sorted(listOfColocatedALonLat.keys()) :
        
        [aLon, aLat, aIndex, aMatches] = listOfColocatedALonLat[aIndex]
        tempMatches = 0
        
        # for each match you found on a given a point
        for matchIndex in sorted(aMatches) :
            
            [bLon, bLat, bIndex, bMatches] = listOfColocatedBLonLat[matchIndex]
            
            # copy the lon/lat info 
            matchedLongitude[currentIndex] = (aLon + bLon) / 2
            matchedLatitide[currentIndex]  = (aLat + bLat) / 2
            
            currentIndex = currentIndex + 1
            tempMatches  = tempMatches  + 1
        
        # update statistics based on the number of matches
        totalValidMatchedPairs = totalValidMatchedPairs + tempMatches
        if tempMatches > 1 :
            multipleMatchesInA = multipleMatchesInA + tempMatches
        elif tempMatches <= 0 :
            unmatchedALatitude.append(aLat)
            unmatchedALongitude.append(aLon)
    
    # gather some additional statistics from the B list
    # go through each b point
    for bIndex in sorted(listOfColocatedBLonLat) :
        
        [bLon, bLat, bIndex, bMatches] = listOfColocatedBLonLat[bIndex]
        tempMatches = len(bMatches)
        
        # update some statistics based on the number of matches
        if tempMatches > 1 :
            multipleMatchesInB = multipleMatchesInB + tempMatches
        elif tempMatches <= 0 :
            unmatchedBLatitude.append(bLat)
            unmatchedBLongitude.append(bLon)
    
    # make the unmatched lists into proper numpy arrays
    unmatchedALatitude  = np.array(unmatchedALatitude,  dtype=aLatitude.dtype)
    unmatchedALongitude = np.array(unmatchedALongitude, dtype=aLongitude.dtype)
    unmatchedBLatitude  = np.array(unmatchedBLatitude,  dtype=bLatitude.dtype)
    unmatchedBLongitude = np.array(unmatchedBLongitude, dtype=bLongitude.dtype)
    
    LOG.debug("Total matched pairs of longitude/latitide: " + str(totalValidMatchedPairs))
    
    return (matchedLongitude,    matchedLatitide, (multipleMatchesInA, multipleMatchesInB)), \
           (unmatchedALongitude, unmatchedALatitude), \
           (unmatchedBLongitude, unmatchedBLatitude)

def create_colocated_data_with_lon_lat_colocation(listOfColocatedALonLat, listOfColocatedBLonLat,
                                                  colocatedLongitude, colocatedLatitude,
                                                  aData, bData,
                                                  missingData=None, altMissingDataInB=None,
                                                  invalidAMask=None, invalidBMask=None) :
    """
    given a pre colocated list of A and B lon/lat info from create_colocation_mapping_within_epsilon,
    match up the valid data in two data sets and return the list of valid data, padded with missing
    values so that it will match the original longitude and latitude
    """
    
    if altMissingDataInB is None :
        altMissingDataInB = missingData
    
    # some general statistics
    multipleMatchesInA      = 0
    multipleMatchesInB      = 0
    totalValidMatchedPairs  = 0
    
    # our final data sets
    matchedAPoints     = np.ones(colocatedLatitude.shape, dtype=aData.dtype) * missingData
    matchedBPoints     = np.ones(colocatedLatitude.shape, dtype=bData.dtype) * altMissingDataInB
    
    # we don't know how many unmatched points we may have
    unmatchedAPoints    = [ ]
    unmatchedBPoints    = [ ]
    unmatchedALongitude = [ ]
    unmatchedALatitude  = [ ]
    unmatchedBLongitude = [ ]
    unmatchedBLatitude  = [ ]
    
    # go through the A list and sort all the valid matches
    currentIndex = 0
    # go through all the a points
    for aIndex in sorted(listOfColocatedALonLat.keys()) :
        
        [aLon, aLat, aIndex, aMatches] = listOfColocatedALonLat[aIndex]
        tempMatches = 0
        isInvalidA = invalidAMask[aIndex]
        
        # for each point that matched to a given a point
        for matchIndex in sorted(aMatches) :
            
            [bLon, bLat, bIndex, bMatches] = listOfColocatedBLonLat[matchIndex]
            isInvalidB = invalidBMask[bIndex]
            
            # if either of our data points is invalid, then the data doesn't match
            if isInvalidA or isInvalidB :
                # fill in missing data in the matches
                matchedAPoints[currentIndex] = missingData
                matchedBPoints[currentIndex] = altMissingDataInB
            else: # we have a valid match!
                tempMatches = tempMatches + 1
                matchedAPoints[currentIndex] = aData[aIndex]
                matchedBPoints[currentIndex] = bData[bIndex]
            
            currentIndex = currentIndex + 1
        
        totalValidMatchedPairs = totalValidMatchedPairs + tempMatches
        if tempMatches > 1 :
            multipleMatchesInA = multipleMatchesInA + tempMatches
        elif (tempMatches <= 0) and (not isInvalidA) :
            unmatchedAPoints.append(aData[aIndex])
            unmatchedALongitude.append(aLon)
            unmatchedALatitude.append(aLat)
    
    # gather some additional statistics from the B list
    # go through all the b points
    for bIndex in sorted(listOfColocatedBLonLat.keys()) :
        
        [bLon, bLat, bIndex, bMatches] = listOfColocatedBLonLat[bIndex]
        tempMatches = 0
        isInvalidB = invalidBMask[bIndex]
        
        # for each point that matched to a given b point
        for matchIndex in sorted(bMatches) :
            
            [aLon, aLat, aIndex, aMatches] = listOfColocatedALonLat[matchIndex]
            isInvalidA = invalidAMask[aIndex]
            
            # if either of our data points is invalid, then the data doesn't match
            if isInvalidB or isInvalidA :
                # we've already built our matched data, so no need to missing it out
                pass
            else: # we have a valid match!
                tempMatches = tempMatches + 1
        
        if tempMatches > 1 :
            multipleMatchesInB = multipleMatchesInB + tempMatches
        elif (tempMatches <= 0) and (not isInvalidB) :
            unmatchedBPoints.append(bData[bIndex])
            unmatchedBLongitude.append(bLon)
            unmatchedBLatitude.append(bLat)
    
    # make the unmatched lists into proper numpy arrays
    unmatchedAPoints    = np.array(unmatchedAPoints,    dtype=aData.dtype)
    unmatchedBPoints    = np.array(unmatchedBPoints,    dtype=bData.dtype)
    unmatchedALongitude = np.array(unmatchedALongitude, dtype=colocatedLongitude.dtype)
    unmatchedALatitude  = np.array(unmatchedALatitude,  dtype=colocatedLatitude.dtype)
    unmatchedBLongitude = np.array(unmatchedBLongitude, dtype=colocatedLongitude.dtype)
    unmatchedBLatitude  = np.array(unmatchedBLatitude,  dtype=colocatedLatitude.dtype)
    
    LOG.debug("Total matched data point pairs found: " + str(totalValidMatchedPairs))
    
    return (matchedAPoints, matchedBPoints, (multipleMatchesInA, multipleMatchesInB)), \
           (unmatchedAPoints, unmatchedALongitude, unmatchedALatitude), \
           (unmatchedBPoints, unmatchedBLongitude, unmatchedBLatitude)

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
