#!/usr/bin/env python
# encoding: utf-8
"""
Routines to do assorted difference and comparison calculations and statistics

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import logging
import numpy as np
from numpy import * # todo, remove this line

from scipy.stats import pearsonr

LOG = logging.getLogger(__name__)

# -------------- generic data manipulation and analysis --------------

def min_with_mask(data, goodMask=None) :
    """
    get the minimum of the values in data,
    inspecting only the values selected in the goodMask
    
    if no mask is passed and a masked array is given the
    mask associated with the array will be used
    if no mask is passed and any other kind of data is
    given, all values are assumed to be good
    """
    
    # if the data array is a masked array and we weren't
    # given a mask, assume that we should use the one
    # associated with the masked array
    if (goodMask is None) and (type(data) is numpy.ma) :
        goodMask = data.mask
    
    # select only the good data values
    goodData = data[goodMask]
    
    # if we have any good data, get the minimum
    toReturn = None
    if goodData.size > 0 :
        toReturn = np.min(goodData)
    
    return toReturn

def max_with_mask(data, goodMask=None) :
    """
    get the maximum of the values in data,
    inspecting only the values selected in the goodMask
    
    if no mask is passed and a masked array is given the
    mask associated with the array will be used
    if no mask is passed and any other kind of data is
    given, all values are assumed to be good
    """
    
    # if the data array is a masked array and we weren't
    # given a mask, assume that we should use the one
    # associated with the masked array
    if (goodMask is None) and (type(data) is numpy.ma) :
        goodMask = data.mask
    
    # select only the good data values
    goodData = data[goodMask]
    
    # if we have any good data, get the maximum
    toReturn = None
    if goodData.size > 0 :
        toReturn = np.max(goodData)
    
    return toReturn

def compute_correlation(xData, yData, goodMask, compute_r_function=pearsonr):
    """
    compute the correlation coefficient of two data sets
    given a mask describing good data values in the sets
    """
    
    # make sure our data sets and mask are the same shape
    assert(xData.shape == yData.shape)
    assert(xData.shape == goodMask.shape)
    
    # pull out just the good data
    good_x_data = xData[goodMask]
    good_y_data = yData[goodMask]
    
    # make sure that there is no remaining bad data
    assert(np.all(np.isfinite(good_x_data)))
    assert(np.all(np.isfinite(good_y_data)))
    
    # if we have enough data, try to build the correlation
    toReturn = np.nan
    if (good_x_data.size >= 2) and (good_y_data.size >= 2) :
        toReturn = compute_r_function(good_x_data, good_y_data)[0]
    
    return toReturn

def calculate_root_mean_square (data, goodMask=None) :
    """
    calculate the root mean square of the data,
    possibly selecting only the points in the given
    goodMask, if no mask is given, all points will
    be used
    """
    
    # get a count of how many good data points we have
    numGoodPoints = data.size
    if goodMask is not None:
        numGoodPoints = np.sum(goodMask)
    
    rootMeanSquare = np.sqrt( np.sum( data[goodMask] ** 2 ) / numGoodPoints )
    
    return rootMeanSquare

# TODO, should the name of this function be changed?
def convert_mag_dir_to_U_V_vector(magnitude_data, direction_data, invalidMask=None, offset_degrees=180):
    """
    This method is intended to convert magnitude and direction data into (U, V) vector data.
    An invalid mask may be given if some of the points in the set should be masked out.
    
    TODO, this method is not fully tested
    """
    
    if invalidMask is None :
        invalidMask = np.zeros(magnitude_data.shape, dtype=bool)
    
    new_direction_data = direction_data[:] + offset_degrees
    
    LOG.debug ("direction data: " + str(new_direction_data[~invalidMask]))
    
    uData = np.zeros(magnitude_data.shape, dtype=float)
    uData[invalidMask]  = np.nan
    uData[~invalidMask] = magnitude_data[~invalidMask] * np.sin (deg2rad(new_direction_data[~invalidMask]))
    
    vData = np.zeros(magnitude_data.shape, dtype=float)
    vData[invalidMask]  = np.nan
    vData[~invalidMask] = magnitude_data[~invalidMask] * np.cos (deg2rad(new_direction_data[~invalidMask]))
    
    return uData, vData

# ------------- bin/tuple related functions --------------------

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

if __name__=='__main__':
    import doctest
    doctest.testmod()
