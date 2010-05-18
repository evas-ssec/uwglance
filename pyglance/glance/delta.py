#!/usr/bin/env python
# encoding: utf-8
"""
Routines to do assorted difference and comparison calculations and statistics

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import logging
import numpy as np
from numpy import *
from scipy.stats import pearsonr

compute_r = pearsonr

LOG = logging.getLogger(__name__)

# TODO, where is this being used?
def _missing(x, missing_value=None):
    if missing_value is not None:
        return isnan(x) | (x==missing_value)
    return isnan(x)

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

if __name__=='__main__':
    import doctest
    doctest.testmod()
