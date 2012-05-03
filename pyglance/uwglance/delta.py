#!/usr/bin/env python
# encoding: utf-8
"""
Routines to do assorted difference and comparison calculations and statistics

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import logging
import math
import numpy as np
from numpy import * # todo, remove this line

from scipy.stats import pearsonr

LOG = logging.getLogger(__name__)

# for calculating the great circle distance, this is the radius well will
# assume that the spherical model of the earth has, in km
SPHERICAL_EARTH_RADIUS = 6373.0

# -------------- generic data manipulation and analysis --------------

# TODO have someone double check the math here
def great_circle_distance (latitudeA, longitudeA, latitudeB, longitudeB) :
    """
    Calculate the great circle distance (in km) between the A and B points
    given in the input parameters, the inputs are expected to be in degrees
    
    note: This method uses the spherical law of cosines, and is best suited
    for smaller distances.
    """
    
    # convert to radians
    latARad = math.radians(latitudeA)
    lonARad = math.radians(longitudeA)
    latBRad = math.radians(latitudeB)
    lonBRad = math.radians(longitudeB)
    
    distToReturn = math.acos(math.sin(latARad) * math.sin(latBRad) +
                             math.cos(latARad) * math.cos(latBRad) *
                             math.cos(lonBRad - lonARad)) * SPHERICAL_EARTH_RADIUS
    
    return distToReturn

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

class BinTupleMapping (object) :
    """
    This class represents a bin / tuple data remapping.
    It encapsulates information about the dimensions that are considered the
    bin and tuple dimensions and is able to transform data into the
    [bin][case][tuple] form. It also allows for the reverse calculation of
    indexes so that you can recreate positioning information in the original
    data set based on the new shape of the case dimension. 
    """
    
    """
    internal instance variables:
    
    bin_dimension_index   - the original index of the bin dimension
    tuple_dimension_index - the original index of the tuple dimension
    
    original_data_shape   - the shape the data was before it was reordered
    new_data_shape        - the data shape after it's been reordered
    
    new_index_order       - a mapping that lists the order of the new dimension indexes
    original_case_shape   - the shape of the case dimension(s) before being flattened
    reverse_case_index    - a reverse index for finding the original positions of
                            flattened case indexes
    
    TODO, in the long run, find a way to get rid of the reverse_case_index
    """
    
    def __init__ (self, dataShape, binIndexNumber=0, tupleIndexNumber=None) :
        """
        Given information on the original data and the desired bin/tuple,
        build the mapping object
        """
        
        # minimally, we need to have a shape
        assert(dataShape is not None)
        
        # get the number of dimensions present in our data
        numberOfDimensions = len(dataShape)
        
        # is our shape ok?
        assert(numberOfDimensions >=2)
        
        self.original_data_shape = dataShape
        
        # set up our tuple if it wasn't selected
        if (tupleIndexNumber is None) :
            tupleIndexNumber = numberOfDimensions - 1
        
        # are the bin and tuple ok?
        assert(binIndexNumber is not None)
        assert(binIndexNumber   >= 0)
        assert(binIndexNumber   < numberOfDimensions)
        assert(tupleIndexNumber >= 0)
        assert(tupleIndexNumber < numberOfDimensions)
        
        self.bin_dimension_index   = binIndexNumber
        self.tuple_dimension_index = tupleIndexNumber
        
        # get the new index ordering for the data
        self.new_index_order     = BinTupleMapping._make_new_index_list(numberOfDimensions,
                                                                        self.bin_dimension_index,
                                                                        self.tuple_dimension_index)
        temp_data_shape = [ ]
        for index in self.new_index_order:
            temp_data_shape = temp_data_shape + [dataShape[index]]
        temp_data_shape = tuple(temp_data_shape)
        """
        temp_data_shape          = np.array(dataShape).transpose(self.new_index_order)
        """
        self.original_case_shape = temp_data_shape[1:-1]
        
        # figure out the new size with the flattened cases
        number_of_cases     = 0
        self.new_data_shape = (temp_data_shape[0], temp_data_shape[-1])
        if len(self.original_case_shape) > 0 :
            number_of_cases = np.multiply.accumulate(self.original_case_shape)[-1]
            self.new_data_shape = (temp_data_shape[0], number_of_cases, temp_data_shape[-1])
        
        # build the reverse index for looking up flat case indexes
        self.reverse_case_index     = None
        if len(self.original_case_shape) > 0 :
            self.reverse_case_index = np.arange(number_of_cases).reshape(self.original_case_shape)
    
    @staticmethod
    def _make_new_index_list(numberOfIndexes, firstIndexNumber, lastIndexNumber) :
        """
        a utility method to make a list of index numbers for reordering a
        multi-dimensional array
        
        the first and last index numbers represent the dimensions you want to be
        first and last (respectively) when the list is reordered; any other indexes
        will retain their relative ordering
        
        Note: This is a private method of the BinTupleMapping class and assumes
        that the index numbers passed to it will have been preverified to be
        acceptable.
        """
        
        # make the new list
        newIndexList = range(numberOfIndexes)
        
        # remove our two "important" indexes, in the correct order
        maxSpecial   = max(firstIndexNumber, lastIndexNumber)
        minSpecial   = min(firstIndexNumber, lastIndexNumber)
        del(newIndexList[maxSpecial])
        del(newIndexList[minSpecial])
        
        # add our two important indexes back into the list in their new places
        newIndexList = [firstIndexNumber] + newIndexList + [lastIndexNumber]
        
        return newIndexList
    
    def reorder_for_bin_tuple (self, data) :
        """
        reorder the data so that the bin index is first, the tuple index is last,
        and any additional dimensions are flattened into a middle "case" index
        
        the reordered data and the shape of flattened case indexes will be returned
        (note: the shape of the data must match the shape with which the BinTupleMatching
        object was originally constructed)
        """
        
        assert(data.shape == self.original_data_shape)
        
        # put the bin and tuple dimensions in the correct places
        newData = data.transpose(self.new_index_order)
        
        # flatten the case dimensions
        newData = newData.reshape(self.new_data_shape)
        
        return newData
    
    def determine_case_indecies (self, flatIndex) :
        """
        determine the original indexes of the case from the flat case index number
        
        Note: this method requires the object to hold a large data structure
        TODO, find a better way of doing this? does numpy guarantee reshaping strategy?
        TODO, can I find information on reshape and do this with pure math?
        """
        
        if self.reverse_case_index is None :
            return None
        
        # find the flat index in our reverse case index
        positionOfIndex = np.where(self.reverse_case_index == flatIndex)
        
        return positionOfIndex

if __name__=='__main__':
    import doctest
    doctest.testmod()
