#!/usr/bin/env python
# encoding: utf-8
"""
Data objects for use in glance

Created by evas Apr 2010.
Copyright (c) 2010 University of Wisconsin SSEC. All rights reserved.
"""

import logging
import numpy as np

import glance.delta as delta

LOG = logging.getLogger(__name__)

class BasicMaskSetObject (object) :
    """
    This class represents a basic set of masks that a data set may have.
    The set must contain an "ignore" mask, and may optionally contain others.
    Note: This item is intended to be read only. If you want to change a mask,
    create a new one with the new masks
    
    ignore_mask - a mask of data that should be ignored for reasons not related to
                  the contents of the actual data set (generally longitude or
                  latitude issues)
    valid_mask - a mask of "good" values (ie. finite and non-missing)
    non_finite_mask - a mask of non-finite values
    missing_mask - a mask of where the data's fill value is present instead of
                   actual data values
    """
    
    def __init__(self, ignoreMask,
                 validMask=None, nonFiniteMask=None, missingMask=None) :
        """
        create the mask set with at least the ignore mask
        (the others are optional)
        """
        self._reset_all_masks()
        
        self.ignore_mask     = ignoreMask
        self.valid_mask      = validMask
        self.non_finite_mask = nonFiniteMask
        self.missing_mask    = missingMask
    
    def _reset_all_masks(self) :
        """
        set all the masks to None
        """
        self.ignore_mask     = None
        self.valid_mask      = None
        self.non_finite_mask = None
        self.missing_mask    = None

class DiffMaskSetObject (BasicMaskSetObject) :
    """
    This class represents a set of masks that related to two or more
    compared data sets. The inherited ignore/valid/non-finite/missing
    masks are used to capture information about where data should be
    ignored/valid/non-finite in the compared data.
    
    Additionally, masks describing trouble points and points outside
    of the epsilon analysis tolerances are included.
    
    trouble_mask - a mask of data points which may indicate issues in the data
    outside_epsilon_mask - a mask of points which did not pass the epsilon
                           tolerance testing
    """
    
    def __init__(self, ignoreMask, validInBothMask, troubleMask, epsilonMask) :
        """
        create a more complex mask, including additional difference information
        """
        self._reset_all_masks()
        
        self.ignore_mask          = ignoreMask
        self.valid_mask           = validInBothMask
        self.trouble_mask         = troubleMask
        self.outside_epsilon_mask = epsilonMask

class DataObject (object) :
    """
    This class represents a data set.
    It may include a multidimentional numpy array of data
    as well as the fill value and a set of masks that apply to this data.
    
    data       - the raw array of data (generally this should be a numpy array)
    fill_value - the fill value used in the data array
    masks      - the set of masks that apply to this data
    """
    
    def __init__(self, dataArray, fillValue=None, ignoreMask=None) :
        """
        Create the data object.
        
        The array of data is expected to be a numpy array.
        The fill value and mask sets are optional.
        If the fill value is provided it is expected to be of the same
        data type as the data array.
        """
        
        # TODO, add some assertions for our expectations
        
        self.data       = dataArray
        self.fill_value = fillValue
        self.masks      = BasicMaskSetObject(ignoreMask)
    
    def self_analysis(self) :
        """
        Gather some basic information about a data set
        """
        
        # hang onto the shape for convenience
        shape = self.data.shape
        
        # if there isn't an ignore mask, make an empty one
        if self.masks.ignore_mask is None :
            self.masks.ignore_mask = np.zeros(shape, dtype=np.bool)
        
        # find the non-finite values
        non_finite_mask = ~np.isfinite(self.data) & ~self.masks.ignore_mask
        
        # find and mark the missing values
        missing_mask    = np.zeros(shape, dtype=np.bool)
        # if the data has a fill value, mark where the missing data is
        if self.fill_value is not None :
            missing_mask[self.data == self.fill_value] = True
            missing_mask[self.masks.ignore_mask]       = False
        
        # define the valid mask as places where the data is not missing,
        # nonfinite, or ignored
        valid_mask = ~ (missing_mask | non_finite_mask | self.masks.ignore_mask)
        
        # set our masks
        self.masks = BasicMaskSetObject(self.masks.ignore_mask, valid_mask,
                                        non_finite_mask, missing_mask)

class DiffInfoObject (object) :
    """
    This class represents the full difference between two data sets.
    
    a_data_object    - data object describing the A data set
    b_data_object    - data object describing the B data set
    diff_data_object - data object describing the raw differences between A and B
    
    epsilon_value    - the epsilon value used for comparison or None
    epsilon_percent  - the percentage (of A) used for epsilon comparisons or None
    (if both a value and percent are present, two epsilon tests will be done)
    """
    
    def __init__(self, aDataObject, bDataObject,
                 epsilonValue=0.0, epsilonPercent=None) :
        """
        analyze the difference between these two data sets at the
        given epsilon values
        """
        
        # set the basic values
        self.a_data_object   = aDataObject
        self.b_data_object   = bDataObject
        self.epsilon_value   = epsilonValue
        self.epsilon_percent = epsilonPercent
        
        # analyze our data and get the difference object
        self.diff_data_object = DiffInfoObject.analyze(aDataObject, bDataObject,
                                                       epsilonValue, epsilonPercent)
    
    # Upcasts to be used in difference computation to avoid overflow. Currently only unsigned
    # ints are upcast.
    # FUTURE: handle uint64s as well (there is no int128, so might have to detect overflow)
    DATATYPE_UPCASTS = {
        np.uint8:  np.int16,
        np.uint16: np.int32,
        np.uint32: np.int64
        }
    
    @staticmethod
    def _get_shared_type_and_fill_value(data1, data2, fill1=None, fill2=None) :
        """
        Figure out a shared type that can be used when adding or subtracting
        the two data sets given (accounting for possible overflow)
        Also returns a fill value that can be used.
        """
        
        # figure out the shared type
        type_to_return = data1.dtype
        if data1.dtype is not data2.dtype:
            type_to_return = np.common_type(data1, data2)
        
        # upcast the type if we need to
        if type_to_return in DiffInfoObject.DATATYPE_UPCASTS :
            type_to_return = DiffInfoObject.DATATYPE_UPCASTS[type_to_return]
            LOG.debug('To prevent overflow, difference data will be upcast from ('
                      + str(data1.dtype) + '/' + str(data2.dtype) + ') to: ' + str(type_to_return))
        
        # figure out the fill value
        fill_value_to_return = None
        
        # if we're looking at float or complex data, use a nan
        if (np.issubdtype(type_to_return, np.float) or
            np.issubdtype(type_to_return, np.complex)) :
            fill_value_to_return = np.nan
        
        # if we're looking at int data, use the minimum value
        elif np.issubdtype(type_to_return, np.int) :
            fill_value_to_return = np.iinfo(type_to_return).min
        
        # if we're looking at unsigned data, use the maximum value
        elif ((type_to_return is np.uint8)  or
              (type_to_return is np.uint16) or
              (type_to_return is np.uint32) or
              (type_to_return is np.uint64)) :
            fill_value_to_return = np.iinfo(type_to_return).max
        
        return type_to_return, fill_value_to_return
    
    @staticmethod
    def analyze(aDataObject, bDataObject,
                epsilonValue=0.0, epsilonPercent=None):
        """
        analyze the differences between the two data sets
        updates the two data objects with additional masks
        and returns data object containing diff data and masks
        """
        shape = aDataObject.data.shape
        assert(bDataObject.data.shape == shape)
        assert(np.can_cast(aDataObject.data.dtype, bDataObject.data.dtype) or
               np.can_cast(bDataObject.data.dtype, aDataObject.data.dtype))
        
        # do some basic analysis on the individual data sets
        aDataObject.self_analysis()
        bDataObject.self_analysis()
        
        # where is the shared valid data?
        valid_in_both  = aDataObject.masks.valid_mask  & bDataObject.masks.valid_mask
        ignore_in_both = aDataObject.masks.ignore_mask | bDataObject.masks.ignore_mask
        
        # get our shared data type and fill value
        sharedType, fill_data_value = DiffInfoObject._get_shared_type_and_fill_value(aDataObject.data,
                                                                                     bDataObject.data,
                                                                                     aDataObject.fill_value,
                                                                                     bDataObject.fill_value)
        
        # construct our diff'ed data set
        raw_diff = np.zeros(shape, dtype=sharedType)
        raw_diff[~valid_in_both] = fill_data_value # throw away invalid data
        # compute difference, using shared type in computation
        raw_diff[valid_in_both] = bDataObject.data[valid_in_both].astype(sharedType) -  \
                                  aDataObject.data[valid_in_both].astype(sharedType)
        
        # the valid data which is too different between the two sets according to the given epsilon
        outside_epsilon_mask = (abs(raw_diff) > epsilonValue) & valid_in_both
        # trouble points = mismatched nans, mismatched missing-values, differences that are too large 
        trouble_pt_mask = ( (aDataObject.masks.non_finite_mask ^ bDataObject.masks.non_finite_mask) |
                            (aDataObject.masks.missing_mask    ^ bDataObject.masks.missing_mask)    |
                            outside_epsilon_mask )
        
        # make our diff data object
        diff_data_object = DataObject(raw_diff, fillValue=fill_data_value)
        diff_data_object.masks = DiffMaskSetObject(ignore_in_both, valid_in_both,
                                                   trouble_pt_mask, outside_epsilon_mask)
        
        return diff_data_object

if __name__=='__main__':
    import doctest
    doctest.testmod()
