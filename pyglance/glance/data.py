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
    
    # TODO, analyze in issolation?

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
        
        # diff the two data sets TODO, this doesn't use epsilon percent yet
        raw_diff, valid_in_both, (valid_in_a_mask, valid_in_b_mask), trouble_pt_mask, outside_epsilon_mask,  \
           (a_not_finite_mask, b_not_finite_mask), (a_missing_mask, b_missing_mask), (ignore_mask_a, ignore_mask_b) = \
                            diff(aDataObject.data, bDataObject.data, epsilonValue,
                                 (aDataObject.fill_value, bDataObject.fill_value),
                                 (aDataObject.masks.ignore_mask, bDataObject.masks.ignore_mask))
        
        # set the various data in our two basic data objects
        aDataObject.masks = BasicMaskSetObject(ignore_mask_a, valid_in_a_mask, a_not_finite_mask, a_missing_mask)
        bDataObject.masks = BasicMaskSetObject(ignore_mask_b, valid_in_b_mask, b_not_finite_mask, b_missing_mask)
        
        # create our diff info object
        self.diff_data_object = DataObject(raw_diff)
        self.diff_data_object.masks = DiffMaskSetObject(ignore_mask_a | ignore_mask_b,
                                                        valid_in_both, trouble_pt_mask, outside_epsilon_mask)

# Upcasts to be used in difference computation to avoid overflow. Currently only unsigned
# ints are upcast.
# FUTURE: handle uint64s as well (there is no int128, so might have to detect overflow)
datatype_upcasts = {
    np.uint8:  np.int16,
    np.uint16: np.int32,
    np.uint32: np.int64
    }

# TODO, rethink how this works
def _select_fill_data(dTypeValue) :
    """
    select a fill data value based on the type of data that is being
    inspected/changed
    """
    
    fill_value_to_return = None
    
    if np.issubdtype(dTypeValue, np.float) or np.issubdtype(dTypeValue, np.complex) :
        fill_value_to_return = np.nan
    elif np.issubdtype(dTypeValue, np.int) :
        fill_value_to_return = np.iinfo(dTypeValue).min
    elif np.issubdtype(dTypeValue, np.bool) :
        fill_value_to_return = True
    elif ((dTypeValue is np.uint8)  or
          (dTypeValue is np.uint16) or
          (dTypeValue is np.uint32) or
          (dTypeValue is np.uint64)) :
        fill_value_to_return = np.iinfo(dTypeValue).max
    
    return fill_value_to_return

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
    assert(np.can_cast(aData.dtype, bData.dtype) or np.can_cast(bData.dtype, aData.dtype))
    
    # if the ignore masks do not exist, set them to include none of the data
    if (ignore_mask_a is None) :
        ignore_mask_a = np.zeros(shape,dtype=bool)
    if (ignore_mask_b is None) :
        ignore_mask_b = np.zeros(shape,dtype=bool)
    
    # deal with the basic masks
    a_not_finite_mask, b_not_finite_mask = ~np.isfinite(aData) & ~ignore_mask_a, ~np.isfinite(bData) & ~ignore_mask_b
    a_missing_mask, b_missing_mask = np.zeros(shape,dtype=bool), np.zeros(shape,dtype=bool)
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
        sharedType = np.common_type(aData, bData)

    # upcast if needed to avoid overflow in difference operation
    if sharedType in datatype_upcasts:
        sharedType = datatype_upcasts[sharedType]

    LOG.debug('Shared data type that will be used for diff comparison: ' + str(sharedType))
    
    # construct our diff'ed array
    raw_diff = np.zeros(shape, dtype=sharedType) #empty_like(aData)
    
    fill_data_value = _select_fill_data(sharedType)
    
    LOG.debug('current fill data value: ' + str(fill_data_value))
    
    raw_diff[~valid_in_both] = fill_data_value # throw away invalid data

    # compute difference, using shared type in computation
    raw_diff[valid_in_both] = bData[valid_in_both].astype(sharedType) - aData[valid_in_both].astype(sharedType)
        
    # the valid data which is too different between the two sets according to the given epsilon
    outside_epsilon_mask = (abs(raw_diff) > epsilon) & valid_in_both
    # trouble points = mismatched nans, mismatched missing-values, differences that are too large 
    trouble_pt_mask = (a_not_finite_mask ^ b_not_finite_mask) | (a_missing_mask ^ b_missing_mask) | outside_epsilon_mask
    
    return raw_diff, valid_in_both, (valid_in_a_mask, valid_in_b_mask), trouble_pt_mask, outside_epsilon_mask,  \
           (a_not_finite_mask, b_not_finite_mask), (a_missing_mask, b_missing_mask), (ignore_mask_a, ignore_mask_b)

if __name__=='__main__':
    import doctest
    doctest.testmod()
