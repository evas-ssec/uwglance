#!/usr/bin/env python
# encoding: utf-8
"""
General filters that can be used to prefilter data before comparison (define in the config file).

Created by Eva Schiffer August 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import numpy      as np

def trim_off_of_top (data, num_elements_to_trim) :
    """
    Remove num_elements_to_trim rows from the top of the data array
    and return a copy of the remaining data
    """
    assert(num_elements_to_trim >= 0)
    
    return data[num_elements_to_trim:, :].copy()

def trim_off_of_bottom (data, num_elements_to_trim) :
    """
    Remove num_elements_to_trim rows from the bottom of the data array
    and return a copy of the remaining data
    """
    assert(num_elements_to_trim >= 0)
    
    return data[:(-1 * num_elements_to_trim), :].copy()

def trim_off_of_right (data, num_elements_to_trim) :
    """
    Remove num_elements_to_trim columns from the right side of the data array
    and return a copy of the remaining data
    """
    assert(num_elements_to_trim >= 0)
    
    return data[:, :(-1 * num_elements_to_trim)].copy()

def trim_off_of_left (data, num_elements_to_trim) :
    """
    Remove num_elements_to_trim columns from the left side of the data array
    and return a copy of the remaining data
    """
    assert(num_elements_to_trim >= 0)
    
    return data[:, num_elements_to_trim:]

def flatten_data_into_bins (data, ranges, new_values, missing_value, return_data_type) :
    """
    Sort the data into the given ranges. Each range should correspond to a value
    given in the list of new_values; this value will be used for all data points
    that fall within that range. Ranges will be treated as a continuious spectrum
    So please list them in increasing order.
    
    Note: Data values that fall directly between two ranges will be grouped with
    the larger range.
    
    Also Note: If a data point is not found to fall within any of the ranges,
    the missing_value will be filled into that spot instead.
    """
    # make sure we have values to match each range, no more and no less
    assert(len(ranges) == (len(new_values) + 1))
    
    new_data = np.zeros(data.shape, dtype=return_data_type)
    
    # find the data in each range
    filled_mask = np.zeros(data.shape, dtype=bool)
    temp_mask   = np.zeros(data.shape, dtype=bool)
    for index in range(len(ranges) - 1) :
        temp_mask = (data >= ranges[index]) & (data <= ranges[index + 1])
        new_data[temp_mask] = new_values[index]
        filled_mask = filled_mask | temp_mask
    
    # clean up anything that didn't get filled
    new_data[~filled_mask] = missing_value
    
    return new_data

def extract_bit_from_packed_mask (data, index_of_bit_to_extract,
                                  (truth_value, false_value),
                                  return_data_type) :
    """
    Extract a one bit boolean mask from a larger packed data set. The bit that
    you wish to extract must be identified by it's index from the lower end of
    the array of bits (with 0 being the bit that would represent a value of 1
    if the mask was an integer; so the 1 indexed bit would represent the integer
    value of 2, the 2 indexed bit a value of 4, and so on).
    Note: It is assumed that the endian-ness of your data was handled correctly
    by whatever code loaded it from the file. 
    
    The return data will be filled with the truth_value and false_value based
    on whether the bit was 1 (true) or 0 (false) and will be of the requested
    return_data_type.
    
    Note: If you wish to examine or compare more than one packed bit, you may
    use this filter multiple times to extract each separately and then program
    your own function to perform any comparisons. 
    """
    # we are only allowing positive indexing due to data typing complexity
    assert(index_of_bit_to_extract >= 0)
    
    # make our mask
    mask_value = 2**index_of_bit_to_extract
    bit_mask = np.zeros(data.shape, dtype=int)
    bit_mask = bit_mask + mask_value
    
    # get the data out and fill in the requested values
    pure_extraction = np.bitwise_and(data, bit_mask)
    new_data = np.zeros(data.shape, dtype=return_data_type)
    new_data[pure_extraction >  0] = truth_value
    new_data[pure_extraction <= 0] = false_value
    
    return new_data

def select_slice_from_3D_last (data, slice_index) :
    """
    Select a slice from a 3 dimensional data set.
    slice_index indicates the index of the slice that you would like returned
    by this function
    
    note: this assumes you wish to slice across the dimention that you would index
    into last
    """
    assert(slice_index >= 0)
    assert(slice_index < data.shape[2])
    
    return data[:, :, slice_index]

def rotate_indexes_right (data) :
    """
    move the order of the indexes in the array to the right in order, taking the last one
    and putting it in the first index spot
    note: at the moment this filter only works with 3 dimentional data sets
    """
    
    # figure out the shapes we have/need
    old_shape = data.shape
    #print ('old shape: ' + str(old_shape))
    new_shape = old_shape[-1:] + old_shape[:-1]
    #print ('new shape: ' + str(new_shape))
    
    # set up our new data
    data_new = np.empty_like(data)
    data_new = data_new.reshape(new_shape)
    
    # move the old data into the new data shape
    for index1 in range(old_shape[0]) :
        for index2 in range(old_shape[1]) :
            for index3 in range(old_shape[2]) :
                data_new[index3, index1, index2] = data[index1, index2, index3]
    
    return data_new

def set_to_value_between_bounds(data, value_to_set_to, bottom_bound_exclusive, top_bound_exclusive) :
    """
    Wherever the data is non-finite or outside the given bounds, set it to the given value.
    """
    
    mask = (data < bottom_bound_exclusive) | (data > top_bound_exclusive) | (~ np.isfinite(data))
    data[mask] = value_to_set_to
    
    return data

def collapse_to_index(data, index, collapsing_function=np.mean,
                      missing_value=None, ignore_below_exclusive=None, ignore_above_exclusive=None) :
    """
    collapse the data given along the selected index, using the collapsing_function
    """
    
    # figure the order to put the index we want first
    new_index_order = range(len(data.shape))
    del(new_index_order[index])
    new_index_order = [index] + new_index_order
    
    # copy our data and put the index we won't collapse first
    new_data = data.copy()
    new_data = new_data.transpose(new_index_order)
    
    # find any bad points that we don't want to collapse
    bad_mask = np.zeros(new_data.shape, dtype=bool)
    if missing_value is not None :
        bad_mask[new_data == missing_value] = True
    if ignore_above_exclusive is not None :
        bad_mask[new_data > ignore_above_exclusive] = True
    if ignore_below_exclusive is not None :
        bad_mask[new_data < ignore_below_exclusive] = True
    
    # collapse any non bad data
    final_num_pts = new_data.shape[0]
    collapsed_data = np.zeros((final_num_pts), dtype=new_data.dtype)
    for index_1_val in range(final_num_pts) :
        values_to_use = new_data[index_1_val][~(bad_mask[index_1_val])]
        if len(values_to_use) > 0 :
            collapsed_data[index_1_val] = collapsing_function(values_to_use)
        else:
            collapsed_data[index_1_val] = missing_value
    
    return collapsed_data

def organize_ipopp_data_into_image(original_ipopp_data, wave_number=None, missing_value=None,
                                   propagate_partial_missing_values=False) :
    """
    organize the ipopp data spatially into an 'image' of sorts
    this basically consists of:
    
                      -> the 30 fields of regard
                  ---------------
                  |             |   |
                  |             |   V
                  |             |  the 4 scan lines
                  |             |
                  ---------------
                  
                  for each field of regard/scan line
                  _______
                  |_|_|_|
                  |_|_|_|
                  |_|_|_|
                  
                  block of the 9 detectors (TODO order?
                  for the moment assuming increasing linear by line then field of regard)
                  
                  for each detector point, if the wave_number was given
                  in the parameters, that specific interferogram data pt will be used
                  if no wave_number was given, the mean of the 717 pts will be used
    """
    
    new_data_image = np.zeros((4 * 3, 30 * 3), dtype=original_ipopp_data.dtype)
    
    # loop to the place in the old array to get each data point
    for scan_line in range(4) :
        for field_of_regard in range(30) :
            
            for detector in range(9) :
                      
                # figure out the value we're moving
                data_pt = None
                data_array = original_ipopp_data[scan_line][field_of_regard][detector]
                if (wave_number is not None):
                    data_pt = data_array[wave_number]
                else:
                    if (propagate_partial_missing_values
                        and
                        sum(data_array == missing_value) > 0) :
                        
                        data_pt = missing_value
                    else:
                        # remember to remove any remaining missing values
                        data_pt = np.mean(data_array[~(data_array == missing_value)])
                
                # figure out where to put the value and put it there
                index1 = scan_line * 3       + (detector / 3)
                index2 = field_of_regard * 3 + (detector % 3)
                new_data_image[index1][index2] = data_pt
    
    return new_data_image

if __name__=='__main__':
    sys.exit(main())