#!/usr/bin/env python
# encoding: utf-8
"""
This module handles statistical analysis of data sets. The code present in
this module is based on previous versions of delta.py.

Created by evas Apr 2010.
Copyright (c) 2010 University of Wisconsin SSEC. All rights reserved.
"""

import glance.data  as dataobj
import glance.delta as delta

import numpy as np

# TODO, finish transitioning to classes

# TODO, I don't like this design, but it's what I could come up
# with for now. Reconsider this again later.
class StatisticalData (object) :
    """
    This class represents a set of statistical data generated from
    the examination of two data sets. This data set is relatively
    abstract. 
    
    All Statistics Data objects should have a title and be able to provide
    a dictionary of their statistics (see dictionary_form function) and
    a dictionary documenting their statistics.
    
    Child classes can include whatever actual statistics they like.
    """
    
    def __init__ (self) :
        """
        a minimal constructor that only sets the title
        """
        
        self.title = None

class MissingValueStatistics (StatisticalData) :
    """
    A class representing information about where fill values are found
    in a pair of data sets.
    
    includes the following statistics:
    
    a_missing_count         -    count of points that are missing in the a data set
    a_missing_fraction      - fraction of points that are missing in the a data set
    b_missing_count         -    count of points that are missing in the b data set
    b_missing_fraction      - fraction of points that are missing in the b data set
    common_missing_count    -    count of points that are missing in both data sets
    common_missing_fraction - fraction of points that are missing in both data sets
    """
    

# --------------------- general statistics methods ------------------

def _get_missing_value_stats(a_missing_mask, b_missing_mask) :
    """
    Get a list of statistics about missing data values in data sets a and b,
    given masks describing the positions of the missing data (such masks
    can be gotten from the diff function), 
    the return value will be a dictionary of statistics
    """
    # calculate information about the missing values
    num_a_missing = np.sum(a_missing_mask)
    num_b_missing = np.sum(b_missing_mask)
    num_common_missing = np.sum(a_missing_mask & b_missing_mask)
    
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

def summarize(a, b, epsilon=0.,
              (a_missing_value, b_missing_value)=(None,None),
              ignoreInAMask=None, ignoreInBMask=None,
              epsilonPercent=None):
    """return dictionary of statistics dictionaries
    stats not including 'nan' in name exclude nans in either arrays
    """
    # diff our two data sets
    aDataObject = dataobj.DataObject(a, fillValue=a_missing_value, ignoreMask=ignoreInAMask)
    bDataObject = dataobj.DataObject(b, fillValue=b_missing_value, ignoreMask=ignoreInBMask)
    diffInfo = dataobj.DiffInfoObject(aDataObject, bDataObject,
                                      epsilonValue=epsilon, epsilonPercent=epsilonPercent) 
    #TODO, for the moment, unpack these values into local variables
    diffData = diffInfo.diff_data_object.data
    finite_mask    = diffInfo.diff_data_object.masks.valid_mask
    finite_a_mask = diffInfo.a_data_object.masks.valid_mask
    finite_b_mask = diffInfo.b_data_object.masks.valid_mask
    trouble = diffInfo.diff_data_object.masks.trouble_mask
    outside_epsilon = diffInfo.diff_data_object.masks.outside_epsilon_mask
    anfin = diffInfo.a_data_object.masks.non_finite_mask
    bnfin = diffInfo.b_data_object.masks.non_finite_mask
    amis   = diffInfo.a_data_object.masks.missing_mask
    bmis   = diffInfo.b_data_object.masks.missing_mask
    ignoreInAMask = diffInfo.a_data_object.masks.ignore_mask
    ignoreInBMask = diffInfo.b_data_object.masks.ignore_mask
    
    general_stats = _get_general_data_stats(a, b, a_missing_value, b_missing_value, epsilon, 
                                            ignoreInAMask, ignoreInBMask, ~finite_a_mask, ~finite_b_mask) 
    additional_statistics = stats(diffData, finite_mask) #*nfo) # grab some additional comparison statistics 
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

def stats(diffData, mask, *etc):
    
    # if there are no values after the mask,
    # we can't do any of these forms of
    # statistical analysis
    if np.sum(mask) <= 0 :
        return { }
    
    absDiffData = abs(diffData)
    root_mean_square_value = delta.calculate_root_mean_square(diffData, mask)
    return {    'rms_diff': root_mean_square_value, 
                'std_diff':       np.std(absDiffData[mask]), 
                'mean_diff':     np.mean(absDiffData[mask]), 
                'median_diff': np.median(absDiffData[mask]),
                'max_diff':       np.max(absDiffData[mask])
                }

def _get_num_perfect(a, b, ignoreMask=None):
    numPerfect = 0
    if not (ignoreMask is None) :
        numPerfect = np.sum(a[~ignoreMask] == b[~ignoreMask])
    else :
        numPerfect = np.sum(a == b)
    return numPerfect

def _get_numerical_data_stats(a, b, diff_data,  data_is_finite_mask,
                              outside_epsilon_mask, trouble_mask,
                              additional_statistics={}) : 
    """
    Get a list of numerical comparison related statistics about a and b,
    given a and b and some other information about them.
    the return value will be a dictionary of statistics
    """
    # calculate our various statistics
    num_finite_values_too_different = np.sum(outside_epsilon_mask)
    num_perfect = _get_num_perfect(a, b, ~data_is_finite_mask)
    r_corr = delta.compute_correlation(a, b, data_is_finite_mask)
    num_trouble = np.sum(trouble_mask)
    
    # we actually want the total number of _finite_ values rather than all the data
    total_num_finite_values = np.sum(data_is_finite_mask)
    
    # no dividing by 0!
    fraction_too_different = 0.0
    fraction_perfect = 0.0
    if total_num_finite_values > 0 :
        fraction_too_different = num_finite_values_too_different / float(total_num_finite_values)
        fraction_perfect = num_perfect / float(total_num_finite_values)
    
    comparison = {  'correlation': r_corr,
                    'r-squared correlation': r_corr * r_corr,
                    'diff_outside_epsilon_count': num_finite_values_too_different,
                    'diff_outside_epsilon_fraction': fraction_too_different,
                    'perfect_match_count': num_perfect,
                    'perfect_match_fraction': fraction_perfect,
                     'trouble_points_count': num_trouble, 
                     'trouble_points_fraction': float(num_trouble) / float(a.size)
                    }
    comparison.update(additional_statistics)
    
    return comparison

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
    num_ignored_in_a = np.sum(spatial_ignore_in_a_mask)
    num_ignored_in_b = np.sum(spatial_ignore_in_b_mask)
    
    # get the number of data points
    total_num_values = a.size
    
    general_stats = {'a_missing_value': a_missing_value,
                     'b_missing_value': b_missing_value,
                     'epsilon': epsilon,
                     'max_a': delta.max_with_mask(a, ~bad_in_a),
                     'max_b': delta.max_with_mask(b, ~bad_in_b),
                     'min_a': delta.min_with_mask(a, ~bad_in_a),
                     'min_b': delta.min_with_mask(b, ~bad_in_b),
                     'num_data_points': total_num_values,
                     'shape': a.shape,
                     'spatially_invalid_pts_ignored_in_a': num_ignored_in_a,
                     'spatially_invalid_pts_ignored_in_b': num_ignored_in_b
                     }
    
    return general_stats

def _get_finite_data_stats(a_is_finite_mask, b_is_finite_mask, common_ignore_mask) :
    """
    Get a list of statistics about finite data values in data sets a and b,
    given masks describing the positions of the finite data in each file (such masks
    can be gotten from the diff function),
    the return value will be a dictionary of statistics
    """
    # calculate information about the finite data
    num_a_finite = np.sum(a_is_finite_mask)
    num_b_finite = np.sum(b_is_finite_mask)
    num_common_finite = np.sum(a_is_finite_mask & b_is_finite_mask)
    num_finite_in_only_one = np.sum((a_is_finite_mask ^ b_is_finite_mask) & (~ common_ignore_mask)) # use an exclusive OR
    
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



def _get_nan_stats(a_nan_mask, b_nan_mask) :
    """
    Get a list of statistics about non-numerical values in data sets a and b,
    the return value will be a dictionary of statistics
    """
    # find the nan values in the data
    num_a_nans = np.sum(a_nan_mask)
    num_b_nans = np.sum(b_nan_mask)
    num_common_nans = np.sum(a_nan_mask & b_nan_mask)
    
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



# -------------------------- documentation -----------------------------

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
                    'r-squared correlation': "the square of the r correlation (see correlation)",
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
