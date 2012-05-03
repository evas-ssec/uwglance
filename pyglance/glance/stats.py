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
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        
        note: child classes should override this method
        """
        return { }
    
    @staticmethod
    def doc_strings( ) :
        """
        get documentation strings that match the
        dictionary form of the statistics this class
        creates
        
        note: child classes should override this method
        """
        return { }

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
    
    _doc_strings = {
                    'a_missing_count':         "number of values flagged missing in A",
                    'a_missing_fraction':      "fraction of values flagged missing in A",
                    'b_missing_count':         "number of values flagged missing in B",
                    'b_missing_fraction':      "fraction of values flagged missing in B",
                    'common_missing_count':    "number of missing values in common between A and B",
                    'common_missing_fraction': "fraction of missing values in common between A and B"
                    }
    
    def __init__(self, diffInfoObject) :
        """
        build our fill value related statistics based on the comparison
        of two data sets
        """
        self.title = 'Missing Value Statistics'
        
        # pull out some masks for later use
        a_missing_mask = diffInfoObject.a_data_object.masks.missing_mask
        b_missing_mask = diffInfoObject.b_data_object.masks.missing_mask
        
        assert(a_missing_mask.shape == b_missing_mask.shape)
        
        # figure out some basic statistics
        self.a_missing_count      = np.sum(a_missing_mask)
        self.b_missing_count      = np.sum(b_missing_mask)
        self.common_missing_count = np.sum(a_missing_mask & b_missing_mask)
        
        # make the assumption that a and b are the same size and only use the size of a's mask
        total_num_values = a_missing_mask.size
        
        # figure out some fraction statistics
        self.a_missing_fraction      = float(self.a_missing_count)      / float(total_num_values)
        self.b_missing_fraction      = float(self.b_missing_count)      / float(total_num_values)
        self.common_missing_fraction = float(self.common_missing_count) / float(total_num_values)
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        
        toReturn = {
                    'a_missing_count':         self.a_missing_count,
                    'a_missing_fraction':      self.a_missing_fraction,
                    'b_missing_count':         self.b_missing_count,
                    'b_missing_fraction':      self.b_missing_fraction,
                    'common_missing_count':    self.common_missing_count,
                    'common_missing_fraction': self.common_missing_fraction
                    }
        
        return toReturn
    
    @staticmethod
    def doc_strings( ) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        return MissingValueStatistics._doc_strings

class MissingValueInspectionStatistics (StatisticalData) :
    """
    A class representing information about where fill values are found
    in a data.
    
    includes the following statistics:
    
    missing_count         -    count of points that are missing in the a data set
    missing_fraction      - fraction of points that are missing in the a data set
    """
    
    _doc_strings = {
                    'missing_count':         "number of values flagged missing",
                    'missing_fraction':      "fraction of values flagged missing",
                    }
    
    def __init__(self, dataObject) :
        """
        build our fill value related statistics based on the data set
        """
        self.title = 'Missing Value Statistics'
        
        # pull out a mask for later use
        missing_mask          = dataObject.masks.missing_mask
        
        # figure out some basic statistics
        self.missing_count    = np.sum(missing_mask)
        self.missing_fraction = float(self.missing_count) / float(missing_mask.size)
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        
        toReturn = {
                    'missing_count':         self.missing_count,
                    'missing_fraction':      self.missing_fraction,
                    }
        
        return toReturn
    
    @staticmethod
    def doc_strings( ) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        return MissingValueInspectionStatistics._doc_strings

class FiniteDataStatistics (StatisticalData) :
    """
    A class representing information about where finite values are found
    in a pair of data sets.
    
    includes the following statistics:
    
    a_finite_count              - the number   of finite data values in the a data set
    a_finite_fraction           - the fraction of finite data values in the a data set
    b_finite_count              - the number   of finite data values in the b data set
    b_finite_fraction           - the fraction of finite data values in the b data set
    common_finite_count         - the number   of finite values the two data sets have in common
    common_finite_fraction      - the fraction of finite values the two data sets have in common
    finite_in_only_one_count    - the number   of points that are finite in only one of the two sets
    finite_in_only_one_fraction - the fraction of points that are finite in only one of the two sets
    """
    
    _doc_strings = {
                    'a_finite_count': "number of finite values in A",
                    'a_finite_fraction': "fraction of finite values in A (out of all data points in A)",
                    'b_finite_count': "number of finite values in B",
                    'b_finite_fraction': "fraction of finite values in B (out of all data points in B)",
                    'common_finite_count': "number of finite values in common between A and B",
                    'common_finite_fraction': "fraction of finite values in common between A and B",
                    'finite_in_only_one_count': "number of values that changed finite-ness between A and B; " +
                                                "only the common spatially valid area is considerd for this statistic",
                    'finite_in_only_one_fraction': "fraction of values that changed finite-ness between A and B; " +
                                                "only the common spatially valid area is considerd for this statistic"
                    }
    
    def __init__(self, diffInfoObject) :
        """
        build our finite data related statistics based on the comparison
        of two data sets
        """
        self.title = 'Finite Data Statistics'
        
        # pull out some data we will use later
        a_is_finite_mask   = diffInfoObject.a_data_object.masks.valid_mask
        b_is_finite_mask   = diffInfoObject.b_data_object.masks.valid_mask
        common_ignore_mask = diffInfoObject.diff_data_object.masks.ignore_mask
        
        assert(a_is_finite_mask.shape == b_is_finite_mask.shape)
        assert(b_is_finite_mask.shape == common_ignore_mask.shape)
        
        # figure out some basic statistics
        self.a_finite_count = np.sum(a_is_finite_mask)
        self.b_finite_count = np.sum(b_is_finite_mask)
        self.common_finite_count = np.sum(a_is_finite_mask & b_is_finite_mask)
        # use an exclusive or to check which points are finite in only one of the two data sets
        self.finite_in_only_one_count = np.sum((a_is_finite_mask ^ b_is_finite_mask) & ~common_ignore_mask)
        
        # make the assumption that a and b are the same size and only use the size of a's mask
        total_num_values = a_is_finite_mask.size
        
        # calculate some fractional statistics
        self.a_finite_fraction           = float(self.a_finite_count)           / float(total_num_values)
        self.b_finite_fraction           = float(self.b_finite_count)           / float(total_num_values)
        self.common_finite_fraction      = float(self.common_finite_count)      / float(total_num_values)
        self.finite_in_only_one_fraction = float(self.finite_in_only_one_count) / float(total_num_values)
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        
        toReturn = {
                    'a_finite_count':              self.a_finite_count,
                    'a_finite_fraction':           self.a_finite_fraction,
                    'b_finite_count':              self.b_finite_count,
                    'b_finite_fraction':           self.b_finite_fraction,
                    'common_finite_count':         self.common_finite_count,
                    'common_finite_fraction':      self.common_finite_fraction,
                    'finite_in_only_one_count':    self.finite_in_only_one_count,
                    'finite_in_only_one_fraction': self.finite_in_only_one_fraction,
                    }
        
        return toReturn
    
    @staticmethod
    def doc_strings( ) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        return FiniteDataStatistics._doc_strings

class FiniteDataInspectionStatistics (StatisticalData) :
    """
    A class representing information about where finite values are found
    in a data set.
    
    includes the following statistics:
    
    finite_count              - the number   of finite data values in the data set
    finite_fraction           - the fraction of finite data values in the data set
    """
    
    _doc_strings = {
                    'finite_count': "number of finite values",
                    'finite_fraction': "fraction of finite values (out of all data points in set)",
                    }
    
    def __init__(self, dataObject) :
        """
        build our finite data related statistics based on the data set
        """
        self.title = 'Finite Data Statistics'
        
        # pull out some data we will use later
        is_finite_mask       = dataObject.masks.valid_mask
        
        # figure out some basic statistics
        self.finite_count    = np.sum(is_finite_mask)
        self.finite_fraction = float(self.finite_count) / float(is_finite_mask.size)
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        
        toReturn = {
                    'finite_count':              self.finite_count,
                    'finite_fraction':           self.finite_fraction,
                    }
        
        return toReturn
    
    @staticmethod
    def doc_strings( ) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        return FiniteDataInspectionStatistics._doc_strings

class NotANumberStatistics (StatisticalData) :
    """
    A class representing information about where non-finite values are found
    in a pair of data sets.
    
    includes the following statistics:
    
    a_nan_count         - the number   of non finite values that are present in the a data set
    a_nan_fraction      - the fraction of non finite values that are present in the a data set
    b_nan_count         - the number   of non finite values that are present in the b data set
    b_nan_fraction      - the fraction of non finite values that are present in the b data set
    common_nan_count    - the number   of non finite values that are shared between the data sets
    common_nan_fraction - the fraction of non finite values that are shared between the data sets
    """
    
    _doc_strings = {
                    'a_nan_count': "number of NaNs in A",
                    'a_nan_fraction': "fraction of NaNs in A",
                    'b_nan_count': "number of NaNs in B",
                    'b_nan_fraction': "fraction of NaNs in B",
                    'common_nan_count': "number of NaNs in common between A and B",
                    'common_nan_fraction': "fraction of NaNs in common between A and B"
                    }
    
    def __init__(self, diffInfoObject) :
        """
        build our nonfinite data related statistics based on the comparison
        of two data sets
        """
        self.title = 'NaN Statistics'
        
        # pull out some masks we will use
        a_nan_mask = diffInfoObject.a_data_object.masks.non_finite_mask
        b_nan_mask = diffInfoObject.b_data_object.masks.non_finite_mask
        
        assert(a_nan_mask.shape == b_nan_mask.shape)
        
        # get some basic statistics
        self.a_nan_count      = np.sum(a_nan_mask)
        self.b_nan_count      = np.sum(b_nan_mask)
        self.common_nan_count = np.sum(a_nan_mask & b_nan_mask)
        
        # make the assumption that a and b are the same size and only use the size of a
        total_num_values = a_nan_mask.size
        
        # calculate some fractional statistics
        self.a_nan_fraction      = float(self.a_nan_count)      / float(total_num_values)
        self.b_nan_fraction      = float(self.b_nan_count)      / float(total_num_values)
        self.common_nan_fraction = float(self.common_nan_count) / float(total_num_values)
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        
        toReturn = {
                    'a_nan_count':         self.a_nan_count,
                    'a_nan_fraction':      self.a_nan_fraction,
                    'b_nan_count':         self.b_nan_count,
                    'b_nan_fraction':      self.b_nan_fraction,
                    'common_nan_count':    self.common_nan_count,
                    'common_nan_fraction': self.common_nan_fraction
                    }
        
        return toReturn
    
    @staticmethod
    def doc_strings( ) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        return NotANumberStatistics._doc_strings

class NotANumberInspectionStatistics (StatisticalData) :
    """
    A class representing information about where non-finite values are found
    in a data set.
    
    includes the following statistics:
    
    nan_count         - the number   of non finite values that are present in the data set
    nan_fraction      - the fraction of non finite values that are present in the data set
    """
    
    _doc_strings = {
                    'nan_count': "number of NaNs",
                    'nan_fraction': "fraction of NaNs",
                    }
    
    def __init__(self, dataObject) :
        """
        build our nonfinite data related statistics based on the data set
        """
        self.title = 'NaN Statistics'
        
        # pull out a mask we will use
        nan_mask = dataObject.masks.non_finite_mask
        
        # get some basic statistics
        self.nan_count = np.sum(nan_mask)
        self.nan_fraction = float(self.nan_count) / float(nan_mask.size)
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        
        toReturn = {
                    'nan_count':         self.nan_count,
                    'nan_fraction':      self.nan_fraction,
                    }
        
        return toReturn
    
    @staticmethod
    def doc_strings( ) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        return NotANumberInspectionStatistics._doc_strings

class GeneralStatistics (StatisticalData) :
    """
    A class representing general information about a pair of data sets.
    
    includes the following statistics:
    
    a_missing_value - the fill data value in the a set
    b_missing_value - the fill data value in the b set
    epsilon         - the fixed epsilon value
    epsilon_percent - the percentage of the a set that will be used for comparison
    max_a           - the maximum value in the a set
    max_b           - the maximum value in the b set
    min_a           - the minimum value in the a set
    min_b           - the minimum value in the b set
    num_data_points - the total number of data points in each of the sets
    shape           - the shape of each of the data sets
    spatially_invalid_pts_ignored_in_a - number of points corresponding to invalid lat/lon in a set
    spatially_invalid_pts_ignored_in_b - number of points corresponding to invalid lat/lon in b set
    """
    
    _doc_strings = {
                    'a_missing_value': 'the value that is considered \"missing\" data when it is found in A',
                    'b_missing_value': 'the value that is considered \"missing\" data when it is found in B',
                    'epsilon': 'amount of difference between matching data points in A and B that is considered acceptable',
                    'epsilon_percent': 'the percentage of difference (of A\'s value) that is acceptable between A and B (optional)',
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
                    }
    
    def __init__(self, diffInfoObject) :
        """
        build our general statistics based on the comparison
        of two data sets
        """
        self.title = 'General Statistics'
        
        # pull out some masks for later use
        a_missing_mask   = diffInfoObject.a_data_object.masks.missing_mask
        b_missing_mask   = diffInfoObject.b_data_object.masks.missing_mask
        ignore_in_a_mask = diffInfoObject.a_data_object.masks.ignore_mask
        ignore_in_b_mask = diffInfoObject.b_data_object.masks.ignore_mask
        good_in_a_mask   = diffInfoObject.a_data_object.masks.valid_mask
        good_in_b_mask   = diffInfoObject.b_data_object.masks.valid_mask
        
        assert(a_missing_mask.shape   ==   b_missing_mask.shape)
        assert(b_missing_mask.shape   == ignore_in_a_mask.shape)
        assert(ignore_in_a_mask.shape == ignore_in_b_mask.shape)
        assert(ignore_in_b_mask.shape ==   good_in_a_mask.shape)
        assert(good_in_a_mask.shape   ==   good_in_b_mask.shape)
        
        # get the number of data points
        total_num_values = a_missing_mask.size
        
        # fill in our statistics
        self.a_missing_value = diffInfoObject.a_data_object.select_fill_value()
        self.b_missing_value = diffInfoObject.b_data_object.select_fill_value()
        self.epsilon         = diffInfoObject.epsilon_value
        self.epsilon_percent = diffInfoObject.epsilon_percent
        self.max_a           = delta.max_with_mask(diffInfoObject.a_data_object.data, good_in_a_mask)
        self.min_a           = delta.min_with_mask(diffInfoObject.a_data_object.data, good_in_a_mask)
        self.max_b           = delta.max_with_mask(diffInfoObject.b_data_object.data, good_in_b_mask)
        self.min_b           = delta.min_with_mask(diffInfoObject.b_data_object.data, good_in_b_mask)
        self.num_data_points = total_num_values
        self.shape           = a_missing_mask.shape
        # also calculate the invalid points
        self.spatially_invalid_pts_ignored_in_a = np.sum(ignore_in_a_mask)
        self.spatially_invalid_pts_ignored_in_b = np.sum(ignore_in_b_mask)
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        
        toReturn = {
                    'a_missing_value': self.a_missing_value,
                    'b_missing_value': self.b_missing_value,
                    'epsilon':         self.epsilon,
                    'epsilon_percent': self.epsilon_percent,
                    'max_a':           self.max_a,
                    'max_b':           self.max_b,
                    'min_a':           self.min_a,
                    'min_b':           self.min_b,
                    'num_data_points': self.num_data_points,
                    'shape':           self.shape,
                    'spatially_invalid_pts_ignored_in_a': self.spatially_invalid_pts_ignored_in_a,
                    'spatially_invalid_pts_ignored_in_b': self.spatially_invalid_pts_ignored_in_b
                    }
        
        return toReturn
    
    @staticmethod
    def doc_strings( ) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        return GeneralStatistics._doc_strings

class GeneralInspectionStatistics (StatisticalData) :
    """
    A class representing general information about a data set.
    
    includes the following statistics:
    
    missing_value   - the fill data value
    max             - the maximum value
    min             - the minimum value
    num_data_points - the total number of data points
    shape           - the shape of the data
    spatially_invalid_pts_ignored - number of points corresponding to invalid lat/lon in the set
                                    (optional if no /lon lat mapped)
    """
    
    _doc_strings = {
                    'missing_value': 'the value that is considered \"missing\" data when it is found in the data',
                    'max': 'the maximum finite, non-missing value found in the data',
                    'min': 'the minimum finite, non-missing value found in the data',
                    'num_data_points': "number of data points (may be valid or invalid data)",
                    'shape': "shape of the data",
                    'spatially_invalid_pts_ignored': 'number of points with invalid latitude/longitude information ' +
                                                     'in the data that were' +
                                                     ' ignored for the purposes of data analysis and presentation',
                    }
    
    def __init__(self, dataObject) :
        """
        build our general statistics based on the data set
        """
        self.title = 'General Statistics'
        
        # pull out some masks for later use
        missing_mask = dataObject.masks.missing_mask
        ignore_mask  = dataObject.masks.ignore_mask
        good_mask    = dataObject.masks.valid_mask
        
        #assert(missing_mask.shape == ignore_mask.shape)
        #assert(ignore_mask.shape  == good_mask.shape  )
        
        # get the number of data points
        total_num_values = missing_mask.size
        
        # fill in our statistics
        self.missing_value   = dataObject.select_fill_value()
        self.max             = delta.max_with_mask(dataObject.data, good_mask)
        self.min             = delta.min_with_mask(dataObject.data, good_mask)
        self.num_data_points = total_num_values
        self.shape           = missing_mask.shape
        # also calculate the invalid points
        self.spatially_invalid_pts_ignored = np.sum(ignore_mask)
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        
        toReturn = {
                    'missing_value':   self.missing_value,
                    'max':             self.max,
                    'max':             self.max,
                    'num_data_points': self.num_data_points,
                    'shape':           self.shape,
                    'spatially_invalid_pts_ignored': self.spatially_invalid_pts_ignored,
                    }
        
        return toReturn
    
    @staticmethod
    def doc_strings( ) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        return GeneralInspectionStatistics._doc_strings

class NumericalComparisonStatistics (StatisticalData) :
    """
    A class representing more complex comparisons between a pair of data sets.
    
    includes the following statistics:
    
    correlation                   - the Pearson correlation r-coefficient from comparing finite values of the sets
    r_squared_correlation         - the square of the correlation
    diff_outside_epsilon_count    - the number   of points that fall outside the acceptable epsilon settings
    diff_outside_epsilon_fraction - the fraction of points that fall outside the acceptable epsilon settings
    perfect_match_count           - the number   of points that match perfectly between the sets
    perfect_match_fraction        - the fraction of points that match perfectly between the sets
    mismatch_points_count         - the number   of points that have possible issues according to the current analysis
    mismatch_points_fraction      - the fraction of points that have possible issues according to the current analysis
    
    It may also contain additional statistics. This is indicated by the does_include_simple boolean.
    The possible additional statistics include:
    
    rms_val      -  the root mean squared of the          difference between the two data sets
    std_val      - the standard deviation of the          difference between the two data sets
    mean_diff    -               the mean of the absolute difference between the two data sets
    median_diff  -             the median of the absolute difference between the two data sets
    max_diff     -            the maximum of the absolute difference between the two data sets
    mean_delta   -               the mean of the          difference between the two data sets
    median_delta -             the median of the          difference between the two data sets
    max_delta    -            the maximum of the          difference between the two data sets
    min_delta    -            the minimum of the          difference between the two data sets
    
    These statistics can also be generated separately in dictionary form by calling the
    basic_analysis method on this class.
    """
    
    _doc_strings = {
                    'correlation': "Pearson correlation r-coefficient (0.0-1.0) for finite values of A and B",
                    'diff_outside_epsilon_count': "number of finite differences falling outside acceptable epsilon definitions; " +
                                            "note: this value includes data excluded by both epsilon and epsilon_percent if " +
                                            "both have been defined",
                    'diff_outside_epsilon_fraction': "fraction of finite differences falling outside acceptable epsilon " +
                                            "definitions (out of common_finite_count)",
                    'max_diff': "maximum absolute valued difference of the finite values",
                    'mean_diff': "mean of the absolute value difference of the finite values",
                    'median_diff': "median of the absolute value difference of the finite values",
                    
                    'mean_delta':      "mean of the subtractive difference of the finite values", 
                    'median_delta':    "median of the subtractive difference of the finite values",
                    'max_delta':       "maximum finite data value from the data set of B file - A file",
                    'min_delta':       "minimum finite data value from the data set of B file - A file",
                    
                    'perfect_match_count': "number of perfectly matched finite data points between A and B",
                    'perfect_match_fraction': "fraction of finite values perfectly matching between A and B (out of common_finite_count)",
                    'rms_val': "root mean square (RMS) difference of finite values",
                    'r-squared correlation': "the square of the r correlation (see correlation)",
                    'std_val': "standard deviation of difference of finite values",
                    'mismatch_points_count': 'number of points that differ in finite/missing status between the input data sets A and B,' +
                                            ' or are unacceptable when compared according to the current epsilon definitions',
                    'mismatch_points_fraction': 'fraction of points that differ in finite/missing status between the input data sets A and B,' +
                                            ' or are unacceptable when compared according to the current epsilon definitions',
                    }
    
    def __init__(self, diffInfoObject, include_basic_analysis=True) :
        """
        build our comparison statistics based on the comparison
        of two data sets
        
        the include_basic_analysis flag indicates whether the statistics generated by the
        basic_analysis method should also be generated
        """
        self.title = 'Numerical Comparison Statistics'
        
        # pull out some info we will use later
        valid_in_both        = diffInfoObject.diff_data_object.masks.valid_mask
        outside_epsilon_mask = diffInfoObject.diff_data_object.masks.outside_epsilon_mask
        mismatch_mask        = diffInfoObject.diff_data_object.masks.mismatch_mask
        aData                = diffInfoObject.a_data_object.data
        bData                = diffInfoObject.b_data_object.data
        
        assert (valid_in_both.shape        == outside_epsilon_mask.shape)
        assert (outside_epsilon_mask.shape == mismatch_mask.shape)
        assert (mismatch_mask.shape        == aData.shape)
        assert (aData.shape                == bData.shape)
        
        # fill in some simple statistics
        self.diff_outside_epsilon_count = np.sum(outside_epsilon_mask)
        self.perfect_match_count        = NumericalComparisonStatistics._get_num_perfect(aData, bData,
                                                                                         goodMask=valid_in_both)
        self.correlation                = delta.compute_correlation(aData, bData, valid_in_both)
        self.r_squared_correlation      = self.correlation * self.correlation
        self.mismatch_points_count      = np.sum(mismatch_mask)
        
        # we actually want the total number of _finite_ values rather than all the data
        total_num_finite_values = np.sum(valid_in_both)
        
        # calculate some more complex statistics
        self.mismatch_points_fraction = float(self.mismatch_points_count) / float(aData.size)
        # be careful not to divide by zero if we don't have finite data
        if total_num_finite_values > 0 :
            self.diff_outside_epsilon_fraction = float(self.diff_outside_epsilon_count) / float(total_num_finite_values)
            self.perfect_match_fraction        = float(self.perfect_match_count)        / float(total_num_finite_values)
        else:
            self.diff_outside_epsilon_fraction = 0.0
            self.perfect_match_fraction        = 0.0
        
        # if desired, do the basic analysis
        self.does_include_simple = include_basic_analysis
        if (include_basic_analysis) :
            basic_dict = NumericalComparisonStatistics.basic_analysis(diffInfoObject.diff_data_object.data,
                                                                      valid_in_both)
            if len(basic_dict) > 0 :
                self.rms_val       = basic_dict['rms_val']
                self.std_val       = basic_dict['std_val']
                self.mean_diff     = basic_dict['mean_diff']
                self.median_diff   = basic_dict['median_diff']
                self.max_diff      = basic_dict['max_diff']
                
                self.mean_delta    = basic_dict['mean_delta']
                self.median_delta  = basic_dict['median_delta']
                self.max_delta     = basic_dict['max_delta']
                self.min_delta     = basic_dict['min_delta']
            else :
                self.rms_val       = np.nan
                self.std_val       = np.nan
                self.mean_diff     = np.nan
                self.median_diff   = np.nan
                self.max_diff      = np.nan
                
                self.mean_delta    = np.nan
                self.median_delta  = np.nan
                self.max_delta     = np.nan
                self.min_delta     = np.nan
            self.temp_analysis = basic_dict
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        
        toReturn = {
                    'correlation':                   self.correlation,
                    'r-squared correlation':         self.r_squared_correlation,
                    'diff_outside_epsilon_count':    self.diff_outside_epsilon_count,
                    'diff_outside_epsilon_fraction': self.diff_outside_epsilon_fraction,
                    'perfect_match_count':           self.perfect_match_count,
                    'perfect_match_fraction':        self.perfect_match_fraction,
                    'mismatch_points_count':         self.mismatch_points_count, 
                    'mismatch_points_fraction':      self.mismatch_points_fraction
                    }
        toReturn.update(self.temp_analysis)
        
        return toReturn
    
    @staticmethod
    def doc_strings( ) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        return NumericalComparisonStatistics._doc_strings
    
    @staticmethod
    def basic_analysis(diffData, valid_mask):
        """
        do some very minimal analysis of the differences
        """
        
        # if all the data is invalid,
        # we can't do any of these forms of statistical analysis
        if np.sum(valid_mask) <= 0 :
            return { }
        
        # calculate our statistics
        root_mean_square_value = delta.calculate_root_mean_square(diffData, valid_mask)
        tempDiffData = diffData[valid_mask]
        absDiffData  = np.abs(tempDiffData)
        return {    'rms_val':       root_mean_square_value, 
                    'std_val':         np.std(tempDiffData),
                    
                    'mean_diff':       np.mean(absDiffData), 
                    'median_diff':   np.median(absDiffData),
                    'max_diff':         np.max(absDiffData),
                    
                    'mean_delta':     np.mean(tempDiffData), 
                    'median_delta': np.median(tempDiffData),
                    'max_delta':       np.max(tempDiffData),
                    'min_delta':       np.min(tempDiffData)
                    }
    
    @staticmethod
    def _get_num_perfect(aData, bData, goodMask=None):
        """
        get the number of data points where
        the value in A perfectly matches the value in B
        """
        numPerfect = 0
        if not (goodMask is None) :
            numPerfect = np.sum(aData[goodMask] == bData[goodMask])
        else :
            numPerfect = np.sum(aData == bData)
        return numPerfect

class StatisticalAnalysis (StatisticalData) :
    """
    This class represents a complete statistical analysis of two data sets.
    
    It includes the following sets of statistics:
    
    general      - a GeneralStatistics object
    comparison   - a NumericalComparisonStatistics object
    notANumber   - a NotANumberStatistics object
    missingValue - a MissingValueStatistics object
    finiteData   - a FiniteDataStatistics object
    
    It can also provide a dictionary form of the statistics or the
    documentation of the statistics.
    """
    
    def __init__ (self) :
        """
        this is a blank constructor to support our new class method creation pattern
        """
        self.title = "Statistical Summary"
    
    @classmethod
    def withSimpleData (in_class,
                        a_data,                b_data,
                        a_missing_value=None,  b_missing_value=None,
                        a_ignore_mask=None,    b_ignore_mask=None,
                        epsilon=0., epsilon_percent=None) :
        """
        do a full statistical analysis of the data, after building the data objects
        """
        
        new_object  = in_class()
        
        aDataObject = dataobj.DataObject(a_data, fillValue=a_missing_value, ignoreMask=a_ignore_mask)
        bDataObject = dataobj.DataObject(b_data, fillValue=b_missing_value, ignoreMask=b_ignore_mask)
        
        diffInfo    = dataobj.DiffInfoObject(aDataObject, bDataObject,
                                             epsilonValue=epsilon, epsilonPercent=epsilon_percent) 
        
        new_object._create_stats(diffInfo)
        
        return new_object
    
    @classmethod
    def withDataObjects (in_class,
                         a_data_object, b_data_object,
                         epsilon=0.,    epsilon_percent=None) :
        """
        do a full statistical analysis of the data, using the given data objects
        """
        
        new_object = in_class()
        
        diffInfo   = dataobj.DiffInfoObject(a_data_object, b_data_object,
                                            epsilonValue=epsilon, epsilonPercent=epsilon_percent) 
        
        new_object._create_stats(diffInfo)
        
        return new_object
    
    def _create_stats(self, diffInfoObject) :
        """
        build and set all of the statistics sets
        """
        
        self.general      = GeneralStatistics(diffInfoObject)
        self.comparison   = NumericalComparisonStatistics(diffInfoObject)
        self.notANumber   = NotANumberStatistics(diffInfoObject)
        self.missingValue = MissingValueStatistics(diffInfoObject)
        self.finiteData   = FiniteDataStatistics(diffInfoObject)
    
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        toReturn = { }
        
        # build a dictionary of all our statistics
        toReturn[self.general.title]      = self.general.dictionary_form()
        toReturn[self.comparison.title]   = self.comparison.dictionary_form()
        toReturn[self.notANumber.title]   = self.notANumber.dictionary_form()
        toReturn[self.missingValue.title] = self.missingValue.dictionary_form()
        toReturn[self.finiteData.title]   = self.finiteData.dictionary_form()
        
        return toReturn
    
    def doc_strings(self) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        return StatisticalAnalysis.doc_strings( )
    
    # TODO, use this method instead of the dictionary at the bottom of this module
    @staticmethod
    def doc_strings( ) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        toReturn = { }
        toReturn.update(GeneralStatistics.doc_strings())
        toReturn.update(NumericalComparisonStatistics.doc_strings())
        toReturn.update(NotANumberStatistics.doc_strings())
        toReturn.update(MissingValueStatistics.doc_strings())
        toReturn.update(FiniteDataStatistics.doc_strings())
        
        return toReturn

class StatisticalInspectionAnalysis (StatisticalData) :
    """
    This class represents a complete statistical analysis of a data set.
    
    It includes the following sets of statistics:
    
    general      - a GeneralInspectionStatistics object
    notANumber   - a NotANumberInspectionStatistics object
    missingValue - a MissingValueInspectionStatistics object
    finiteData   - a FiniteDataInspectionStatistics object
    
    It can also provide a dictionary form of the statistics or the
    documentation of the statistics.
    """
    
    def __init__ (self) :
        """
        this is a blank constructor to support our new class method creation pattern
        """
        self.title = "Statistical Summary"
    
    @classmethod
    def withSimpleData (in_class,
                        dataSet,
                        missingValue=None,
                        ignoreMask=None) :
        """
        do a full statistical analysis of the data, after building the data object
        """
        
        new_object  = in_class()
        
        dataObject = dataobj.DataObject(dataSet, fillValue=missingValue, ignoreMask=ignoreMask)
        dataObject.self_analysis()
        
        new_object._create_stats(dataObject)
        
        return new_object
    
    @classmethod
    def withDataObjects (in_class,
                         dataObject) :
        """
        do a full statistical analysis of the data, using the given data object
        """
        
        new_object = in_class()
        
        dataObject.self_analysis()
        new_object._create_stats(dataObject)
        
        return new_object
    
    def _create_stats(self, dataObject) :
        """
        build and set all of the statistics sets
        """
        
        self.general      = GeneralInspectionStatistics(dataObject)
        self.notANumber   = NotANumberInspectionStatistics(dataObject)
        self.missingValue = MissingValueInspectionStatistics(dataObject)
        self.finiteData   = FiniteDataInspectionStatistics(dataObject)
    
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        toReturn = { }
        
        # build a dictionary of all our statistics
        toReturn[self.general.title]      = self.general.dictionary_form()
        toReturn[self.notANumber.title]   = self.notANumber.dictionary_form()
        toReturn[self.missingValue.title] = self.missingValue.dictionary_form()
        toReturn[self.finiteData.title]   = self.finiteData.dictionary_form()
        
        return toReturn
    
    def doc_strings(self) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        return StatisticalAnalysis.doc_strings( )
    
    # TODO, use this method instead of the dictionary at the bottom of this module
    @staticmethod
    def doc_strings( ) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        toReturn = { }
        toReturn.update(GeneralInspectionStatistics.doc_strings())
        toReturn.update(NotANumberInspectionStatistics.doc_strings())
        toReturn.update(MissingValueInspectionStatistics.doc_strings())
        toReturn.update(FiniteDataInspectionStatistics.doc_strings())
        
        return toReturn

# -------------------------- documentation -----------------------------

# TODO, can this be moved?
STATISTICS_DOC_STR = '\n'.join( '%s:\n    %s' % x for x in sorted(list(StatisticalAnalysis.doc_strings().items())) ) + '\n'

if __name__=='__main__':
    import doctest
    doctest.testmod()
