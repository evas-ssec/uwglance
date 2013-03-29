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
    the examination of data sets. What form of data is accepted for
    analysis is relatively abstract. 
    
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
    def doc_strings(inspect=False) :
        """
        get documentation strings that match the
        dictionary form of the statistics this class
        creates
        
        note: child classes should override this method
        """
        return { }
    
    def make_prefix_and_suffix (self, descriptionText) :
        """
        given text describing a statistic (or none)
        return an appropriate prefix and suffix
        """
        
        prefix = "" if descriptionText is None else str(descriptionText) + '_'
        suffix = "" if descriptionText is None else '_' + str(descriptionText)
        
        return prefix, suffix

class MissingValueStatistics (StatisticalData) :
    """
    A class representing information about where fill values are found
    in data. It can analyze either a pair of data sets encapsulated in a
    glance.data.DiffInfoObject or a single data set in a glance.data.DataObject.
    
    if a DiffInfoObject is given it will produce the following statistics:
    
    common_missing_count    -    count of points that are missing in both data sets
    common_missing_fraction - fraction of points that are missing in both data sets
    
    it will also include the following intermediary objects with stats about the
    individual data sets in the DiffInfoObject:
    
    a_missing_stats         - a MissingValueStatistics object specific to the a data set
    b_missing_stats         - a MissingValueStatistics object specific to the b data set
    
    when turned into a dictionary these become:
    
    a_missing_count         -    count of points that are missing in the a data set
    a_missing_fraction      - fraction of points that are missing in the a data set
    b_missing_count         -    count of points that are missing in the b data set
    b_missing_fraction      - fraction of points that are missing in the b data set
    
    if it is only given a DataObject it will produce the following :
    
    <data set descrption>missing_count    -    count of points that are missing in the data set
    <data set descrption>missing_fraction - fraction of points that are missing in the data set
    """
    
    _doc_strings = \
                    {
                    'a_missing_count':         "number of values flagged missing in A",
                    'a_missing_fraction':      "fraction of values flagged missing in A",
                    'b_missing_count':         "number of values flagged missing in B",
                    'b_missing_fraction':      "fraction of values flagged missing in B",
                    'common_missing_count':    "number of missing values in common between A and B",
                    'common_missing_fraction': "fraction of missing values in common between A and B"
                    }
    
    _doc_strings_inspection = \
                    {
                    'missing_count':         "number of values flagged missing",
                    'missing_fraction':      "fraction of values flagged missing",
                    }
    
    def __init__(self, diffInfoObject=None, dataObject=None, dataSetDescription=None) :
        """
        build our fill value related statistics
        
        diffInfoObject is assumed to be a glance.data.DiffInfoObject
        dataObject     is assumed to be a glance.data.DataObject
        
        Either the diffInfoObject or the dataObject must be passed in. If the
        diffInfoObject is passed the dataObject will be ignored and the
        a_data_object and b_data_object associated with the diffInfoObject
        will be analyzed.
        
        If only dataObject is analysed dataSetDescription will be used in labeling
        the resulting dictionary form statistics. 
        """
        self.title           = 'Missing Value Statistics'
        self.is_one_data_set = False
        
        # if we don't have comparison information and we do have a single data set
        if (diffInfoObject is None) and (dataObject is not None) :
            
            # we have one data set and should save the prefix information
            self.is_one_data_set = True
            self.desc_text     = dataSetDescription
            
            # figure out some basic statistics
            self.missing_count    = np.sum(dataObject.masks.missing_mask)
            self.missing_fraction = float(self.missing_count) / float(dataObject.data.size)
            
        # if we have a comparison object analyze the data associated with that comparison
        elif diffInfoObject is not None :
            
            # analyze each of the original data sets that are being compared
            self.a_missing_stats = MissingValueStatistics(dataObject=diffInfoObject.a_data_object, dataSetDescription="a")
            self.b_missing_stats = MissingValueStatistics(dataObject=diffInfoObject.b_data_object, dataSetDescription="b")
            
            # common statistics
            self.common_missing_count    = np.sum(diffInfoObject.a_data_object.masks.missing_mask & diffInfoObject.b_data_object.masks.missing_mask)
            self.common_missing_fraction = float(self.common_missing_count)          / float(diffInfoObject.a_data_object.data.size)
            
        else :
            raise ValueError ("No data set was given when requesting statistical analysis of missing values.")
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        
        toReturn = { }
        
        # if we only have stats for one data set
        if self.is_one_data_set :
            temp_prefix, _ = self.make_prefix_and_suffix(self.desc_text)
            toReturn = {
                        temp_prefix + 'missing_count':      self.missing_count,
                        temp_prefix + 'missing_fraction':   self.missing_fraction,
                        }
        
        # otherwise we must have stats for a comparison
        else :
            toReturn = {
                        'common_missing_count':    self.common_missing_count,
                        'common_missing_fraction': self.common_missing_fraction,
                        }
            a_dict = self.a_missing_stats.dictionary_form()
            toReturn.update(a_dict)
            b_dict = self.b_missing_stats.dictionary_form()
            toReturn.update(b_dict)
        
        return toReturn
    
    @staticmethod
    def doc_strings(inspect=False) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        return MissingValueStatistics._doc_strings if not inspect else MissingValueStatistics._doc_strings_inspection

class FiniteDataStatistics (StatisticalData) :
    """
    A class representing information about where finite values are found
    in data. It can analyze either a pair of data sets encapsulated in a
    glance.data.DiffInfoObject or a single data set in a glance.data.DataObject.
    
    when a single data set is analyzed the following stats are produced:
    
    <data prefix>finite_count    - the number   of finite data values in the data set
    <data prefix>finite_fraction - the fraction of finite data values in the data set
    
    if a DiffInfoObject is given for analysis the following statistics are produced:
    
    common_finite_count         - the number   of finite values the two data sets have in common
    common_finite_fraction      - the fraction of finite values the two data sets have in common
    finite_in_only_one_count    - the number   of points that are finite in only one of the two sets
    finite_in_only_one_fraction - the fraction of points that are finite in only one of the two sets
    
    it will also include the following intermediary objects with stats about the
    individual data sets in the DiffInfoObject:
    
    a_finite_stats               - a FiniteDataStatistics object with further stats on the a data set
    b_finite_stats               - a FiniteDataStatistics object with further stats on the b data set
    
    and the dictionary form will includes the following statistics:
    
    a_finite_count              - the number   of finite data values in the a data set
    a_finite_fraction           - the fraction of finite data values in the a data set
    b_finite_count              - the number   of finite data values in the b data set
    b_finite_fraction           - the fraction of finite data values in the b data set
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
    
    _doc_strings_inspection = \
                    {
                    'finite_count': "number of finite values",
                    'finite_fraction': "fraction of finite values (out of all data points in set)",
                    }
    
    def __init__(self, diffInfoObject=None, dataObject=None, dataSetDescription=None) :
        """
        build our finite data related statistics 
        
        diffInfoObject is assumed to be a glance.data.DiffInfoObject
        dataObject     is assumed to be a glance.data.DataObject
        
        Either the diffInfoObject or the dataObject must be passed in. If the
        diffInfoObject is passed the dataObject will be ignored and the
        a_data_object and b_data_object associated with the diffInfoObject
        will be analyzed.
        
        If only dataObject is analysed dataSetDescription will be used in labeling
        the resulting dictionary form statistics. 
        """
        self.title           = 'Finite Data Statistics'
        self.is_one_data_set = False
        
        # if we don't have comparison information and we do have a single data set
        if (diffInfoObject is None) and (dataObject is not None) :
            
            # we have one data set and should save the prefix information
            self.is_one_data_set = True
            self.desc_text       = dataSetDescription
            
            # figure out some basic statistics
            self.finite_count    = np.sum(dataObject.masks.valid_mask)
            self.finite_fraction = float(self.finite_count) / float(dataObject.data.size)
            
        # if we have a comparison object analyze the data associated with that comparison
        elif diffInfoObject is not None :
            
            # analyze each of the original data sets that are being compared
            self.a_finite_stats = FiniteDataStatistics(dataObject=diffInfoObject.a_data_object, dataSetDescription="a")
            self.b_finite_stats = FiniteDataStatistics(dataObject=diffInfoObject.b_data_object, dataSetDescription="b")
            
            # calculate some common statistics
            self.common_finite_count = np.sum(diffInfoObject.a_data_object.masks.valid_mask & diffInfoObject.b_data_object.masks.valid_mask)
            # use an exclusive or to check which points are finite in only one of the two data sets
            self.finite_in_only_one_count = np.sum((diffInfoObject.a_data_object.masks.valid_mask ^ diffInfoObject.b_data_object.masks.valid_mask) \
                                                    & ~diffInfoObject.diff_data_object.masks.ignore_mask)
            self.common_finite_fraction      = float(self.common_finite_count)      / float(diffInfoObject.a_data_object.data.size)
            self.finite_in_only_one_fraction = float(self.finite_in_only_one_count) / float(diffInfoObject.a_data_object.data.size)
            
        else:
            raise ValueError ("No data set was given when requesting statistical analysis of finite values.")
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        
        toReturn = { }
        
        # if we only have stats for one data set
        if self.is_one_data_set :
            temp_prefix, _ = self.make_prefix_and_suffix(self.desc_text)
            toReturn = {
                        temp_prefix + 'finite_count':    self.finite_count,
                        temp_prefix + 'finite_fraction': self.finite_fraction,
                        }
        
        # otherwise we must have stats for a comparison
        else :
            toReturn = {
                        'common_finite_count':         self.common_finite_count,
                        'common_finite_fraction':      self.common_finite_fraction,
                        'finite_in_only_one_count':    self.finite_in_only_one_count,
                        'finite_in_only_one_fraction': self.finite_in_only_one_fraction,
                        }
            a_dict = self.a_finite_stats.dictionary_form()
            toReturn.update(a_dict)
            b_dict = self.b_finite_stats.dictionary_form()
            toReturn.update(b_dict)
        
        return toReturn
    
    @staticmethod
    def doc_strings(inspect=False) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        return FiniteDataStatistics._doc_strings if not inspect else FiniteDataStatistics._doc_strings_inspection

class NotANumberStatistics (StatisticalData) :
    """
    A class representing information about where non-numerical values are found
    in data. It can analyze either a pair of data sets encapsulated in a
    glance.data.DiffInfoObject or a single data set in a glance.data.DataObject.
    
    when a single data set is analyzed the following stats are produced:
    
    nan_count         - the number   of non finite values that are present in the data set
    nan_fraction      - the fraction of non finite values that are present in the data set
    
    if a DiffInfoObject is given for analysis the following statistics are produced:
    
    common_nan_count    - the number   of non finite values that are shared between the data sets
    common_nan_fraction - the fraction of non finite values that are shared between the data sets
    
    if a DiffInfoObject is given the object will also have:
    
    a_finite_stats               - a NotANumberStatistics object with further stats on the a data set
    b_finite_stats               - a NotANumberStatistics object with further stats on the b data set
    
    and the dictionary form will includes the following statistics:
    
    a_nan_count         - the number   of non finite values that are present in the a data set
    a_nan_fraction      - the fraction of non finite values that are present in the a data set
    b_nan_count         - the number   of non finite values that are present in the b data set
    b_nan_fraction      - the fraction of non finite values that are present in the b data set
    """
    
    _doc_strings = {
                    'a_nan_count': "number of NaNs in A",
                    'a_nan_fraction': "fraction of NaNs in A",
                    'b_nan_count': "number of NaNs in B",
                    'b_nan_fraction': "fraction of NaNs in B",
                    'common_nan_count': "number of NaNs in common between A and B",
                    'common_nan_fraction': "fraction of NaNs in common between A and B"
                    }
    
    _doc_strings_inspection = \
                    {
                    'nan_count': "number of NaNs",
                    'nan_fraction': "fraction of NaNs",
                    }
    
    def __init__(self, diffInfoObject=None, dataObject=None, dataSetDescription=None) :
        """
        build our nonfinite data related statistics
        
        diffInfoObject is assumed to be a glance.data.DiffInfoObject
        dataObject     is assumed to be a glance.data.DataObject
        
        Either the diffInfoObject or the dataObject must be passed in. If the
        diffInfoObject is passed the dataObject will be ignored and the
        a_data_object and b_data_object associated with the diffInfoObject
        will be analyzed.
        
        If only dataObject is analysed dataSetDescription will be used in labeling
        the resulting dictionary form statistics. 
        """
        self.title           = 'NaN Statistics'
        self.is_one_data_set = False
        
        # if we don't have comparison information and we do have a single data set
        if (diffInfoObject is None) and (dataObject is not None) :
            
            # we have one data set and should save the prefix information
            self.is_one_data_set = True
            self.desc_text       = dataSetDescription
            
            # get some basic statistics
            self.nan_count = np.sum(dataObject.masks.non_finite_mask)
            self.nan_fraction = float(self.nan_count) / float(dataObject.data.size)
            
        # if we have a comparison object analyze the data associated with that comparison
        elif diffInfoObject is not None :
            
            # analyze each of the original data sets that are being compared
            self.a_nan_stats = NotANumberStatistics(dataObject=diffInfoObject.a_data_object, dataSetDescription="a")
            self.b_nan_stats = NotANumberStatistics(dataObject=diffInfoObject.b_data_object, dataSetDescription="b")
            
            # calculate some common statistics
            self.common_nan_count = np.sum(diffInfoObject.a_data_object.masks.non_finite_mask & diffInfoObject.b_data_object.masks.non_finite_mask)
            self.common_nan_fraction = float(self.common_nan_count) / float(diffInfoObject.a_data_object.data.size)
            
        else:
            raise ValueError ("No data set was given when requesting statistical analysis of NaN values.")
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        
        toReturn = { }
        
        # if we only have stats for one data set
        if self.is_one_data_set :
            temp_prefix, _ = self.make_prefix_and_suffix(self.desc_text)
            toReturn = {
                        temp_prefix + 'nan_count':         self.nan_count,
                        temp_prefix + 'nan_fraction':      self.nan_fraction,
                        }
        
        # otherwise we must have stats for a comparison
        else :
            toReturn = {
                        'common_nan_count':    self.common_nan_count,
                        'common_nan_fraction': self.common_nan_fraction
                        }
            a_dict = self.a_nan_stats.dictionary_form()
            toReturn.update(a_dict)
            b_dict = self.b_nan_stats.dictionary_form()
            toReturn.update(b_dict)
        
        return toReturn
    
    @staticmethod
    def doc_strings(inspect=False) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        return NotANumberStatistics._doc_strings if not inspect else NotANumberStatistics._doc_strings_inspection

class GeneralStatistics (StatisticalData) :
    """
    A class representing general information about data. It can analyze either a
    pair of data sets encapsulated in a glance.data.DiffInfoObject or a single
    data set in a glance.data.DataObject.
    
    if a single DataObject is given the following will be produced:
    (some of these are labeled with any dataSetDescription given in the
    constructor)
    
    missing_value                 - the fill data value
    max                           - the maximum value
    min                           - the minimum value
    num_data_points               - the total number of data points
    shape                         - the shape of the data
    spatially_invalid_pts_ignored - number of points corresponding to invalid lat/lon in the set
                                    (optional if no /lon lat mapped)
    mean                          - the mean of the data values
    median                        - the median of the data values
    std_val                       - the standard deviation of the data values
    
    if a DiffInfoObject is given these comparison stats will be produced:
    
    epsilon         - the fixed epsilon value
    epsilon_percent - the percentage of the a set that will be used for comparison
    num_data_points - the number of data points in each of the sets
    shape           - the shape of each of the data sets
    
    it will also have the following self owned variables:
    
    a_gen_stats      - a GeneralStatistics object with further stats on the a data set
    b_gen_stats      - a GeneralStatistics object with further stats on the b data set
    
    in dictionary form those objects will produce:
    
    a_missing_value - the fill data value in the a set
    b_missing_value - the fill data value in the b set
    max_a           - the maximum value in the a set
    max_b           - the maximum value in the b set
    min_a           - the minimum value in the a set
    min_b           - the minimum value in the b set
    """
    
    _doc_strings = {
                    'a_missing_value': 'the value that is considered \"missing\" or \"fill\" data when it is found in A',
                    'b_missing_value': 'the value that is considered \"missing\" or \"fill\" data when it is found in B',
                    'epsilon': 'amount of difference between matching data points in A and B that is considered acceptable',
                    'epsilon_percent': 'the percentage of difference (of A\'s value) that is acceptable between A and B (optional)',
                    'max_a': 'the maximum finite, non-missing value found in A',
                    'max_b': 'the maximum finite, non-missing value found in B',
                    'min_a': 'the minimum finite, non-missing value found in A',
                    'min_b': 'the minimum finite, non-missing value found in B',
                    'num_data_points': "number of data values in A",
                    'shape': "shape of A",
                    'spatially_invalid_pts_ignored_a': 'number of points with invalid latitude/longitude information in A that were' +
                                                            ' ignored for the purposes of data analysis and presentation',
                    'spatially_invalid_pts_ignored_b': 'number of points with invalid latitude/longitude information in B that were' +
                                                            ' ignored for the purposes of data analysis and presentation',
                    # these are new!
                    'mean_a': "the mean of all finite, non-missing values found in A",
                    'mean_b': "the mean of all finite, non-missing values found in B",
                    'median_a': "the median of all finite, non-missing values in A",
                    'median_b': "the median of all finite, non-missing values in B",
                    'std_val_a': "the standard deviation of all finite, non-missing values in A",
                    'std_val_b': "the standard deviation of all finite, non-missing values in B",
                    }
    
    _doc_strings_inspect = \
                    {
                    'missing_value': 'the value that is considered \"missing\" or \"fill\" data in this data set',
                    'max': 'the maximum finite, non-missing value found in the data',
                    'min': 'the minimum finite, non-missing value found in the data',
                    'num_data_points': "number of data points (may be valid or invalid data)",
                    'shape': "shape of the data",
                    'spatially_invalid_pts_ignored': 'number of points with invalid latitude/longitude information ' +
                                                     'in the data that were' +
                                                     ' ignored for the purposes of data analysis and presentation',
                    'mean': "the mean of all finite, non-missing values in the data",
                    'median': "the median of all finite, non-missing values in the data",
                    'std_val': "the standard deviation of all finite, non-missing values in the data",
                    }
    
    def __init__(self, diffInfoObject=None, dataObject=None,
                 doExtras=False, dataSetDescription=None) :
        """
        build our general statistics based on the comparison of two data sets
        
        diffInfoObject is assumed to be a glance.data.DiffInfoObject
        dataObject     is assumed to be a glance.data.DataObject
        
        Either the diffInfoObject or the dataObject must be passed in. If the
        diffInfoObject is passed the dataObject will be ignored and the
        a_data_object and b_data_object associated with the diffInfoObject
        will be analyzed.
        
        If only dataObject is analyzed dataSetDescription will be
        used in labeling the resulting dictionary form statistics.
        
        If you are passing a single dataObject and would like shape and size
        statistics reported as well, pass doExtras as True (otherwise these
        stats will be omitted).
        """
        self.title           = 'General Statistics'
        self.is_one_data_set = False
        
        # if we don't have comparison information and we do have a single data set
        if (diffInfoObject is None) and (dataObject is not None) :
            
            # we have one data set and should save the prefix/suffix information
            self.is_one_data_set = True
            self.do_extras       = doExtras
            self.desc_text       = dataSetDescription
            
            # grab the valid data for some calculations
            tempGoodData = dataObject.data[dataObject.masks.valid_mask]
            
            # fill in our statistics
            self.missing_value   = dataObject.select_fill_value()
            self.max             =    np.max(tempGoodData)
            self.min             =    np.min(tempGoodData)
            self.mean            =   np.mean(tempGoodData)
            self.median          = np.median(tempGoodData)
            self.std_val         =    np.std(tempGoodData)
            # also calculate the invalid points
            self.spatially_invalid_pts_ignored = np.sum(dataObject.masks.ignore_mask)
            
            # if we should also do extra stats, do so
            if (doExtras) :
                self.num_data_points = dataObject.masks.missing_mask.size
                self.shape           = dataObject.masks.missing_mask.shape
            
        # if we have a comparison object analyze the data associated with that comparison
        elif diffInfoObject is not None :
            
            # analyze each of the original data sets that are being compared
            self.a_gen_stats = GeneralStatistics(dataObject=diffInfoObject.a_data_object, dataSetDescription="a")
            self.b_gen_stats = GeneralStatistics(dataObject=diffInfoObject.b_data_object, dataSetDescription="b")
            
            # fill in our statistics
            self.epsilon         = diffInfoObject.epsilon_value
            self.epsilon_percent = diffInfoObject.epsilon_percent
            self.num_data_points = diffInfoObject.a_data_object.masks.missing_mask.size
            self.shape           = diffInfoObject.a_data_object.masks.missing_mask.shape
            # also calculate the invalid points
            self.spatially_invalid_pts_ignored_in_a = np.sum(diffInfoObject.a_data_object.masks.ignore_mask)
            self.spatially_invalid_pts_ignored_in_b = np.sum(diffInfoObject.b_data_object.masks.ignore_mask)
            
        else:
            raise ValueError ("No data set was given when requesting general statistical analysis.")
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        
        toReturn = { }
        
        # if we only have stats for one data set
        if self.is_one_data_set :
            temp_prefix, temp_suffix = self.make_prefix_and_suffix(self.desc_text)
            toReturn = {
                        temp_prefix + 'missing_value':                 self.missing_value,
                        'max'                           + temp_suffix: self.max,
                        'min'                           + temp_suffix: self.min,
                        'mean'                          + temp_suffix: self.mean,
                        'median'                        + temp_suffix: self.median,
                        'std_val'                       + temp_suffix: self.std_val,
                        'spatially_invalid_pts_ignored' + temp_suffix: self.spatially_invalid_pts_ignored,
                        }
            
            if self.do_extras :
                toReturn['num_data_points'] = self.num_data_points
                toReturn['shape']           = self.shape
        
        # otherwise we must have stats for a comparison
        else :
            toReturn = {
                        'epsilon':         self.epsilon,
                        'epsilon_percent': self.epsilon_percent,
                        'num_data_points': self.num_data_points,
                        'shape':           self.shape,
                        }
            a_dict = self.a_gen_stats.dictionary_form()
            toReturn.update(a_dict)
            b_dict = self.b_gen_stats.dictionary_form()
            toReturn.update(b_dict)
        
        return toReturn
    
    @staticmethod
    def doc_strings(inspect=False) :
        """
        get documentation strings that match the
        dictionary form of the statistics
        """
        
        return GeneralStatistics._doc_strings if not inspect else GeneralStatistics._doc_strings_inspect

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
        valid_in_both           = diffInfoObject.diff_data_object.masks.valid_mask
        aData                   = diffInfoObject.a_data_object.data
        bData                   = diffInfoObject.b_data_object.data
        total_num_finite_values = np.sum(valid_in_both) # just the finite values, not all data
        
        # fill in some simple statistics
        self.diff_outside_epsilon_count = np.sum(diffInfoObject.diff_data_object.masks.outside_epsilon_mask)
        self.perfect_match_count        = NumericalComparisonStatistics._get_num_perfect(aData, bData,
                                                                                         goodMask=valid_in_both)
        self.correlation                = delta.compute_correlation(aData, bData, valid_in_both)
        self.r_squared_correlation      = self.correlation * self.correlation
        self.mismatch_points_count      = np.sum(diffInfoObject.diff_data_object.masks.mismatch_mask)
        
        # calculate some more complex statistics, be careful not to divide by zero
        self.mismatch_points_fraction      = float(self.mismatch_points_count)      / float(aData.size)              if (aData.size > 0)              else 0.0
        self.diff_outside_epsilon_fraction = float(self.diff_outside_epsilon_count) / float(total_num_finite_values) if (total_num_finite_values > 0) else 0.0
        self.perfect_match_fraction        = float(self.perfect_match_count)        / float(total_num_finite_values) if (total_num_finite_values > 0) else 0.0
        
        # if desired, do the basic analysis
        self.temp_analysis = NumericalComparisonStatistics.basic_analysis(diffInfoObject.diff_data_object.data, valid_in_both) if include_basic_analysis else { }
        self.rms_val       = self.temp_analysis['rms_val']      if (len(self.temp_analysis) > 0) else np.nan
        self.std_val       = self.temp_analysis['std_val']      if (len(self.temp_analysis) > 0) else np.nan
        self.mean_diff     = self.temp_analysis['mean_diff']    if (len(self.temp_analysis) > 0) else np.nan
        self.median_diff   = self.temp_analysis['median_diff']  if (len(self.temp_analysis) > 0) else np.nan
        self.max_diff      = self.temp_analysis['max_diff']     if (len(self.temp_analysis) > 0) else np.nan
        self.mean_delta    = self.temp_analysis['mean_delta']   if (len(self.temp_analysis) > 0) else np.nan
        self.median_delta  = self.temp_analysis['median_delta'] if (len(self.temp_analysis) > 0) else np.nan
        self.max_delta     = self.temp_analysis['max_delta']    if (len(self.temp_analysis) > 0) else np.nan
        self.min_delta     = self.temp_analysis['min_delta']    if (len(self.temp_analysis) > 0) else np.nan
    
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
        
        # if everything's invalid, stop now
        if np.sum(valid_mask) <= 0 :
            return { }
        
        # calculate and return statistics
        root_mean_square_value = delta.calculate_root_mean_square(diffData, valid_mask)
        tempDiffData           = diffData[valid_mask]
        absDiffData            = np.abs(tempDiffData)
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
        
        if goodMask is None :
            numPerfect = np.sum(aData == bData)
        else :
            numPerfect = np.sum(aData[goodMask] == bData[goodMask])
        
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
    
    It can also provide a dictionary form of the statistics and
    documentation for the statistics.
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
        
        self.general      = GeneralStatistics            (diffInfoObject=diffInfoObject)
        self.comparison   = NumericalComparisonStatistics(diffInfoObject)
        self.notANumber   = NotANumberStatistics         (diffInfoObject=diffInfoObject)
        self.missingValue = MissingValueStatistics       (diffInfoObject=diffInfoObject)
        self.finiteData   = FiniteDataStatistics         (diffInfoObject=diffInfoObject)
    
    def check_pass_or_fail(self,
                           epsilon_failure_tolerance   =np.nan, epsilon_failure_tolerance_default   =None,
                           non_finite_data_tolerance   =np.nan, non_finite_data_tolerance_default   =None,
                           total_data_failure_tolerance=np.nan, total_data_failure_tolerance_default=None,
                           min_acceptable_r_squared    =np.nan, min_acceptable_r_squared_default    =None
                           ) :
        """
        Check whether the variable passed analysis, failed analysis, or
        did not need to be quantitatively tested
        
        also returns information about the fractions of failure
        """
        
        passValues = [ ]
        
        # test the epsilon value tolerance
        
        # get the tolerance for failures compared to epsilon
        epsilonTolerance = epsilon_failure_tolerance if epsilon_failure_tolerance is not np.nan else epsilon_failure_tolerance_default
        
        # did we fail based on the epsilon?
        failed_fraction = self.comparison.diff_outside_epsilon_fraction
        passed_epsilon  = None if (epsilonTolerance is None) else (failed_fraction <= epsilonTolerance)
        passValues.append(passed_epsilon)
        
        # test the nonfinite tolerance
        
        # get the tolerance for failures in amount of nonfinite data (in spatially valid areas)
        nonfiniteTolerance = non_finite_data_tolerance if non_finite_data_tolerance is not np.nan else non_finite_data_tolerance_default
        
        # did we fail based on nonfinite data
        non_finite_diff_fraction = self.finiteData.finite_in_only_one_fraction
        passed_nonfinite         = None if (nonfiniteTolerance is None) else (non_finite_diff_fraction <= nonfiniteTolerance)
        passValues.append(passed_nonfinite)
        
        # test if the total failed percentage is acceptable
        
        # get the total percentage of failed data that is acceptable
        totalFailTolerance = total_data_failure_tolerance if total_data_failure_tolerance is not np.nan else total_data_failure_tolerance_default
        
        # did we fail based on all data failures?
        passed_all_percentage = None if (totalFailTolerance is None) else ((non_finite_diff_fraction + failed_fraction) <= totalFailTolerance)
        passValues.append(passed_all_percentage)
        
        # test the r-squared correlation coefficent
        
        # get the minimum acceptable r-squared correlation coefficient
        min_r_squared = min_acceptable_r_squared if (min_acceptable_r_squared is not np.nan) else min_acceptable_r_squared_default
        
        # did we fail based on the r-squared correlation coefficient?
        r_squared_value  = None if (min_r_squared is None) else self.comparison.r_squared_correlation
        passed_r_squared = None if (min_r_squared is None) else (r_squared_value >= min_r_squared)
        passValues.append(passed_r_squared)
        
        # figure out the overall pass/fail result
        didPass = None
        for passValue in passValues :
            # if passValue isn't none, we need to update didPass
            if passValue is not None :
                if didPass is not None :
                    didPass = passValue and didPass
                else :
                    didPass = passValue
        
        return didPass, failed_fraction, non_finite_diff_fraction, r_squared_value
    
    def dictionary_form(self) :
        """
        get a dictionary form of the statistics
        """
        toReturn = { }
        
        # build a dictionary of all our statistics
        toReturn[self.general.title]      =      self.general.dictionary_form()
        toReturn[self.comparison.title]   =   self.comparison.dictionary_form()
        toReturn[self.notANumber.title]   =   self.notANumber.dictionary_form()
        toReturn[self.missingValue.title] = self.missingValue.dictionary_form()
        toReturn[self.finiteData.title]   =   self.finiteData.dictionary_form()
        
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
        toReturn.update(            GeneralStatistics.doc_strings())
        toReturn.update(NumericalComparisonStatistics.doc_strings())
        toReturn.update(         NotANumberStatistics.doc_strings())
        toReturn.update(       MissingValueStatistics.doc_strings())
        toReturn.update(         FiniteDataStatistics.doc_strings())
        
        return toReturn

class StatisticalInspectionAnalysis (StatisticalData) :
    """
    This class represents a complete statistical analysis of a data set.
    
    It includes the following sets of statistics:
    
    general      - a GeneralStatistics object
    notANumber   - a NotANumberStatistics object
    missingValue - a MissingValueStatistics object
    finiteData   - a FiniteDataStatistics object
    
    It can also provide a dictionary form of the statistics and
    documentation for the statistics.
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
        
        self.general      = GeneralStatistics(     dataObject=dataObject,
                                                           doExtras=True)
        self.notANumber   = NotANumberStatistics(  dataObject=dataObject)
        self.missingValue = MissingValueStatistics(dataObject=dataObject)
        self.finiteData   = FiniteDataStatistics(  dataObject=dataObject)
    
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
        toReturn.update(     GeneralStatistics.doc_strings(inspect=True))
        toReturn.update(  NotANumberStatistics.doc_strings(inspect=True))
        toReturn.update(MissingValueStatistics.doc_strings(inspect=True))
        toReturn.update(  FiniteDataStatistics.doc_strings(inspect=True))
        
        return toReturn

# -------------------------- documentation -----------------------------

# TODO, can this be moved?
STATISTICS_DOC_STR      = '\n'.join( '%s:\n    %s' % x for x in sorted(list(          StatisticalAnalysis.doc_strings().items())) ) + '\n'
INSP_STATISTICS_DOC_STR = '\n'.join( '%s:\n    %s' % x for x in sorted(list(StatisticalInspectionAnalysis.doc_strings().items())) ) + '\n'

if __name__=='__main__':
    import doctest
    doctest.testmod()
