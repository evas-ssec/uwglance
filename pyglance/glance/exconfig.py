#!/usr/bin/env python
# encoding: utf-8
"""
An example config file to input report generation parameters.
Please read the descriptions of the various values below, many must be
defined for your configuration file to work correctly.

The meanings of these values are described below.

Created by Eva Schiffer Jun 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

# various general settings to control how reports are created
settings = {}
# whether or not images should be generated and shown in the report
settings['shouldIncludeImages'] = True
# should we create multiple processes to make more than one image at a time?
# turning on this option can cause glance to use a very large amount of system
# resources (if your data set is particularly large, your machine may not have
# enough), but will speed up image generation in cases where your data set is
# relatively small or your machine is very powerful
settings['doFork'] = True

# the names of the latitude and longitude variables that will be used
lat_lon_info = {}
lat_lon_info['longitude'] = 'imager_prof_retr_abi_r4_generic1' # 'pixel_longitude'# the name of the longitude variable
lat_lon_info ['latitude'] = 'imager_prof_retr_abi_r4_generic2' # 'pixel_latitude' # the name of the latitude variable

"""
# the following two functions can be defined in order to filter the longitude and latitude data, for example, these could be used to compensate
# for differening data types (like ints/floats or float32/float64) or to handle slicing out only a subset of the data for analysis
# note: these two filters will only be applied to the longitude and latitude data in file A
lat_lon_info['data_filter_function_lon_in_a'] = (insert lambda function here)
lat_lon_info['data_filter_function_lat_in_a'] = (insert lambda function here)

# the following two values are optional and only need to be set if the the latitude and longitude have
# different names in file A and file B
lat_lon_info['longitude_alt_name_in_b'] = 'resampled_longitude' # the alternate name of the longitude in file B
lat_lon_info ['latitude_alt_name_in_b'] = 'resampled_latitude'  # the alternate name of the latitude in file A

# the following two functions can be defined in order to filter the longitude and latitude data, for example, these could be used to compensate
# for differening data types (like ints/floats or float32/float64) or to handle slicing out only a subset of the data for analysis
# note: these two filters will only be applied to the longitude and latitude data in file B
lat_lon_info['data_filter_function_lon_in_b'] = (insert lambda function here)
lat_lon_info['data_filter_function_lat_in_b'] = (insert lambda function here)
"""
# this value can be used to control how similar the longitude and latitude must be to be considered matching
# Note: this value is only intended to allow you to avoid very small floating point errors that would make glance
# think that your data is disparate, when really it is very close together. If you put a large epsilon in here
# the various comparison plots may contain misleading data
lat_lon_info['lon_lat_epsilon'] = 0.0001

# per variable defaults
# these default variables will only apply if you don't define them in a given variable
# description in the setOfVariables
#
# for instance, if you defined a variable "pixel_longitude" in the set of variables and
# failed to give it an 'epsilon' entry, the epsilon entry in the default values would be
# used when "pixel_longitude was analyzed
#
# if any of these defaults are not defined in your config file, glance will fall back
# on it's internal defaults.
defaultValues = {'epsilon': 0.0,                        # the acceptable difference between file A and file B
                 
                 'missing_value': -999,                 # the value to be interpreted as "missing" data
                 
                 'epsilon_failure_tolerance': None,     # the allowed fraction of epsilon comparison failure
                                                        # None indicates that variables should not be tested
                                                        # on nearness of epsilon comparison
                                                        
                 'nonfinite_data_tolerance': None       # the allowed fraction of non-finite data
                                                        # None indicates that variables should not be tested
                                                        # on amount of non-finite data
                 # the following two functions can be defined in order to filter the variable data,
                 # for example, these could be used to compensate
                 # for differening data types (like ints/floats or float32/float64)
                 # or to handle slicing out only a subset of the data for analysis
    #            'data_filter_function_a': (insert lambda function here), # note: will only be applied to file A data
    #            'data_filter_function_b': (insert lambda function here)  # note: will only be applied to file B data
                 }

# a list of all the variables to analyze, all of the details are optional,
# but at a minimum you must define variable entries in the form
#   setOfVariables['variable_name'] = {}
#
# if you chose to omit display_name, epsilon, epsilon_failure_tolerance, or
# nonfinite_data_tolerance for a given variable, the default values will
# be used for that variable
#
# If you wish to analyze all variables that match your longitude and latitude
# variables in size, do not place any entries in the setOfVariables (ie. just
# leave it as an empty dictionary)
setOfVariables = {}

setOfVariables['Total Totals'] = {           # the key of your variable entry will be used as the display name for your variable
                                  
                                  'variable_name': 'imager_prof_retr_abi_total_totals_index',       
                                                                        # this should match the variable name in your file
                                                                        
                                  'epsilon': 1.0,                       # the acceptable difference between file A and file B
                                  
                                  'missing_value': -999,                # the value to be interpreted as "missing" data
                                  
                                  'epsilon_failure_tolerance': 0.02,    # the allowed fraction of epsilon comparison failure
                                                                        # None indicates that this variable should not be tested
                                                                        # on nearness of epsilon comparison
                                                                        # note, this setting overrides the default
                                  
                                  'nonfinite_data_tolerance': None,     # the allowed fraction of non-finite data
                                                                        # None indicates that this variable should not be tested
                                                                        # on amount of non-finite data
                                                                        # note, this setting overrides the default
                                  
                                  'display_ranges':         [13.0,   14.0,  15.0,  20.0,  32.0,  40.0,  54.0,  60.0],
                                                                        # this array of ranges can be defined in order to control
                                                                        # a custom display of color ranges in any graphs produced
                                                                        # for this variable, ranges will fall between the numbers
                                                                        # in the list (inclusive of the beginning and endding numbers)
                                                                        # please list your range in increasing value if you choose
                                                                        # to use this feature
                                  
                                  'display_range_names':    [    'a',    'b',   'c',   'd',   'e',   'f',   'g'],
                                                                        # this array of labels can be defined in order to control
                                                                        # the display labels associated with your range boundaries
                                                                        # if less labels are given than boundaries, the boundaries
                                                                        # will be labeled starting with the lowest valued boundary
                                                                        # and working upwards until there are no more labels, any
                                                                        # remaining boundaries will be left unlabeled
                                                                        # TODO In the future this array may have the ability to
                                                                        # let you label your ranges (ie, the space between two
                                                                        # boundaries) or your boundaries.
                                                                        
                                  'display_colors': ('m',    'b',   'c',   'r',   'k',   'y',   'g'),
                                                                        # the colors that should be used to display the ranges
                                                                        # color definition information can be found at
                                                                        # http://matplotlib.sourceforge.net/api/colors_api.html
                                  
                                  # data filters can be defined/overridden on a variable by variable basis
#                                  'data_filter_function_a': (insert lambda function here), # note: will only be applied to file A data
#                                  'data_filter_function_b': (insert lambda function here)  # note: will only be applied to file B data
                                  
                                  }
setOfVariables['Total Precipitable Water, High'] = {
                                  'variable_name': 'imager_prof_retr_abi_total_precipitable_water_high',
                                  'epsilon': 3.0,
                                  
                                  'alternate_name_in_B': 'imager_prof_retr_abi_total_precipitable_water_low',
                                                                        # this represents an alternate name that would be
                                                                        # the equivalent variable in file B (if this is
                                                                        # defined, our primary variable name will be expected
                                                                        # to appear only in file A)
                                  
                                  'shouldIncludeImages': False          # this variable will not be plotted or compared to
                                                                        # the longitude and latitude variables
                                  }

setOfVariables['Total Precipitable Water'] = {
                                  'variable_name': 'imager_prof_retr_abi_total_precipitable_water',
                                  'epsilon': 3.0
                                  # example:
                                  # because missing, and the two tolerances are not defined here,
                                  # this variable would use the defaultValues for those 
                                  }
setOfVariables['Total Precipitable Water, Low'] = {
                                  'variable_name': 'imager_prof_retr_abi_total_precipitable_water_low',
                                  'epsilon': 3.0
                                  }
setOfVariables['Total Precipitable Water, Mid'] = {
                                  'variable_name': 'imager_prof_retr_abi_total_precipitable_water_mid',
                                  'epsilon': 3.0
                                  }
setOfVariables['Land Surface Temperature'] = {
                                  'variable_name': 'imager_prof_retr_abi_land_surface_temperature',
                                  'epsilon': 5.0
                                  }
setOfVariables['K-Index'] = {
                                  'variable_name': 'imager_prof_retr_abi_k_index',
                                  'epsilon': 2.0
                                  }
setOfVariables['Temperature Profile'] = {
                                  'variable_name': 'imager_prof_retr_abi_tprof',
                                  'epsilon': 0.1
                                  }
setOfVariables['CAPE'] = {
                                  'variable_name': 'imager_prof_retr_abi_cape',
                                  'epsilon': 1000
                                  }
setOfVariables['Lifted Index'] = {
                                  'variable_name': 'imager_prof_retr_abi_lifted_index',
                                  'epsilon': 2.0
                                  }
setOfVariables['Showalter Index'] = {
                                  'variable_name': 'imager_prof_retr_abi_showalter_index',
                                  'epsilon': 2.0
                                  }
