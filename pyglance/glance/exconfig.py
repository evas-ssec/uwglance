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

import glance.constants as constants

# various general settings to control how reports are created
settings = {}
# whether or not images should be generated and shown in the report
settings[constants.DO_MAKE_IMAGES_KEY] = True
# whether or not separate threads should be spawned to create each image
# for the purpose of controlling python's memory usage.
# this feature should not be used in Mac OSX owing to a bug in multithreading
# but at appears to work well in other Unix systems
settings[constants.DO_CLEAR_MEM_THREADED_KEY] = False
# should we create multiple processes to make more than one image at a time?
# turning on this option can cause glance to use a very large amount of system
# resources (if your data set is particularly large, your machine may not have
# enough), but will speed up image generation in cases where your data set is
# relatively small or your machine is very powerful
settings[constants.DO_MAKE_FORKS_KEY] = False
# should the two original data sets for a variable be plotted in the same range?
# by default each data set will be plotted in it's own range, if you set this
# value to True, then the maximum of the two ranges will be used to plot both
settings[constants.USE_SHARED_ORIG_RANGE_KEY] = False

# the names of the latitude and longitude variables that will be used
lat_lon_info = {}
lat_lon_info[constants.LONGITUDE_NAME_KEY] = 'imager_prof_retr_abi_r4_generic1' # 'pixel_longitude'# the name of the longitude variable
lat_lon_info [constants.LATITUDE_NAME_KEY] = 'imager_prof_retr_abi_r4_generic2' # 'pixel_latitude' # the name of the latitude variable

"""
# the following two functions can be defined in order to filter the longitude and latitude data, for example, these could be used to compensate
# for differening data types (like ints/floats or float32/float64) or to handle slicing out only a subset of the data for analysis
# note: these two filters will only be applied to the longitude and latitude data in file A
lat_lon_info[constants.LON_FILTER_FUNCTION_A_KEY] = (insert lambda function here)
lat_lon_info[constants.LAT_FILTER_FUNCTION_A_KEY] = (insert lambda function here)

# the following two values are optional and only need to be set if the the latitude and longitude have
# different names in file A and file B
lat_lon_info[constants.LON_ALT_NAME_IN_B_KEY] = 'resampled_longitude' # the alternate name of the longitude in file B
lat_lon_info[constants.LAT_ALT_NAME_IN_B_KEY] = 'resampled_latitude'  # the alternate name of the latitude in file A

# the following two functions can be defined in order to filter the longitude and latitude data, for example, these could be used to compensate
# for differening data types (like ints/floats or float32/float64) or to handle slicing out only a subset of the data for analysis
# note: these two filters will only be applied to the longitude and latitude data in file B
lat_lon_info[constants.LON_FILTER_FUNCTION_B_KEY] = (insert lambda function here)
lat_lon_info[constants.LAT_FILTER_FUNCTION_B_KEY] = (insert lambda function here)
"""

"""
# if you wish to load the longitude and latitude from a file, use these values to specify what file
# the other longitude and latitude values specifying formatting functions, variable names, etc. will be
# used on the specified file.
lat_lon_info[constants.LONLAT_ALT_FILE_A_KEY] = '/path/to/alternate/file/for/lonlat/to/use/in/a'
lat_lon_info[constants.LONLAT_ALT_FILE_B_KEY] = '/path/to/alternate/file/for/lonlat/to/use/in/b'
"""

# this value can be used to control how similar the longitude and latitude must be to be considered matching
# Note: this value is only intended to allow you to avoid very small floating point errors that would make glance
# think that your data is disparate, when really it is very close together. If you put a large epsilon in here
# the various comparison plots may contain misleading data
lat_lon_info[constants.LON_LAT_EPSILON_KEY] = 0.0001

# per variable defaults
# these default variables will only apply if you don't define them in a given variable
# description in the setOfVariables
#
# for instance, if you defined a variable "pixel_longitude" in the set of variables
# and failed to give it a constants.EPSILON_KEY entry, the epsilon entry in the default
# values would be used when "pixel_longitude" was analyzed
#
# if any of these defaults are not defined in your config file, glance will fall back
# on it's internal defaults.
defaultValues = {constants.EPSILON_KEY: 0.0,            # the acceptable difference between file A and file B
                 
                 constants.FILL_VALUE_KEY: -999,        # the value to be interpreted as fill data
                 
                 constants.EPSILON_FAIL_TOLERANCE_KEY: None,
                                                        # the allowed fraction of epsilon comparison failure
                                                        # None indicates that variables should not be tested
                                                        # on nearness of epsilon comparison
                                                        
                 constants.NONFINITE_TOLERANCE_KEY: None,
                                                        # the allowed fraction of non-finite data that
                                                        # that differs between the two files
                                                        # None indicates that variables should not be tested
                                                        # on amount of non-finite data
                 
                 constants.DO_IMAGES_ONLY_ON_FAIL_KEY: True
                                                        # only create the variable images if the variable
                                                        # fails it's tolerance tests
                                                        # Note: This setting can be overriden if you set the 
                                                        # perVariable version of DO_MAKE_IMAGES_KEY to false
                                                        # (then plots will never be created for that variable)
                 
                 # the following two functions can be defined in order to filter the variable data,
                 # for example, these could be used to compensate
                 # for differening data types (like ints/floats or float32/float64)
                 # or to handle slicing out only a subset of the data for analysis
    #            constants.FILTER_FUNCTION_A_KEY: (insert lambda function here), # note: will only be applied to file A data
    #            constants.FILTER_FUNCTION_B_KEY: (insert lambda function here)  # note: will only be applied to file B data
                 }

# a list of all the variables to analyze, all of the details are optional,
# but at a minimum you must define variable entries in the form
#   setOfVariables["display name for the variable"] =   {
#                                                        constants.VARIABLE_TECH_NAME_KEY: "name of the variable in the file",
#                                                       }
#
# if you chose to omit constants.EPSILON_KEY,
# constants.EPSILON_FAIL_TOLERANCE_KEY, or constants.NONFINITE_TOLERANCE_KEY
# for a given variable, the default values will be used for that variable
#
# If you wish to analyze all variables that match your longitude and latitude
# variables in size, do not place any entries in the setOfVariables (ie. just
# leave it as an empty dictionary)
setOfVariables = { }

setOfVariables['Total Totals'] = {           # the key of your variable entry will be used as the display name for your variable
                                  
                                  constants.VARIABLE_TECH_NAME_KEY: 'imager_prof_retr_abi_total_totals_index',       
                                                                        # this should match the variable name in your file
                                                                        
                                  constants.EPSILON_KEY: 1.0,           # the acceptable difference between file A and file B
                                  
                                  constants.FILL_VALUE_KEY: -999,       # the value to be interpreted as "missing" data
                                  
                                  constants.EPSILON_FAIL_TOLERANCE_KEY: 0.02,
                                                                        # the allowed fraction of epsilon comparison failure
                                                                        # None indicates that this variable should not be tested
                                                                        # on nearness of epsilon comparison
                                                                        # note, this setting overrides the default
                                  
                                  constants.NONFINITE_TOLERANCE_KEY: None,
                                                                        # the allowed fraction of non-finite data
                                                                        # None indicates that this variable should not be tested
                                                                        # on amount of non-finite data
                                                                        # note, this setting overrides the default
                                  
                                  constants.DISPLAY_RANGES_KEY:         [13.0,   14.0,  15.0,  20.0,  32.0,  40.0,  54.0,  60.0],
                                                                        # this array of ranges can be defined in order to control
                                                                        # a custom display of color ranges in any graphs produced
                                                                        # for this variable, ranges will fall between the numbers
                                                                        # in the list (inclusive of the beginning and endding numbers)
                                                                        # please list your range in increasing value if you choose
                                                                        # to use this feature
                                  
                                  constants.DISPLAY_RANGE_NAMES_KEY:    [    'a',    'b',   'c',   'd',   'e',   'f',   'g'],
                                                                        # this array of labels can be defined in order to control
                                                                        # the display labels associated with your range boundaries
                                                                        # if less labels are given than boundaries, the boundaries
                                                                        # will be labeled starting with the lowest valued boundary
                                                                        # and working upwards until there are no more labels, any
                                                                        # remaining boundaries will be left unlabeled
                                                                        # TODO In the future this array may have the ability to
                                                                        # let you label your ranges (ie, the space between two
                                                                        # boundaries) or your boundaries.
                                                                        
                                  constants.DISPLAY_RANGE_COLORS_KEY: ('m',    'b',   'c',   'r',   'k',   'y',   'g'),
                                                                        # the colors that should be used to display the ranges
                                                                        # color definition information can be found at
                                                                        # http://matplotlib.sourceforge.net/api/colors_api.html
                                  
                                  # data filters can be defined/overridden on a variable by variable basis
#                                  constants.FILTER_FUNCTION_A_KEY: (insert lambda function here), # note: will only be applied to file A data
#                                  constants.FILTER_FUNCTION_B_KEY: (insert lambda function here)  # note: will only be applied to file B data
                                  
                                  }
setOfVariables['Total Precipitable Water, High'] = {
                                  constants.VARIABLE_TECH_NAME_KEY: 'imager_prof_retr_abi_total_precipitable_water_high',
                                  constants.EPSILON_KEY: 3.0,
                                  
                                  constants.VARIABLE_B_TECH_NAME_KEY: 'imager_prof_retr_abi_total_precipitable_water_low',
                                                                        # this represents an alternate name that would be
                                                                        # the equivalent variable in file B (if this is
                                                                        # defined, our primary variable name will be expected
                                                                        # to appear only in file A)
                                  
                                  constants.DO_MAKE_IMAGES_KEY: False   # this variable will not be plotted or compared to
                                                                        # the longitude and latitude variables
                                  }

setOfVariables['Total Precipitable Water'] = {
                                  constants.VARIABLE_TECH_NAME_KEY: 'imager_prof_retr_abi_total_precipitable_water',
                                  constants.EPSILON_KEY: 3.0
                                  # example:
                                  # because the fill value, and the two tolerances are not defined here,
                                  # this variable would use the defaultValues for those 
                                  }
setOfVariables['Total Precipitable Water, Low'] = {
                                  constants.VARIABLE_TECH_NAME_KEY: 'imager_prof_retr_abi_total_precipitable_water_low',
                                  constants.EPSILON_KEY: 3.0,
                                  
                                  # the following settings allow you to control plotting on an individual plot basis
                                  # if these valuse are not defined or set to True, the associated plots will be plotted
                                  # or not according to the other settings for the variable
                                  # if any of these are set to false, the associated plot will not be plotted
                                  
                                  # should the two original images be plotted?
                                  constants.DO_PLOT_ORIGINALS_KEY: True,
                                  # should the absolute value of the difference between the two data sets be plotted?
                                  constants.DO_PLOT_ABS_DIFF_KEY:  True,
                                  # should the (subtractive) difference between the two data sets be plotted?
                                  constants.DO_PLOT_SUB_DIFF_KEY:  True,
                                  # should the scatter plot be plotted?
                                  constants.DO_PLOT_SCATTER_KEY:   True,
                                  # should the density scatter plot be plotted?
                                  constants.DO_PLOT_HEX_KEY:       True,
                                  # should the histogram be plotted?
                                  constants.DO_PLOT_HISTOGRAM_KEY: True,
                                  # should the mismatch plot be plotted?
                                  constants.DO_PLOT_MISMATCH_KEY:  True
                                  }
setOfVariables['Total Precipitable Water, Mid'] = {
                                  constants.VARIABLE_TECH_NAME_KEY: 'imager_prof_retr_abi_total_precipitable_water_mid',
                                  constants.EPSILON_KEY: 3.0
                                  }
setOfVariables['Land Surface Temperature'] = {
                                  constants.VARIABLE_TECH_NAME_KEY: 'imager_prof_retr_abi_land_surface_temperature',
                                  constants.EPSILON_KEY: 5.0
                                  }
setOfVariables['K-Index'] = {
                                  constants.VARIABLE_TECH_NAME_KEY: 'imager_prof_retr_abi_k_index',
                                  constants.EPSILON_KEY: 2.0
                                  }
setOfVariables['Temperature Profile'] = {
                                  constants.VARIABLE_TECH_NAME_KEY: 'imager_prof_retr_abi_tprof',
                                  constants.EPSILON_KEY: 0.1
                                  }
setOfVariables['CAPE'] = {
                                  constants.VARIABLE_TECH_NAME_KEY: 'imager_prof_retr_abi_cape',
                                  constants.EPSILON_KEY: 1000
                                  }
setOfVariables['Lifted Index'] = {
                                  constants.VARIABLE_TECH_NAME_KEY: 'imager_prof_retr_abi_lifted_index',
                                  constants.EPSILON_KEY: 2.0
                                  }
setOfVariables['Showalter Index'] = {
                                  constants.VARIABLE_TECH_NAME_KEY: 'imager_prof_retr_abi_showalter_index',
                                  constants.EPSILON_KEY: 2.0
                                  }
