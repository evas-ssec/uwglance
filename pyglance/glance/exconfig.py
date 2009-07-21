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

# whether or not images should be generated and shown in the report
shouldIncludeImages = False

# the names of the latitude and longitude variables that will be used
lat_lon_info = {}
lat_lon_info['longitude'] = 'imager_prof_retr_abi_r4_generic1' # 'pixel_longitude'# the name of the longitude variable
lat_lon_info ['latitude'] = 'imager_prof_retr_abi_r4_generic2' # 'pixel_latitude' # the name of the latitude variable
# the following two values are optional and only need to be set if the the latitude and longitude have
# different names in file A and file B
"""
lat_lon_info['longitude_alt_name_in_b'] = 'resampled_longitude' # the alternate name of the longitude in file B
lat_lon_info ['latitude_alt_name_in_b'] = 'resampled_latitude'  # the alternate name of the latitude in file A
"""
# this value can be used to control how similar the longitude and latitude must be to be considered matching
# if all of your longitude and latitude do not match under this epsilon, most of the comparison report will
# not be generated, since the data would not correlate spatially
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

setOfVariables['imager_prof_retr_abi_total_totals_index'] = {           # this should match the variable name in your files
                                  
                                  'display_name': 'Total Totals',       # this entry is totally optional, it's used to label
                                                                        # the report and some of the plots, if ommitted
                                                                        # the variable name will be used in those places instead
                                                                        
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
                                  }
setOfVariables['imager_prof_retr_abi_total_precipitable_water_high'] = {
                                  'display_name': 'Total Precipitable Water, High',
                                  'epsilon': 3.0,
                                  
                                  'alternate_name_in_B': 'imager_prof_retr_abi_total_precipitable_water_low',
                                                                        # this represents an alternate name that would be
                                                                        # the equivalent variable in file B (if this is
                                                                        # defined, our primary variable name will be expected
                                                                        # to appear only in file A)
                                  }

setOfVariables['imager_prof_retr_abi_total_precipitable_water'] = {
                                  'display_name': 'Total Precipitable Water',
                                  'epsilon': 3.0
                                  # example:
                                  # because missing, and the two tolerances are not defined here,
                                  # this variable would use the defaultValues for those 
                                  }
setOfVariables['imager_prof_retr_abi_total_precipitable_water_low'] = {
                                  'display_name': 'Total Precipitable Water, Low',
                                  'epsilon': 3.0
                                  }
setOfVariables['imager_prof_retr_abi_total_precipitable_water_mid'] = {
                                  'display_name': 'Total Precipitable Water, Mid',
                                  'epsilon': 3.0
                                  }
setOfVariables['imager_prof_retr_abi_land_surface_temperature'] = {
                                  'display_name': 'Land Surface Temperature',
                                  'epsilon': 5.0
                                  }
setOfVariables['imager_prof_retr_abi_k_index'] = {
                                  'display_name': 'K-Index',
                                  'epsilon': 2.0
                                  }
setOfVariables['imager_prof_retr_abi_tprof'] = {
                                  'display_name': 'Temperature Profile',
                                  'epsilon': 0.1
                                  }
setOfVariables['imager_prof_retr_abi_cape'] = {
                                  'display_name': 'CAPE',
                                  'epsilon': 1000
                                  }
setOfVariables['imager_prof_retr_abi_lifted_index'] = {
                                  'display_name': 'Lifted Index',
                                  'epsilon': 2.0
                                  }
setOfVariables['imager_prof_retr_abi_showalter_index'] = {
                                  'display_name': 'Showalter Index',
                                  'epsilon': 2.0
                                  }