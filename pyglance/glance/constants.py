#!/usr/bin/env python
# encoding: utf-8
"""
Constant strings for use in glance. Most of these are keys for various dictionaries.

Created by evas Dec 2012.
Copyright (c) 2012 University of Wisconsin SSEC. All rights reserved.
"""

# most of these are keys that allow you to map to values
# in glance's informational dictionaries; the comments
# generally describe the value that the key will lead to

# constants related to controlling large parts of the run

DO_MAKE_IMAGES_KEY         = 'shouldIncludeImages'
DO_MAKE_REPORT_KEY         = 'shouldIncludeReport'
DO_COLOCATION_KEY          = 'doColocate'
DO_MAKE_FORKS_KEY          = 'doFork'
DO_CLEAR_MEM_THREADED_KEY  = 'useThreadsToControlMemory'
USE_SHARED_ORIG_RANGE_KEY  = 'useSharedRangeForOriginal'
DO_TEST_PASSFAIL_KEY       = 'usePassFail'
DO_IMAGES_ONLY_ON_FAIL_KEY = 'only_plot_on_fail'
USE_NO_LON_OR_LAT_VARS_KEY = 'noLonLatVars'
SHORT_CIRCUIT_DIFFS_KEY    = 'short_circuit_diffs'
USE_CUSTOM_PROJ_KEY        = 'use_custom_projection'

# constants related to storing information from the run

MACHINE_INFO_KEY           = 'machine'
USER_INFO_KEY              = 'user'
GLANCE_VERSION_INFO_KEY    = 'version'
TIME_INFO_KEY              = 'time'

DID_VARIABLE_PASS_KEY      = 'did_pass'

# the base directory where the variable is
VARIABLE_DIRECTORY_KEY     = 'variable_dir'
# the path to the variable report
VAR_REPORT_PATH_KEY        = 'variable_report_path_escaped'
# the path to the documentation report
DOCUMENTATION_PATH_KEY     = 'doc_path'
# the path where a copy of the config file will be
CONFIG_FILE_PATH_KEY       = 'config_file_path'
# the name of the configuration file
CONFIG_FILE_NAME_KEY       = 'config_file_name'

# high level pass fail settings

EPSILON_FAIL_TOLERANCE_KEY = 'epsilon_failure_tolerance'
NONFINITE_TOLERANCE_KEY    = 'nonfinite_data_tolerance'
TOTAL_FAIL_TOLERANCE_KEY   = 'total_data_failure_tolerance'
MIN_OK_R_SQUARED_COEFF_KEY = 'minimum_acceptable_squared_correlation_coefficient'

# constants related to lon/lat information during a run

LONGITUDE_NAME_KEY         = 'longitude'
LATITUDE_NAME_KEY          = 'latitude'
LON_ALT_NAME_IN_B_KEY      = 'longitude_alt_name_in_b'
LAT_ALT_NAME_IN_B_KEY      = 'latitude_alt_name_in_b'

LON_LAT_EPSILON_KEY        = 'lon_lat_epsilon'

LON_FILTER_FUNCTION_A_KEY  = 'data_filter_function_lon_in_a'
LAT_FILTER_FUNCTION_A_KEY  = 'data_filter_function_lat_in_a'
LON_FILTER_FUNCTION_B_KEY  = 'data_filter_function_lon_in_b'
LAT_FILTER_FUNCTION_B_KEY  = 'data_filter_function_lat_in_b'

LONLAT_ALT_FILE_A_KEY      = 'a_lon_lat_from_alt_file'
LONLAT_ALT_FILE_B_KEY      = 'b_lon_lat_from_alt_file'

# constants related to variable information during a run

# the name of the variable as it should appear in the data file
VARIABLE_TECH_NAME_KEY     = 'variable_name'
# if there is an alternate version of the variable name in the b file, use this key
VARIABLE_B_TECH_NAME_KEY   = 'alternate_name_in_B'
# internally used to hold the name that will be shown to the user
DISPLAY_NAME_KEY           = 'display_name'

# if the data is a vector, the magnitude and direction are needed
MAGNITUDE_VAR_NAME_KEY     = 'magnitudeName'
DIRECTION_VAR_NAME_KEY     = 'directionName'
MAGNITUDE_B_VAR_NAME_KEY   = 'magnitudeBName'
DIRECTION_B_VAR_NAME_KEY   = 'directionBName'

# if we're supposed to do bin/tuple analysis, info on the bin and tuple
BIN_INDEX_KEY              = 'binIndex'
TUPLE_INDEX_KEY            = 'tupleIndex'
BIN_NAME_KEY               = 'binName'
TUPLE_NAME_KEY             = 'tupleName'

# the variable's fill value
FILL_VALUE_KEY             = 'missing_value'
# if there's an alternate fill value in b, it'll be here
FILL_VALUE_ALT_IN_B_KEY    = 'missing_value_alt_in_b'

# the units for the variable (usually gotten from the file)
VAR_UNITS_A_KEY            = 'units_a'
VAR_UNITS_B_KEY            = 'units_b'

# a list of numerical ranges that the data should be displayed in (discrete chunks)
DISPLAY_RANGES_KEY         = 'display_ranges'
# a list of names to be displayed with the display ranges
DISPLAY_RANGE_NAMES_KEY    = 'display_range_names'
# a list of colors to be used for the display ranges
DISPLAY_RANGE_COLORS_KEY   = 'display_colors'

# the range the histogram should be displayed in (data outside the range will not be shown)
HISTOGRAM_RANGE_KEY        = 'histogram_range'

# the epsilon for comparing the variable
EPSILON_KEY                = 'epsilon'
# another way to define epsilon: the % of A's value that A and B can be different
EPSILON_PERCENT_KEY        = 'epsilon_percent'

# filter functions to use to filter the data
FILTER_FUNCTION_A_KEY      = 'data_filter_function_a'
FILTER_FUNCTION_B_KEY      = 'data_filter_function_b'
# filter functions to filter the data based on a different variable
VAR_FILTER_FUNCTION_A_KEY  = 'variable_based_filter_a'
VAR_FILTER_NAME_A_KEY      = 'variable_to_filter_on_a'
VAR_FILTER_ALT_FILE_A_KEY  = 'variable_to_filter_alt_file_a'
VAR_FILTER_FUNCTION_B_KEY  = 'variable_based_filter_b'
VAR_FILTER_NAME_B_KEY      = 'variable_to_filter_on_b'
VAR_FILTER_ALT_FILE_B_KEY  = 'variable_to_filter_alt_file_b'

# constants related to plotting settings

DO_PLOT_HISTOGRAM_KEY      = 'do_plot_histogram'
DO_PLOT_SCATTER_KEY        = 'do_plot_scatter'
DO_PLOT_HEX_KEY            = 'do_plot_hex'
DO_PLOT_ORIGINALS_KEY      = 'do_plot_originals'
DO_PLOT_ABS_DIFF_KEY       = 'do_plot_abs_diff'
DO_PLOT_SUB_DIFF_KEY       = 'do_plot_sub_diff'
DO_PLOT_MISMATCH_KEY       = 'do_plot_mismatch'

DETAIL_DPI_KEY             = 'detail_DPI'
THUMBNAIL_DPI_KEY          = 'thumb_DPI'

# constants for keying the functions in the plot function list

HIST_FUNCTION_KEY          = 'histogram'
SCATTER_FUNCTION_KEY       = 'scatter'
DENSITY_SCATTER_FN_KEY     = 'density-scatter'
MULTI_SCATTER_FUNCTION_KEY = 'multi-scatter'
HEX_PLOT_FUNCTION_KEY      = 'scatterD'
ORIG_FUNCTION_KEY          = 'original'
ORIG_A_FUNCTION_KEY        = 'originalA'
ORIG_B_FUNCTION_KEY        = 'originalB'
ABS_DIFF_FUNCTION_KEY      = 'diffAbs'
SUB_DIFF_FUNCTION_KEY      = 'diffSub'
MISMATCH_FUNCTION_KEY      = 'mismatch'

# constants for the lon/lat structure and the paths structure

A_FILE_KEY                 = 'a'
B_FILE_KEY                 = 'b'
COMMON_KEY                 = 'common'
LON_KEY                    = 'lon'
LAT_KEY                    = 'lat'
INVALID_MASK_KEY           = 'inv_mask'
OUT_FILE_KEY               = 'out'
LON_FILL_VALUE_KEY         = 'lon_fill'
LAT_FILL_VALUE_KEY         = 'lat_fill'

# constants for the files structure

A_FILE_TITLE_KEY           = 'file A'
B_FILE_TITLE_KEY           = 'file B'
PATH_KEY                   = 'path'
LAST_MODIFIED_KEY          = 'lastModifiedTime'
MD5SUM_KEY                 = 'md5sum'

# TEMP these are for the stats file structure that will eventually be removed
FILE_OBJECT_KEY            = 'fileObject'
FILE_VARIABLE_NAMES_KEY    = 'varNames'
COMMON_VAR_NAMES_KEY       = 'commonVarNames'

# constants related to name analysis

# this is only used when a single file is analyzed
POSSIBLE_NAMES_KEY         = 'possibleNames'
# these are more general name analysis constants
SHARED_VARIABLE_NAMES_KEY  = 'sharedVars'
VAR_NAMES_UNIQUE_TO_A_KEY  = 'uniqueToAVars'
VAR_NAMES_UNIQUE_TO_B_KEY  = 'uniqueToBVars'

# high level constant sets saved for the report generation

VARIABLE_RUN_INFO_KEY      = 'variable_run_info'
# related to pass/fail testing and display
PASSED_EPSILON_PERCENT_KEY = 'pass_epsilon_percent'
FINITE_SIMILAR_PERCENT_KEY = 'finite_similar_percent'
R_SQUARED_COEFF_VALUE_KEY  = 'r_squared_correlation'

# image types

ORIGINAL_IMAGES_KEY        = 'original'
COMPARED_IMAGES_KEY        = 'compared' # TODO, got to here in searching to make sure I replaced all uses

# keys used for spatial information

TOTAL_NUM_INVALID_PTS_KEY  = 'totNumInvPts'
PERCENT_INVALID_PTS_KEY    = 'perInvPts'
NUMBER_INVALID_PTS_KEY     = 'numInvPts'
PERCENT_INV_PTS_SHARED_KEY = 'perInvPtsInBoth'

LONLAT_NOT_EQUAL_COUNT_KEY = 'lon_lat_not_equal_points_count'
LONLAT_NOT_EQ_PERCENT_KEY  = 'lon_lat_not_equal_points_percent'

# values only used in the parsing of the command line options structure
# FUTURE, eventually phase out direct use of this structure
# and only have these constants in the config_organizer?

OPTIONS_FILL_VALUE_KEY     = 'missing'
OPTIONS_OUTPUT_PATH_KEY    = 'outputpath'
OPTIONS_CONFIG_FILE_KEY    = 'configFile'
OPTIONS_NO_REPORT_KEY      = 'imagesOnly'
OPTIONS_NO_IMAGES_KEY      = 'htmlOnly'
OPTIONS_LAT_VAR_NAME_KEY   = 'latitudeVar'
OPTIONS_LON_VAR_NAME_KEY   = 'longitudeVar'
OPTIONS_LONLAT_EPSILON_KEY = 'lonlatepsilon'

# values used by the reports

RUN_INFO_DICT_KEY          = 'runInfo'
FILES_INFO_DICT_KEY        = 'files'
STATS_INFO_DICT_KEY        = 'statGroups'
SPATIAL_INFO_DICT_KEY      = 'spatial'
IMAGE_NAME_INFO_DICT_KEY   = 'imageNames'
VARIABLE_NAMES_DICT_KEY    = 'varNames'
VARIABLE_RUN_INFO_DICT_KEY = 'variables'

DEFINITIONS_INFO_KEY       = 'definitions'

if __name__=='__main__':
    pass
