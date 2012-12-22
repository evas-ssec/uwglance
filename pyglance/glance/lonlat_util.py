#!/usr/bin/env python
# encoding: utf-8
"""
This module contains utility functions specifically related to
longitude and latitude.

Created by evas Dec 2012.
Copyright (c) 2012 University of Wisconsin SSEC. All rights reserved.
"""

import numpy

import glance.data   as dataobj
import glance.plot   as plot
from   glance.util      import get_percentage_from_mask
from   glance.constants import *

# TODO, this comparison needs to encorporate epsilon percent as well
def check_lon_lat_equality(longitudeADataObject, latitudeADataObject,
                           longitudeBDataObject, latitudeBDataObject,
                           llepsilon, doMakeImages, outputPath,
                           fullDPI=None, thumbDPI=None) :
    """
    check to make sure the longitude and latitude are equal everywhere that's not in the ignore masks
    if they are not and doMakeImages was passed as True, generate appropriate figures to show where
    return the number of points where they are not equal (0 would mean they're the same)
    
    If the latitude or longitude cannot be compared, this may raise a VariableComparisonError.
    """
    # first of all, if the latitude and longitude are not the same shape, then things can't ever be "equal"
    if (longitudeADataObject.data.shape != longitudeBDataObject.data.shape) :
        raise VariableComparisonError ("Unable to compare longitue variables due to different sizes (" + str(longitudeADataObject.data.shape) +
                                       ") and (" + str(longitudeBDataObject.data.shape) +").")
    if (latitudeADataObject.data.shape  !=  latitudeBDataObject.data.shape) :
        raise VariableComparisonError ("Unable to compare latitude variables due to different sizes (" + str(latitudeADataObject.data.shape) +
                                       ") and (" + str(latitudeBDataObject.data.shape) +").")
    
    # get information about how the latitude and longitude differ
    longitudeDiffInfo = dataobj.DiffInfoObject(longitudeADataObject, longitudeBDataObject, epsilonValue=llepsilon)
    latitudeDiffInfo  = dataobj.DiffInfoObject(latitudeADataObject,  latitudeBDataObject,  epsilonValue=llepsilon)
    
    # how much difference is there between the two sets?
    lon_lat_not_equal_mask           = longitudeDiffInfo.diff_data_object.masks.mismatch_mask | latitudeDiffInfo.diff_data_object.masks.mismatch_mask
    lon_lat_not_equal_points_count   = numpy.sum(lon_lat_not_equal_mask)
    lon_lat_not_equal_points_percent = (float(lon_lat_not_equal_points_count) / float(lon_lat_not_equal_mask.size)) * 100.0
    
    # if we have unequal points, create user legible info about the problem
    if (lon_lat_not_equal_points_count > 0) :
        LOG.warn("Possible mismatch in values stored in file a and file b longitude and latitude values."
                 + " Depending on the degree of mismatch, some data value comparisons may be "
                 + "distorted or spacially nonsensical.")
        # if we are making images, make two showing the invalid lons/lats
        if (doMakeImages) :
            
            if ((len(longitudeADataObject.data[~longitudeADataObject.masks.ignore_mask]) > 0) and
                (len( latitudeADataObject.data[~ latitudeADataObject.masks.ignore_mask]) > 0)) :
                plot.plot_and_save_spacial_mismatch(longitudeADataObject, latitudeADataObject,
                                                   lon_lat_not_equal_mask,
                                                   "A", "Lon./Lat. Points Mismatched between A and B\n" +
                                                   "(Shown in A)",
                                                   "LonLatMismatch",
                                                   outputPath, True,
                                                   fullDPI=fullDPI, thumbDPI=thumbDPI, units="degrees")
            
            if ((len(longitudeBDataObject.data[~longitudeBDataObject.masks.ignore_mask]) > 0) and
                (len( latitudeBDataObject.data[~ latitudeBDataObject.masks.ignore_mask]) > 0)) :
                plot.plot_and_save_spacial_mismatch(longitudeBDataObject, latitudeBDataObject,
                                                   lon_lat_not_equal_mask,
                                                   "B", "Lon./Lat. Points Mismatched between A and B\n" +
                                                   "(Shown in B)",
                                                   "LonLatMismatch",
                                                   outputPath, True,
                                                   fullDPI=fullDPI, thumbDPI=thumbDPI, units="degrees")
    
    # setup our return data
    returnInfo = {}
    returnInfo[LONLAT_NOT_EQUAL_COUNT_KEY] = lon_lat_not_equal_points_count
    returnInfo[LONLAT_NOT_EQ_PERCENT_KEY]  = lon_lat_not_equal_points_percent
    
    return returnInfo

def compare_spatial_invalidity(longitude_a_object, longitude_b_object,
                               latitude_a_object,  latitude_b_object,
                               spatial_info, do_include_images, output_path,
                               fullDPI=None, thumbDPI=None) :
    """ 
    Given information about where the two files are spatially invalid, figure
    out what invalidity they share and save information or plots for later use
    also build a shared longitude/latitude based on A but also including valid
    points in B
    """
    # make our common invalid masks
    invalid_in_a_mask = longitude_a_object.masks.ignore_mask | latitude_a_object.masks.ignore_mask
    invalid_in_b_mask = longitude_b_object.masks.ignore_mask | latitude_b_object.masks.ignore_mask
    invalid_in_common_mask = invalid_in_a_mask | invalid_in_b_mask
    
    # make a "common" longitude/latitude based on A
    longitude_common = longitude_a_object.data.copy()
    latitude_common  =  latitude_a_object.data.copy()
    
    # compare our spacialy invalid info
    spatial_info[PERCENT_INV_PTS_SHARED_KEY] = spatial_info[A_FILE_TITLE_KEY][PERCENT_INVALID_PTS_KEY]
            # set a default that will hold if the two files have the same spatially invalid pts
    if not numpy.all(invalid_in_a_mask.ravel() == invalid_in_b_mask.ravel()) : 
        LOG.info("Mismatch in number of spatially invalid points. " +
                 "Files may not have corresponding data where expected.")
        
        # figure out which points are only valid in one of the two files
        valid_only_in_mask_a = (~invalid_in_a_mask) & invalid_in_b_mask
        spatial_info[A_FILE_TITLE_KEY][NUMBER_INVALID_PTS_KEY] = numpy.sum(valid_only_in_mask_a.ravel())
        valid_only_in_mask_b = (~invalid_in_b_mask) & invalid_in_a_mask
        spatial_info[B_FILE_TITLE_KEY][NUMBER_INVALID_PTS_KEY] = numpy.sum(valid_only_in_mask_b.ravel())
        
        # so how many do they have together?
        spatial_info[PERCENT_INV_PTS_SHARED_KEY] = get_percentage_from_mask(invalid_in_common_mask)[0]
        # make a "clean" version of the lon/lat
        longitude_common[valid_only_in_mask_a] = longitude_a_object.data[valid_only_in_mask_a]
        longitude_common[valid_only_in_mask_b] = longitude_b_object.data[valid_only_in_mask_b]
        latitude_common [valid_only_in_mask_a] =  latitude_a_object.data[valid_only_in_mask_a]
        latitude_common [valid_only_in_mask_b] =  latitude_b_object.data[valid_only_in_mask_b]
        
        # plot the points that are only valid one file and not the other
        if ((spatial_info[A_FILE_TITLE_KEY][NUMBER_INVALID_PTS_KEY] > 0) and (do_include_images) and
            (len(longitude_a_object.data[~invalid_in_a_mask]) > 0) and
            (len( latitude_a_object.data[~invalid_in_a_mask]) > 0)) :
            plot.plot_and_save_spacial_mismatch(longitude_a_object, latitude_a_object,
                                               valid_only_in_mask_a,
                                               "A", "Points only valid in\nFile A\'s longitude & latitude",
                                               "SpatialMismatch",
                                               output_path, True,
                                               fullDPI=fullDPI, thumbDPI=thumbDPI, units="degrees")
        if ((spatial_info[B_FILE_TITLE_KEY][NUMBER_INVALID_PTS_KEY] > 0) and (do_include_images) and
            (len(longitude_b_object.data[~invalid_in_b_mask]) > 0) and
            (len( latitude_b_object.data[~invalid_in_b_mask]) > 0)
            ) :
            plot.plot_and_save_spacial_mismatch(longitude_b_object, latitude_b_object,
                                               valid_only_in_mask_b,
                                               "B", "Points only valid in\nFile B\'s longitude & latitude",
                                               "SpatialMismatch",
                                               output_path, True,
                                               fullDPI=fullDPI, thumbDPI=thumbDPI, units="degrees")
    
    return invalid_in_common_mask, spatial_info, longitude_common, latitude_common

class VariableComparisonError(Exception):
    """
    The exception raised when a variable could not be compared.
    
        msg  -- explanation of which variable could be compared (and, if possible, why)
    """
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

if __name__=='__main__':
    pass
