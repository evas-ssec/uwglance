#!/usr/bin/env python
# encoding: utf-8
"""
This module handles loading data for some of the top level comparison
routines.

Created by evas Dec 2012.
Copyright (c) 2012 University of Wisconsin SSEC. All rights reserved.
"""

import logging
#from pycdf import CDFError
import numpy

import glance.delta  as delta
import glance.data   as dataobj
import glance.io     as io
from glance.util        import get_percentage_from_mask
from glance.lonlat_util import check_lon_lat_equality, compare_spatial_invalidity
from glance.constants   import *

LOG = logging.getLogger(__name__)

def _get_and_analyze_lon_lat (fileObject,
                              latitudeVariableName, longitudeVariableName,
                              latitudeDataFilterFn=None, longitudeDataFilterFn=None,
                              alternateFilePath=None, fileDescriptior="") :
    """
    get the longitude and latitude data from the given file, assuming they are in the given variable names
    and analyze them to identify spacially invalid data (ie. data that would fall off the earth)
    
    This may result in a ValueError if the variable cannot be loaded.
    """
    
    # for the a file, do we have an alternate?
    fileToUse = fileObject
    if (alternateFilePath is not None) :
        LOG.info("Loading alternate file (" + alternateFilePath
                 + ") for file " + fileDescriptior + " longitude/latitude.")
        fileToUse = dataobj.FileInfo(alternateFilePath)
    
    # get the longitude
    LOG.info ('longitude name: ' + longitudeVariableName)
    lonObject = load_data_object (fileToUse, longitudeVariableName,
                                  rangeMin=-180,
                                  rangeMax=360,
                                  forceDType=numpy.float,
                                  dataFilter=longitudeDataFilterFn)
    
    # get the latitude
    LOG.info ('latitude name: '  + latitudeVariableName)
    latObject = load_data_object (fileToUse, latitudeVariableName,
                                  rangeMin=-90,
                                  rangeMax=90,
                                  forceDType=numpy.float,
                                  dataFilter=latitudeDataFilterFn)
    
    # we are going to have issues with our comparision if they aren't the same shape
    LOG.debug('latitude  shape: ' + str(latObject.data.shape))
    LOG.debug('longitude shape: ' + str(lonObject.data.shape))
    assert (latObject.data.shape == lonObject.data.shape)
    
    # build a mask of our spacially invalid data
    spaciallyInvalidMask = latObject.masks.ignore_mask | lonObject.masks.ignore_mask
    
    # analyze our spacially invalid data
    percentageOfSpaciallyInvalidPts, numberOfSpaciallyInvalidPts = get_percentage_from_mask(spaciallyInvalidMask)
    
    spatialStatInfo = {
                       TOTAL_NUM_INVALID_PTS_KEY: numberOfSpaciallyInvalidPts,
                       PERCENT_INVALID_PTS_KEY:    percentageOfSpaciallyInvalidPts
                       }
    
    return lonObject, latObject, spatialStatInfo

def handle_lon_lat_info (lon_lat_settings, a_file_object, b_file_object, output_path,
                         should_make_images=False, should_check_equality=True,
                         fullDPI=None, thumbDPI=None) :
    """
    Manage loading and comparing longitude and latitude information for two files
    
    This may result in a ValueError if the longitude or latitude cannot be loaded.
    This may result in a VariableComparisonError if the longitude or latitude cannot be compared due to size.
    
    """
    # a place to save some general stats about our lon/lat data
    spatialInfo = { }
    
    # if there is no lon/lat specified, stop now
    if ( (LONGITUDE_NAME_KEY not in lon_lat_settings) or (LATITUDE_NAME_KEY not in lon_lat_settings)
        or ((USE_NO_LON_OR_LAT_VARS_KEY in lon_lat_settings) and lon_lat_settings[USE_NO_LON_OR_LAT_VARS_KEY]) ) :
        return { }, { }
    
    # if we should not be comparing against the logitude and latitude, stop now
    LOG.debug ('lon_lat_settings: ' + str(lon_lat_settings))
    
    # figure out the names to be used for the longitude and latitude variables
    a_longitude_name = lon_lat_settings[LONGITUDE_NAME_KEY]
    a_latitude_name  = lon_lat_settings[LATITUDE_NAME_KEY]
    b_longitude_name = lon_lat_settings[LON_ALT_NAME_IN_B_KEY] if LON_ALT_NAME_IN_B_KEY in lon_lat_settings else a_longitude_name
    b_latitude_name  = lon_lat_settings[LAT_ALT_NAME_IN_B_KEY] if LAT_ALT_NAME_IN_B_KEY in lon_lat_settings else a_latitude_name
    
    # if we need to load our lon/lat from different files, open those files
    longitude_a_object, latitude_a_object, spatialInfo[A_FILE_TITLE_KEY] = \
                          _get_and_analyze_lon_lat (a_file_object,
                                                    a_latitude_name, a_longitude_name,
                              latitudeDataFilterFn=lon_lat_settings[LAT_FILTER_FUNCTION_A_KEY],
                              longitudeDataFilterFn=lon_lat_settings[LON_FILTER_FUNCTION_A_KEY],
                                                    alternateFilePath=lon_lat_settings[LONLAT_ALT_FILE_A_KEY] if (LONLAT_ALT_FILE_A_KEY in lon_lat_settings) else None,
                                                    fileDescriptior="a")
    longitude_b_object, latitude_b_object, spatialInfo[B_FILE_TITLE_KEY] = \
                          _get_and_analyze_lon_lat (b_file_object,
                                                    b_latitude_name, b_longitude_name,
                                                    latitudeDataFilterFn=lon_lat_settings[LAT_FILTER_FUNCTION_B_KEY],
                                                    longitudeDataFilterFn=lon_lat_settings[LON_FILTER_FUNCTION_B_KEY],
                                                    alternateFilePath=lon_lat_settings[LONLAT_ALT_FILE_B_KEY] if (LONLAT_ALT_FILE_B_KEY in lon_lat_settings) else None,
                                                    fileDescriptior="b")
    
    # if we need to, test the level of equality of the "valid" values in our lon/lat
    if should_check_equality :
        
        moreSpatialInfo = check_lon_lat_equality(longitude_a_object, latitude_a_object,
                                                 longitude_b_object, latitude_b_object,
                                                 lon_lat_settings[LON_LAT_EPSILON_KEY],
                                                 should_make_images, output_path,
                                                 fullDPI=fullDPI, thumbDPI=thumbDPI)
        # update our existing spatial information
        spatialInfo.update(moreSpatialInfo)
        
        # compare our spatially invalid info to see if the two files have invalid longitudes and latitudes in the same places
        spaciallyInvalidMask, spatialInfo, longitude_common, latitude_common = \
                                compare_spatial_invalidity(longitude_a_object, longitude_b_object,
                                                           latitude_a_object,  latitude_b_object,
                                                           spatialInfo, should_make_images, output_path,
                                                           fullDPI=fullDPI, thumbDPI=thumbDPI)
    else:
        spaciallyInvalidMask = None
        longitude_common     = None
        latitude_common      = None
    
    # FUTURE, return the lon/lat objects instead?
    return {
            A_FILE_KEY:
                      {
                       LON_KEY:             longitude_a_object.data,
                       LAT_KEY:             latitude_a_object.data,
                       INVALID_MASK_KEY:    longitude_a_object.masks.ignore_mask,
                       LON_FILL_VALUE_KEY:  longitude_a_object.fill_value,
                       LAT_FILL_VALUE_KEY:  latitude_a_object.fill_value
                      },
            B_FILE_KEY:
                      {
                       LON_KEY:             longitude_b_object.data,
                       LAT_KEY:             latitude_b_object.data,
                       INVALID_MASK_KEY:    longitude_b_object.masks.ignore_mask,
                       LON_FILL_VALUE_KEY:  longitude_b_object.fill_value,
                       LAT_FILL_VALUE_KEY:  latitude_b_object.fill_value
                      },
            COMMON_KEY:
                      {
                       LON_KEY:             longitude_common,
                       LAT_KEY:             latitude_common,
                       INVALID_MASK_KEY:    spaciallyInvalidMask
                      }
            }, \
           spatialInfo

def handle_lon_lat_info_for_one_file (lon_lat_settings, file_object) :
    """
    Manage loading longitude and latitude information for a file
    
    This may result in a ValueError if the longitude or latitude cannot be loaded.
    
    """
    
    # if there is no lon/lat specified, stop now
    if ( (LONGITUDE_NAME_KEY not in lon_lat_settings) or (LATITUDE_NAME_KEY not in lon_lat_settings)
        or ((USE_NO_LON_OR_LAT_VARS_KEY in lon_lat_settings) and lon_lat_settings[USE_NO_LON_OR_LAT_VARS_KEY]) ) :
        return { }, { }
    
    # print our settings for debugging purposes
    LOG.debug ('lon_lat_settings: ' + str(lon_lat_settings))
    
    # figure out the names to be used for the longitude and latitude variables
    lon_name = lon_lat_settings[LONGITUDE_NAME_KEY]
    lat_name = lon_lat_settings[LATITUDE_NAME_KEY ]
    
    # load our lon/lat data
    
    lon_object, lat_object, spatialInfo = \
                        _get_and_analyze_lon_lat (file_object,
                                                  lat_name, lon_name,
                                                  latitudeDataFilterFn=lon_lat_settings[LAT_FILTER_FUNCTION_A_KEY],
                                                  longitudeDataFilterFn=lon_lat_settings[LON_FILTER_FUNCTION_A_KEY],
                                                  alternateFilePath=lon_lat_settings[LONLAT_ALT_FILE_A_KEY] if (LONLAT_ALT_FILE_A_KEY in lon_lat_settings) else None)
    
    # FUTURE, return the lon/lat objects instead?
    return {
            LON_KEY:             lon_object.data,
            LAT_KEY:             lat_object.data,
            INVALID_MASK_KEY:    lon_object.masks.ignore_mask | lat_object.masks.ignore_mask,
            LON_FILL_VALUE_KEY:  lon_object.fill_value,
            LAT_FILL_VALUE_KEY:  lat_object.fill_value
            }, \
           spatialInfo

def open_and_process_files (fileNames):
    """
    open files listed in the args and get information about the variables in them
    """
    # open all the files & get their variable names
    files = {}
    commonNames = None
    for fileName in fileNames:
        LOG.info("opening %s" % fileName)
        files[fileName] = {}
        tempFileObject = (io.open(fileName))
        files[fileName][FILE_OBJECT_KEY] = tempFileObject
        tempNames = set(tempFileObject())
        LOG.debug ('variable names for ' + fileName + ': ' + str(tempNames)) 
        files[fileName][FILE_VARIABLE_NAMES_KEY] = tempNames
        commonNames = tempNames if commonNames is None else commonNames.intersection(tempNames)
    files[COMMON_VAR_NAMES_KEY] = commonNames
    
    return files

def load_variable_data(fileObject, variableNameInFile,
                       forceDType=None,
                       dataFilter=None,
                       variableToFilterOn=None,
                       variableBasedFilter=None,
                       altVariableFileObject=None,
                       fileDescriptionForDisplay="file",
                       correctForAWIPS=False) :
    """
    load data for a variable from a file
    optionally filter the variable data based on a data filter or another variable
    
    dataFilter must be in the form of (lambda data: some manipulation returning the new data)
    variableBasedFilter must be in the form of (lambda data, filterData: some manipulation returning the new data))
    """
    
    variableData     = None
    exceptionToRaise = None
    
    # get the data for the variable
    LOG.debug("loading basic data for variable " + variableNameInFile + " from " + fileDescriptionForDisplay)
    if fileObject is None :
        exceptionToRaise = ValueError("File was not properly opened so variable '" + variableNameInFile + "' could not be loaded.")
    else :
        try :
            variableData = numpy.array(fileObject[variableNameInFile]) if forceDType is None else numpy.array(fileObject[variableNameInFile], dtype=forceDType)
            variableData = variableData.astype(numpy.uint8) if correctForAWIPS else variableData
        except Exception, ex :
            import traceback
            exceptionToRaise = ValueError('Unable to retrieve ' + variableNameInFile + ' data. The variable name ' + 
                      ' may not exist in this file or an error may have occured while attempting to' +
                      ' access the data. Details of file access error observed: ' + str(ex))
    
    # if we ended up with an exception, raise that now
    if exceptionToRaise is not None :
        raise exceptionToRaise
    
    # apply the basic filter if there is one
    if dataFilter is not None :
        LOG.debug ("applying filter function to data from " + fileDescriptionForDisplay + " for variable " + variableNameInFile)
        variableData = dataFilter(variableData)
    
    # if we've got another variable to filter on, do that
    if (variableToFilterOn is not None) and (variableBasedFilter is not None) :
        LOG.debug ("filtering data from " + fileDescriptionForDisplay + " for variable " + variableNameInFile
                   + " based on additional data from variable " + variableToFilterOn)
        
        fileToUseTemp = fileObject
        if altVariableFileObject is not None :
            fileToUseTemp = altVariableFileObject # TODO, is this the right kind of object?
        
        dataToFilterOn = fileToUseTemp[variableToFilterOn]
        variableData   = variableBasedFilter(variableData, dataToFilterOn)
    
    return variableData

def load_data_object (fileObject, variableNameInFile,
                      rangeMin=None,
                      rangeMax=None,
                      forceDType=None,
                      dataFilter=None,
                      variableToFilterOn=None,
                      variableBasedFilter=None,
                      altVariableFileObject=None,
                      fileDescriptionForDisplay="file") :
    """
    load the data and put it into an appropriate DataObject
    
    Note: the rangeMin and rangeMax are the minimum and maximum acceptable data
    values used to calculate the invalid data mask
    """
    # get the data from the file
    rawData = load_variable_data(fileObject.file_object,
                                 variableNameInFile,
                                 forceDType=forceDType,
                                 dataFilter=dataFilter,
                                 variableToFilterOn=variableToFilterOn,
                                 variableBasedFilter=variableBasedFilter,
                                 altVariableFileObject=altVariableFileObject,
                                 fileDescriptionForDisplay=fileDescriptionForDisplay)
    
    # get the fill value
    fillValue = fileObject.file_object.missing_value(variableNameInFile)
    
    # build a mask of invalid data
    invalidMask = ~numpy.isfinite(rawData)
    if rangeMin  is not None :
        invalidMask |= rawData < rangeMin
    if rangeMax  is not None :
        invalidMask |= rawData > rangeMax
    if fillValue is not None :
        invalidMask |= rawData == fillValue
    
    return dataobj.DataObject(rawData, fillValue=fillValue, ignoreMask=invalidMask)

def get_UV_info_from_magnitude_direction_info(fileObject, magnitudeName, directionName, invalidMask=None) :
    """
    If there are magnitude and direction names, load that information and calculate the u and v that correspond to it
    """
    
    # if we don't have magnitude and direction, we can't calculate the U and V values
    if (magnitudeName is None) or (directionName is None) :
        return None, None
    
    # load the magnitude and direction data sets
    magnitude = load_variable_data(fileObject, magnitudeName)
    direction = load_variable_data(fileObject, directionName)
    
    # convert the magnitude and direction data into u and v vectors
    uData, vData = delta.convert_mag_dir_to_U_V_vector(magnitude, direction, invalidMask=invalidMask)
    
    return uData, vData

if __name__=='__main__':
    pass
