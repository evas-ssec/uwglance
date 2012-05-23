#!/usr/bin/env python
# encoding: utf-8
"""

Top-level routines to compare two files.


Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging, re, subprocess, datetime
import imp as imp
from pprint import pprint, pformat
from numpy import *
import pkg_resources
from pycdf import CDFError
from subprocess import check_call as sh
from urllib import quote

import matplotlib
# this is a hack to keep glance from needing pyqt unless you run the gui
if "gui" in sys.argv[1:] :
    try :
        matplotlib.use('Qt4Agg')
        import glance.gui_controller as gui_control
    except ImportError :
        print ("*** Unable to import PyQt4. Please install PyQt4 and add it to your PYTHONPATH in order to use the Glance GUI. ***")
        raise
else :
    matplotlib.use('Agg')

import glance.io     as io
import glance.delta  as delta
import glance.data   as dataobj
import glance.plot   as plot
import glance.report as report
import glance.stats  as statistics
import glance.plotcreatefns as plotcreate
import glance.collocation   as collocation

LOG = logging.getLogger(__name__)

# these are the built in defaults for the settings
glance_setting_defaults = {'shouldIncludeReport':       True,
                           'shouldIncludeImages':       False,
                           'doFork':                    False,
                           'useThreadsToControlMemory': False,
                           'useSharedRangeForOriginal': False,
                           'noLonLatVars':              False,
                           'detail_DPI':                150,
                           'thumb_DPI':                 50}

# these are the built in longitude/latitude defaults
glance_lon_lat_defaults = {'longitude': 'pixel_longitude',
                           'latitude':  'pixel_latitude',
                           'lon_lat_epsilon': 0.0,
                           'data_filter_function_lon_in_a': None,
                           'data_filter_function_lat_in_a': None,
                           'data_filter_function_lon_in_b': None,
                           'data_filter_function_lat_in_b': None
                           }

# these are the built in default settings for the variable analysis
glance_analysis_defaults = {'epsilon': 0.0,
                            'epsilon_percent': None,
                            'missing_value': None,
                            'epsilon_failure_tolerance': 0.0,
                            'nonfinite_data_tolerance':  0.0,
                            'total_data_failure_tolerance': None,
                            'minimum_acceptable_squared_correlation_coefficient': None,
                            'only_plot_on_fail': False
                            }

def _clean_path(string_path) :
    """
    Return a clean form of the path without any '.', '..', or '~'
    """
    clean_path = None
    if string_path is not None :
        clean_path = os.path.abspath(os.path.expanduser(string_path))
    
    return clean_path

def _parse_varnames(names, terms, epsilon=0.0, missing=None):
    """filter variable names and substitute default epsilon and missing settings if none provided
    returns (variable name, epsilon, missing) triples
    
    >>> _parse_varnames( ['foo','bar', 'baz', 'zoom', 'cat'], ['f..:0.5:-999', 'ba.*:0.001', 'c.t::-9999'], 1e-7 )
    set([('foo', 0.5, -999.0), ('cat', 9.9999999999999995e-08, -9999.0), ('bar', 0.001, None), ('baz', 0.001, None)])
    
    names   - all the variable names in the file (ie. names that should be considered valid)
    terms   - variable selection terms given from the command line
    epsilon - a default epsilon to be used for all variables that do not have a specific epsilon given
    missing - a default fill value to be used for all variables that do not have a specific fill value given
    """
    terms = [x.split(':') for x in terms]
    terms = [(re.compile(x[0]).match,x[1:]) for x in terms]
    def _cvt_em(eps=None, mis=None):
        eps = float(eps) if eps else epsilon
        mis = float(mis) if mis else missing
        return eps, mis
    sel = [ ((x,)+_cvt_em(*em)) for x in names for (t,em) in terms if t(x) ]
    return set(sel)

def _check_file_names(fileAObject, fileBObject) :
    """
    get information about the names in the two files and how they compare to each other
    """
    # get information about the variables stored in the files
    aNames = set(fileAObject())
    bNames = set(fileBObject())
    
    # get the variable names they have in common
    commonNames = aNames.intersection(bNames)
    # which names are unique to only one of the two files?
    uniqueToANames = aNames - commonNames
    uniqueToBNames = bNames - commonNames
    
    return _check_shared_names(set(fileAObject()), set(fileBObject()))

def _check_shared_names (nameSetA, nameSetB) :
    """
    compare the names in the two sets
    """
    
    # what names do they have in common?
    commonNames = nameSetA.intersection(nameSetB)
    # what names are unique to each set?
    uniqueToANames = nameSetA - commonNames
    uniqueToBNames = nameSetB - commonNames
    
    return {'sharedVars': commonNames,  'uniqueToAVars': uniqueToANames, 'uniqueToBVars': uniqueToBNames}

def _resolve_names(fileAObject, fileBObject, defaultValues,
                   requestedNames, usingConfigFileFormat=False) :
    """
    figure out which names the two files share and which are unique to each file, as well as which names
    were requested and are in both sets
    
    usingConfigFileFormat signals whether the requestedNames parameter will be in the form of the inputed
    names from the command line or a more complex dictionary holding information about the names read in
    from a configuration file
    
    Note: if we ever need a variable with different names in file A and B to be comparable, this logic
    will need to be changed.
    """
    # look at the names present in the two files and compare them
    nameComparison = _check_file_names(fileAObject, fileBObject)
    
    # figure out which set should be selected based on the user requested names
    fileCommonNames = nameComparison['sharedVars']
    finalNames = {}
    if (usingConfigFileFormat) :
        
        # if the user didn't ask for any, try everything
        if (len(requestedNames) is 0) :
            finalFromCommandLine = _parse_varnames(fileCommonNames, ['.*'],
                                                   defaultValues['epsilon'], defaultValues['missing_value'])
            for name, epsilon, missing in finalFromCommandLine :
                # we'll use the variable's name as the display name for the time being
                finalNames[name] = {}
                # make sure we pick up any other controlling defaults
                finalNames[name].update(defaultValues) 
                # but override the values that would have been determined by _parse_varnames
                finalNames[name]['variable_name'] = name
                finalNames[name]['epsilon'] = epsilon
                
                # load the missing value if it was not provided
                missing, missing_b = _get_missing_values_if_needed((fileAObject, fileBObject), name,
                                                                   missing_value_A=missing, missing_value_B=missing)
                finalNames[name]['missing_value'] = missing 
                finalNames[name]['missing_value_alt_in_b'] = missing_b
                
                # get any information about the units listed in the files
                finalNames[name]['units_a'] = fileAObject.get_attribute(name, io.UNITS_CONSTANT)
                finalNames[name]['units_b'] = fileBObject.get_attribute(name, io.UNITS_CONSTANT)
                
        # otherwise just do the ones the user asked for
        else : 
            # check each of the names the user asked for to see if it is either in the list of common names
            # or, if the user asked for an alternate name mapping in file B, if the two mapped names are in
            # files A and B respectively
            for dispName in requestedNames :
                
                # hang on to info on the current variable
                currNameInfo = requestedNames[dispName] 
                
                # get the variable name 
                if 'variable_name' in currNameInfo :
                    name = currNameInfo['variable_name']
                    name_b = name
                    
                    if ('alternate_name_in_B' in currNameInfo) :
                        name_b = currNameInfo['alternate_name_in_B']
                    
                    if ( (name in fileCommonNames) and (not currNameInfo.has_key('alternate_name_in_B')) ) or \
                            ( (currNameInfo.has_key('alternate_name_in_B') and
                              ((name   in nameComparison['uniqueToAVars']) or (name   in fileCommonNames))  and
                              ((name_b in nameComparison['uniqueToBVars']) or (name_b in fileCommonNames))) ) :
                        finalNames[dispName] = defaultValues.copy() 
                        finalNames[dispName]['display_name'] = dispName
                        finalNames[dispName].update(currNameInfo)
                        
                        # load the missing value if it was not provided
                        missing = finalNames[dispName]['missing_value']
                        if ('missing_value_alt_in_b' in finalNames[dispName]) :
                            missing_b = finalNames[dispName]['missing_value_alt_in_b']
                        else :
                            missing_b = missing
                        finalNames[dispName]['missing_value'], finalNames[dispName]['missing_value_alt_in_b'] = \
                                    _get_missing_values_if_needed((fileAObject, fileBObject), name, name_b,
                                                                  missing, missing_b)
                        
                        # get any information about the units listed in the files
                        finalNames[dispName]['units_a'] = fileAObject.get_attribute(name,   io.UNITS_CONSTANT)
                        finalNames[dispName]['units_b'] = fileBObject.get_attribute(name_b, io.UNITS_CONSTANT)
                        
                else :
                    LOG.warn('No technical variable name was given for the entry described as "' + dispName + '". ' +
                             'Skipping this variable.')
    else:
        # format command line input similarly to the stuff from the config file
        #print (requestedNames)
        finalFromCommandLine = _parse_varnames(fileCommonNames, requestedNames,
                                               defaultValues['epsilon'], defaultValues['missing_value'])
        for name, epsilon, missing in finalFromCommandLine :
            ## we'll use the variable's name as the display name for the time being
            finalNames[name] = {}
            # make sure we pick up any other controlling defaults
            finalNames[name].update(defaultValues) 
            # but override the values that would have been determined by _parse_varnames
            finalNames[name]['variable_name'] = name
            finalNames[name]['epsilon'] = epsilon
            
            # load the missing value if it was not provided
            missing, missing_b = _get_missing_values_if_needed((fileAObject, fileBObject), name,
                                                               missing_value_A=missing, missing_value_B=missing)
            finalNames[name]['missing_value'] = missing 
            finalNames[name]['missing_value_alt_in_b'] = missing_b
            
            # get any information about the units listed in the files
            finalNames[name]['units_a'] = fileAObject.get_attribute(name, io.UNITS_CONSTANT)
            finalNames[name]['units_b'] = fileBObject.get_attribute(name, io.UNITS_CONSTANT)
    
    LOG.debug("Final selected set of variables to analyze:")
    LOG.debug(str(finalNames))
    
    return finalNames, nameComparison

def _resolve_names_one_file(fileObject, defaultValues,
                   requestedNames, usingConfigFileFormat=False) :
    """
    sort out which names to examine based on a file that contains names and the names
    the caller asked for, then fill in information on missing values based on the
    caller requests, possible config file, and defaults
    """
    # look at the names present in the file
    possibleNames = set(fileObject())
    
    # figure out which names should be selected based on the user requested names
    finalNames = {}
    if (usingConfigFileFormat) :
        
        # if the user didn't ask for any, try everything
        if (len(requestedNames) is 0) :
            finalFromCommandLine = _parse_varnames(possibleNames, ['.*'],
                                                   None, defaultValues['missing_value'])
            for name, _, missing in finalFromCommandLine :
                # we'll use the variable's name as the display name for the time being
                finalNames[name] = {}
                # make sure we pick up any other controlling defaults
                finalNames[name].update(defaultValues) 
                # but override the values that would have been determined by _parse_varnames
                finalNames[name]['variable_name'] = name
                
                # load the missing value if it was not provided
                if missing is None :
                    missing = fileObject.missing_value(name)
                finalNames[name]['missing_value'] = missing
                
                # get any information about the units listed in the file
                finalNames[name]['units'] = fileObject.get_attribute(name, io.UNITS_CONSTANT)
                
        # otherwise just do the ones the user asked for
        else : 
            # check each of the names the user asked for to see if it's among the possible names
            for dispName in requestedNames :
                
                # hang on to info on the current variable
                currNameInfo = requestedNames[dispName] 
                
                # get the variable name 
                if 'variable_name' in currNameInfo :
                    name = currNameInfo['variable_name']
                    
                    if (name in possibleNames) :
                        finalNames[dispName] = defaultValues.copy() 
                        finalNames[dispName]['display_name'] = dispName
                        finalNames[dispName].update(currNameInfo)
                        
                        # load the missing value if it was not provided
                        missing = finalNames[dispName]['missing_value']
                        if missing is None :
                            missing = fileObject.missing_value(name)
                        finalNames[dispName]['missing_value'] = missing
                        
                        # get any information about the units listed in the file
                        finalNames[dispName]['units'] = fileObject.get_attribute(name, io.UNITS_CONSTANT)
                        
                else :
                    LOG.warn('No technical variable name was given for the entry described as "' + dispName + '". ' +
                             'Skipping this variable.')
    else:
        # format command line input similarly to the stuff from the config file
        #print (requestedNames)
        finalFromCommandLine = _parse_varnames(possibleNames, requestedNames,
                                               None, defaultValues['missing_value'])
        for name, _, missing in finalFromCommandLine :
            ## we'll use the variable's name as the display name for the time being
            finalNames[name] = {}
            # make sure we pick up any other controlling defaults
            finalNames[name].update(defaultValues) 
            # but override the values that would have been determined by _parse_varnames
            finalNames[name]['variable_name'] = name
            
            # load the missing value if it was not provided
            if missing is None :
                missing = fileObject.missing_value(name)
            finalNames[name]['missing_value'] = missing
            
            # get any information about the units listed in the file
            finalNames[name]['units'] = fileObject.get_attribute(name, io.UNITS_CONSTANT)
    
    LOG.debug("Final selected set of variables to inspect:")
    LOG.debug(str(finalNames))
    
    return finalNames, possibleNames

def _get_missing_values_if_needed((fileA, fileB),
                                  var_name, alt_var_name=None, 
                                  missing_value_A=None, missing_value_B=None) :
    """
    get the missing values for two files based on the variable name(s)
    if the alternate variable name is passed it will be used for the
    second file in place of the primary variable name
    """
    # if we don't have an alternate variable name, use the existing one
    if alt_var_name is None :
        alt_var_name = var_name
    
    if missing_value_A is None :
        missing_value_A = fileA.missing_value(var_name)
    
    if missing_value_B is None :
        missing_value_B = fileB.missing_value(alt_var_name)
    
    return missing_value_A, missing_value_B

def _load_config_or_options(aPath, bPath, optionsSet, requestedVars = [ ]) :
    """
    load information on how the user wants to run the command from a dictionary of options 
    and info on the files and variables to compare
    note: the options may include a configuration file, which will override many of the
    settings in the options
    """
    
    # basic defaults for stuff we will need to return
    runInfo = {}
    runInfo.update(glance_setting_defaults) # get the default settings
    if ('noLonLatVars' not in optionsSet) or (not optionsSet['noLonLatVars']):
        runInfo.update(glance_lon_lat_defaults) # get the default lon/lat info
    
    # by default, we don't have any particular variables to analyze
    desiredVariables = { }
    # use the built in default values, to start with
    defaultsToUse = glance_analysis_defaults.copy()
    
    requestedNames = None
    
    # set up the paths, they can only come from the command line
    paths = {}
    paths['a']   = aPath
    if bPath is not None:
        paths['b']   = bPath
    paths['out'] = optionsSet['outputpath']
    
    # the colocation selection can only come from the command line options
    # note: since this is really only coming from the user's selection of the call,
    # this is ok for the moment, may want to reconsider later (FUTURE)
    runInfo['doColocate'] = ('doColocate' in optionsSet) and (optionsSet['doColocate'])
    
    # check to see if the user wants to use a config file and if the path exists
    requestedConfigFile = optionsSet['configFile']
    usedConfigFile = False
    
    if (requestedConfigFile is not None) and (requestedConfigFile != "") :
        if not os.path.exists(requestedConfigFile) :
            LOG.warn("Could not open config file: \"" + requestedConfigFile + "\"")
            LOG.warn("Unable to continue analysis without selected configuration file.")
            sys.exit(1)
            
        else :
            
            LOG.info ("Using Config File Settings")
            
            # this will handle relative paths
            requestedConfigFile = os.path.abspath(os.path.expanduser(requestedConfigFile))
            
            # split out the file base name and the file path
            (filePath, fileName) = os.path.split(requestedConfigFile)
            splitFileName = fileName.split('.')
            fileBaseName = fileName[:-3] # remove the '.py' from the end
            
            # hang onto info about the config file for later
            runInfo['config_file_name'] = fileName
            runInfo['config_file_path'] = requestedConfigFile
            
            # load the file
            LOG.debug ('loading config file: ' + str(requestedConfigFile))
            glanceRunConfig = imp.load_module(fileBaseName, file(requestedConfigFile, 'U'),
                                              filePath, ('.py' , 'U', 1))
            
            # this is an exception, since it is not advertised to the user we don't expect it to be in the file
            # (at least not at the moment, it could be added later and if they did happen to put it in the
            # config file, it would override this line)
            runInfo['shouldIncludeReport'] = not optionsSet['imagesOnly'] if 'imagesOnly'   in optionsSet else False
            runInfo['noLonLatVars']        = optionsSet['noLonLatVars']   if 'noLonLatVars' in optionsSet else False
            
            # get everything from the config file
            runInfo.update(glanceRunConfig.settings)
            if ('noLonLatVars' not in runInfo) or (not runInfo['noLonLatVars']) :
                runInfo.update(glanceRunConfig.lat_lon_info) # get info on the lat/lon variables
            
            # get any requested names
            requestedNames = glanceRunConfig.setOfVariables.copy()
            # user selected defaults, if they omit any we'll still be using the program defaults
            defaultsToUse.update(glanceRunConfig.defaultValues)
            
            usedConfigFile = True
    
    # if we didn't get the info from the config file for some reason
    # (the user didn't want to, we couldn't, etc...) get it from the command line options
    if not usedConfigFile:
        
        LOG.info ('Using Command Line Settings')
        
        # so get everything from the options directly
        runInfo['shouldIncludeReport'] = not optionsSet['imagesOnly']
        runInfo['shouldIncludeImages'] = not optionsSet['htmlOnly']
        runInfo['doFork'] = optionsSet['doFork']
        
        # only record these if we are using lon/lat
        runInfo['noLonLatVars']       = optionsSet['noLonLatVars']
        if not runInfo['noLonLatVars'] :
            runInfo['latitude']        = optionsSet['latitudeVar']  or runInfo['latitude']
            runInfo['longitude']       = optionsSet['longitudeVar'] or runInfo['longitude']
            runInfo['lon_lat_epsilon'] = optionsSet['lonlatepsilon'] if 'lonlatepsilon' in optionsSet else None
        
        # get any requested names from the command line
        requestedNames = requestedVars or ['.*'] 
        
        # user selected defaults
        defaultsToUse['epsilon']         = optionsSet['epsilon'] if 'epsilon' in optionsSet else None
        defaultsToUse['missing_value']   = optionsSet['missing']
        
        # note: there is no way to set the tolerances from the command line
    
    return paths, runInfo, defaultsToUse, requestedNames, usedConfigFile

class VariableLoadError(Exception):
    """
    The exception raised when a variable could not be loaded.
    
        msg  -- explanation of which variable could be loaded (and, if possible, why)
    """
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

def _get_variable_from_file(fileObject, variableName, dataType, filter=None) :
    """
    load a variable, using the given data type and applying a filter if one is given
    
    This may throw a VariableLoadError if the variable cannot be loaded.
    """
    
    dataToReturn = None
    exceptionToRaise = None
    
    # get the data from the file
    if fileObject.file_object is None :
        exceptionToRaise = VariableLoadError("File was not properly opened so variable '" + variableName + "' could not be loaded.")
    else :
        try :
            dataToReturn = array(fileObject.file_object[variableName], dtype=dataType)
        except CDFError :
            exceptionToRaise = VariableLoadError('Unable to retrieve ' + variableName + ' data. The variable name ' + 
                      ' may not exist in this file or an error may have occured while attempting to' +
                      ' access the data. Details of file access error observed: ' + str(CDFError))
    
    if (exceptionToRaise is not None) :
        raise exceptionToRaise
    
    if (filter is not None) and (dataToReturn is not None) :
        dataToReturn = filter(dataToReturn)
    
    return dataToReturn

def _get_and_analyze_lon_lat (fileObject,
                              latitudeVariableName, longitudeVariableName,
                              latitudeDataFilterFn=None, longitudeDataFilterFn=None) :
    """
    get the longitude and latitude data from the given file, assuming they are in the given variable names
    and analyze them to identify spacially invalid data (ie. data that would fall off the earth)
    
    This may result in a VariableLoadError if the variable cannot be loaded.
    """
    # get the data from the file
    
    # get the longitude
    LOG.info ('longitude name: ' + longitudeVariableName)
    # TODO, should this dtype be a float?
    longitudeData = _get_variable_from_file(fileObject, longitudeVariableName,
                                            float, filter=longitudeDataFilterFn)
    # get the latitude
    LOG.info ('latitude name: '  + latitudeVariableName)
    # TODO, should this dtype be a float?
    latitudeData  = _get_variable_from_file(fileObject, latitudeVariableName,
                                            float, filter=latitudeDataFilterFn)
    
    # we are going to have issues with our comparision if they aren't the same shape
    LOG.debug('latitude  shape: ' + str(latitudeData.shape))
    LOG.debug('longitude shape: ' + str(longitudeData.shape))
    assert (latitudeData.shape == longitudeData.shape)
    
    # build a mask of our spacially invalid data
    invalidLatitude  = (latitudeData < -90)     | (latitudeData > 90)   | ~isfinite(latitudeData)
    invalidLongitude = (longitudeData < -180)   | (longitudeData > 360) | ~isfinite(longitudeData)
    spaciallyInvalidMask = invalidLatitude | invalidLongitude
    
    # get the missing value as well
    longitudeMissingVal = fileObject.file_object.missing_value(longitudeVariableName)
    latitudeMissingVal  = fileObject.file_object.missing_value( latitudeVariableName)
    
    # analyze our spacially invalid data
    percentageOfSpaciallyInvalidPts, numberOfSpaciallyInvalidPts = _get_percentage_from_mask(spaciallyInvalidMask)
    
    spatialStatInfo = {
                       'totNumInvPts': numberOfSpaciallyInvalidPts,
                       'perInvPts':    percentageOfSpaciallyInvalidPts
                       }
    
    return dataobj.DataObject(longitudeData, fillValue=longitudeMissingVal, ignoreMask=invalidLongitude), \
           dataobj.DataObject(latitudeData,  fillValue=latitudeMissingVal,  ignoreMask=invalidLatitude), spatialStatInfo

def _get_percentage_from_mask(dataMask) :
    """
    given a mask that marks the elements we want the percentage of as True (and is the size of our original data),
    figure out what percentage of the whole they are
    """
    numMarkedDataPts = sum(dataMask)
    totalDataPts = dataMask.size
    
    # avoid dividing by 0
    if totalDataPts is 0 :
        return 0.0, 0
    
    percentage = 100.0 * float(numMarkedDataPts) / float(totalDataPts)
    
    return percentage, numMarkedDataPts

# TODO, this comparison needs to encorporate epsilon percent as well
def _check_lon_lat_equality(longitudeADataObject, latitudeADataObject,
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
    lon_lat_not_equal_mask = longitudeDiffInfo.diff_data_object.masks.mismatch_mask | latitudeDiffInfo.diff_data_object.masks.mismatch_mask
    lon_lat_not_equal_points_count = sum(lon_lat_not_equal_mask)
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
    returnInfo['lon_lat_not_equal_points_count']   = lon_lat_not_equal_points_count
    returnInfo['lon_lat_not_equal_points_percent'] = lon_lat_not_equal_points_percent
    
    return returnInfo

def _compare_spatial_invalidity(longitude_a_object, longitude_b_object,
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
    spatial_info['perInvPtsInBoth'] = spatial_info['file A']['perInvPts']
            # a default that will hold if the two files have the same spatially invalid pts
    if not all(invalid_in_a_mask.ravel() == invalid_in_b_mask.ravel()) : 
        LOG.info("Mismatch in number of spatially invalid points. " +
                 "Files may not have corresponding data where expected.")
        
        # figure out which points are only valid in one of the two files
        valid_only_in_mask_a = (~invalid_in_a_mask) & invalid_in_b_mask
        spatial_info['file A']['numInvPts'] = sum(valid_only_in_mask_a.ravel())
        valid_only_in_mask_b = (~invalid_in_b_mask) & invalid_in_a_mask
        spatial_info['file B']['numInvPts'] = sum(valid_only_in_mask_b.ravel())
        
        # so how many do they have together?
        spatial_info['perInvPtsInBoth'] = _get_percentage_from_mask(invalid_in_common_mask)[0]
        # make a "clean" version of the lon/lat
        longitude_common[valid_only_in_mask_a] = longitude_a_object.data[valid_only_in_mask_a]
        longitude_common[valid_only_in_mask_b] = longitude_b_object.data[valid_only_in_mask_b]
        latitude_common [valid_only_in_mask_a] =  latitude_a_object.data[valid_only_in_mask_a]
        latitude_common [valid_only_in_mask_b] =  latitude_b_object.data[valid_only_in_mask_b]
        
        # plot the points that are only valid one file and not the other
        if ((spatial_info['file A']['numInvPts'] > 0) and (do_include_images) and
            (len(longitude_a_object.data[~invalid_in_a_mask]) > 0) and
            (len( latitude_a_object.data[~invalid_in_a_mask]) > 0)) :
            plot.plot_and_save_spacial_mismatch(longitude_a_object, latitude_a_object,
                                               valid_only_in_mask_a,
                                               "A", "Points only valid in\nFile A\'s longitude & latitude",
                                               "SpatialMismatch",
                                               output_path, True,
                                               fullDPI=fullDPI, thumbDPI=thumbDPI, units="degrees")
        if ((spatial_info['file B']['numInvPts'] > 0) and (do_include_images) and
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

def _load_lon_lat (lon_name, lat_name, file_object, file_descriptior="",
                   alt_file=None, lat_filter=None, lon_filter=None) :
    """
    load the longitude and latitude from a file, filtering as needed
    """
    
    # for the a file, do we have an alternate?
    file_to_use = file_object
    if (alt_file is not None) :
        LOG.info("Loading alternate file (" + alt_file
                 + ") for file " + file_descriptior + " longitude/latitude.")
        file_to_use = dataobj.FileInfo(alt_file)
    
    # load our longitude and latitude and do some analysis on them
    lon_object, lat_object, spatial_info = \
        _get_and_analyze_lon_lat (file_to_use,
                                  lat_name,   lon_name, 
                                  lat_filter, lon_filter)
    
    return lon_object, lat_object, spatial_info

def _handle_lon_lat_info (lon_lat_settings, a_file_object, b_file_object, output_path,
                          should_make_images=False, should_check_equality=True,
                          fullDPI=None, thumbDPI=None) :
    """
    Manage loading and comparing longitude and latitude information for two files
    
    This may result in a VariableLoadError if the longitude or latitude cannot be loaded.
    This may result in a VariableComparisonError if the longitude or latitude cannot be compared due to size.
    
    """
    # a place to save some general stats about our lon/lat data
    spatialInfo = { }
    
    # if there is no lon/lat specified, stop now
    if ( ('longitude' not in lon_lat_settings) or ('latitude' not in lon_lat_settings)
        or (('noLonLatVars' in lon_lat_settings) and lon_lat_settings['noLonLatVars']) ) :
        return { }, spatialInfo
    
    # if we should not be comparing against the logitude and latitude, stop now
    LOG.debug ('lon_lat_settings: ' + str(lon_lat_settings))
    
    # figure out the names to be used for the longitude and latitude variables
    a_longitude_name = lon_lat_settings['longitude']
    a_latitude_name =  lon_lat_settings['latitude']
    b_longitude_name = a_longitude_name
    b_latitude_name =  a_latitude_name
    # if we have alternate b names, use those for b instead
    if ('longitude_alt_name_in_b' in lon_lat_settings) :
        b_longitude_name = lon_lat_settings['longitude_alt_name_in_b']
    if ( 'latitude_alt_name_in_b' in lon_lat_settings):
        b_latitude_name  = lon_lat_settings['latitude_alt_name_in_b']
        
    # if we need to load our lon/lat from different files, open those files
    
    longitude_a_object, latitude_a_object, spatialInfo['file A'] = \
                        _load_lon_lat (a_longitude_name, a_latitude_name, a_file_object, file_descriptior="a",
                                       alt_file=lon_lat_settings['a_lon_lat_from_alt_file']
                                       if ('a_lon_lat_from_alt_file' in lon_lat_settings) else None,
                                       lat_filter=lon_lat_settings['data_filter_function_lat_in_a'],
                                       lon_filter=lon_lat_settings['data_filter_function_lon_in_a'])
    longitude_b_object, latitude_b_object, spatialInfo['file B'] = \
                        _load_lon_lat (b_longitude_name, b_latitude_name, b_file_object, file_descriptior="b",
                                       alt_file=lon_lat_settings['b_lon_lat_from_alt_file']
                                       if ('b_lon_lat_from_alt_file' in lon_lat_settings) else None,
                                       lat_filter=lon_lat_settings['data_filter_function_lat_in_b'],
                                       lon_filter=lon_lat_settings['data_filter_function_lon_in_b'])
    
    # if we need to, test the level of equality of the "valid" values in our lon/lat
    if should_check_equality :
        
        moreSpatialInfo = _check_lon_lat_equality(longitude_a_object, latitude_a_object,
                                                  longitude_b_object, latitude_b_object,
                                                  lon_lat_settings['lon_lat_epsilon'],
                                                  should_make_images, output_path,
                                                  fullDPI=fullDPI, thumbDPI=thumbDPI)
        # update our existing spatial information
        spatialInfo.update(moreSpatialInfo)
        
        # compare our spatially invalid info to see if the two files have invalid longitudes and latitudes in the same places
        spaciallyInvalidMask, spatialInfo, longitude_common, latitude_common = \
                                _compare_spatial_invalidity(longitude_a_object, longitude_b_object,
                                                            latitude_a_object,  latitude_b_object,
                                                            spatialInfo, should_make_images, output_path,
                                                            fullDPI=fullDPI, thumbDPI=thumbDPI)
    else:
        spaciallyInvalidMask = None
        longitude_common     = None
        latitude_common      = None
    
    # FUTURE, return the lon/lat objects instead?
    return {
            'a':      {
                       "lon":       longitude_a_object.data,
                       "lat":       latitude_a_object.data,
                       "inv_mask":  longitude_a_object.masks.ignore_mask,
                       "lon_fill":  longitude_a_object.fill_value,
                       "lat_fill":  latitude_a_object.fill_value
                       },
            'b':      {
                       "lon":       longitude_b_object.data,
                       "lat":       latitude_b_object.data,
                       "inv_mask":  longitude_b_object.masks.ignore_mask,
                       "lon_fill":  longitude_b_object.fill_value,
                       "lat_fill":  latitude_b_object.fill_value
                       },
            'common': {
                       "lon":       longitude_common,
                       "lat":       latitude_common,
                       "inv_mask":  spaciallyInvalidMask
                       }
            }, \
           spatialInfo

def _handle_lon_lat_info_for_one_file (lon_lat_settings, file_object) :
    """
    Manage loading longitude and latitude information for a file
    
    This may result in a VariableLoadError if the longitude or latitude cannot be loaded.
    
    """
    
    # if there is no lon/lat specified, stop now
    if ( ('longitude' not in lon_lat_settings) or ('latitude' not in lon_lat_settings)
        or (('noLonLatVars' in lon_lat_settings) and lon_lat_settings['noLonLatVars']) ) :
        return { }, { }
    
    # print our settings for debugging purposes
    LOG.debug ('lon_lat_settings: ' + str(lon_lat_settings))
    
    # figure out the names to be used for the longitude and latitude variables
    lon_name = lon_lat_settings['longitude']
    lat_name = lon_lat_settings['latitude' ]
    
    # load our lon/lat data
    
    lon_object, lat_object, spatialInfo = \
                        _load_lon_lat (lon_name, lat_name, file_object,
                                       alt_file=lon_lat_settings['a_lon_lat_from_alt_file']
                                       if ('a_lon_lat_from_alt_file' in lon_lat_settings) else None,
                                       lat_filter=lon_lat_settings['data_filter_function_lat_in_a'],
                                       lon_filter=lon_lat_settings['data_filter_function_lon_in_a'])
    
    # FUTURE, return the lon/lat objects instead?
    return {
            "lon":       lon_object.data,
            "lat":       lat_object.data,
            "inv_mask":  lon_object.masks.ignore_mask | lat_object.masks.ignore_mask,
            "lon_fill":  lon_object.fill_value,
            "lat_fill":  lat_object.fill_value
            }, \
           spatialInfo

def _open_and_process_files (args, numFilesExpected):
    """
    open files listed in the args and get information about the variables in them
    """
    # get all the file names
    fileNames = args[:numFilesExpected]
    # open all the files & get their variable names
    files = {}
    commonNames = None
    for fileName in fileNames:
        LOG.info("opening %s" % fileName)
        files[fileName] = {}
        tempFileObject = (io.open(fileName))
        files[fileName]['fileObject'] = tempFileObject
        tempNames = set(tempFileObject())
        LOG.debug ('variable names for ' + fileName + ': ' + str(tempNames)) 
        files[fileName]['varNames'] = tempNames
        if commonNames is None :
            commonNames = tempNames
        else :
            commonNames = commonNames.intersection(tempNames)
    files['commonVarNames'] = commonNames
    
    return files

def _check_pass_or_fail(varRunInfo, variableStats, defaultValues) :
    """
    Check whether the variable passed analysis, failed analysis, or
    did not need to be quantitatively tested
    
    also returns information about the fractions of failure
    """
    
    passValues = [ ]
    
    # test the epsilon value tolerance
    
    # get the tolerance for failures compared to epsilon
    epsilonTolerance = None
    if ('epsilon_failure_tolerance' in varRunInfo) :
        epsilonTolerance = varRunInfo['epsilon_failure_tolerance']
    else :
        epsilonTolerance = defaultValues['epsilon_failure_tolerance']
    
    # did we fail based on the epsilon?
    failed_fraction = variableStats['Numerical Comparison Statistics']['diff_outside_epsilon_fraction']
    passed_epsilon  = None
    if epsilonTolerance is not None :
        passed_epsilon = failed_fraction <= epsilonTolerance
    passValues.append(passed_epsilon)
    
    # test the nonfinite tolerance
    
    # get the tolerance for failures in amount of nonfinite data (in spatially valid areas)
    nonfiniteTolerance = None
    if ('nonfinite_data_tolerance'  in varRunInfo) :
        nonfiniteTolerance = varRunInfo['nonfinite_data_tolerance']
    else :
        nonfiniteTolerance = defaultValues['nonfinite_data_tolerance']
    
    # did we fail based on nonfinite data
    non_finite_diff_fraction = variableStats['Finite Data Statistics']['finite_in_only_one_fraction']
    passed_nonfinite         = None
    if nonfiniteTolerance is not None :
        passed_nonfinite = non_finite_diff_fraction <= nonfiniteTolerance
    passValues.append(passed_nonfinite)
    
    # test if the total failed percentage is acceptable
    
    # get the total percentage of failed data that is acceptable
    totalFailTolerance = None
    if ('total_data_failure_tolerance' in varRunInfo) :
        totalFailTolerance = varRunInfo['total_data_failure_tolerance']
    
    # did we fail based on all data failures?
    passed_all_percentage = None
    if totalFailTolerance is not None :
        passed_all_percentage = (non_finite_diff_fraction + failed_fraction) <= totalFailTolerance
    passValues.append(passed_all_percentage)
    
    # test the r-squared correlation coefficent
    
    # get the minimum acceptable r-squared correlation coefficient
    min_r_squared = None
    if ('minimum_acceptable_squared_correlation_coefficient' in varRunInfo) :
        min_r_squared = varRunInfo['minimum_acceptable_squared_correlation_coefficient']
    else :
        min_r_squared = defaultValues['minimum_acceptable_squared_correlation_coefficient']
    
    # did we fail based on the r-squared correlation coefficient?
    r_squared_value  = None
    passed_r_squared = None
    if min_r_squared is not None :
        r_squared_value  = variableStats['Numerical Comparison Statistics']['r-squared correlation']
        passed_r_squared = r_squared_value >= min_r_squared
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

def _get_run_identification_info( ) :
    """
    get info about what user/machine/version of glance is being used
    """
    info_to_return = { }
    
    # get info on who's doing the run and where
    info_to_return['machine'] = os.uname()[1]      # the name of the machine running the report
    info_to_return['user'] = os.getenv("LOGNAME")  #os.getlogin() # the name of the user running the report
    info_to_return['version'] = _get_glance_version_string()
    
    return info_to_return

def _get_glance_version_string() :
    version_num = pkg_resources.require('uwglance')[0].version
    
    return "glance, version " + str(version_num) 

def _get_name_info_for_variable(original_display_name, variable_run_info) :
    """
    based on the variable run info, figure out the various names for
    the variable and return them
    
    the various names are:
    
    technical_name -            the name the variable is listed under in the file
    b_variable_technical_name - the name the variable is listed under in the b file (may be the same as technical_name)
    explanation_name -          the more verbose name that will be shown to the user to identify the variable
    original_display_name -     the display name given by the user to describe the variable
    """
    
    # figure out the various name related info
    technical_name = variable_run_info['variable_name']
    explanation_name = technical_name # for now, will add to this later
    
    # if B has an alternate variable name, figure that out
    b_variable_technical_name = technical_name
    if 'alternate_name_in_B' in variable_run_info :
        b_variable_technical_name = variable_run_info['alternate_name_in_B']
        # put both names in our explanation
        explanation_name = explanation_name + " / " + b_variable_technical_name
    
    # show both the display and current explanation names if they differ
    if not (original_display_name == explanation_name) :
        explanation_name = original_display_name + ' (' + explanation_name + ')'
    
    return technical_name, b_variable_technical_name, explanation_name

def _load_variable_data(fileObject, variableNameInFile,
                        dataFilter=None,
                        variableToFilterOn=None,
                        variableBasedFilter=None,
                        fileDescriptionForDisplay="file") :
    """
    load data for a variable from a file
    optionally filter the variable data based on a data filter or another variable
    
    dataFilter must be in the form of (lambda data: some manipulation returning the new data)
    variableBasedFilter must be in the form of (lambda data, filterData: some manipulation returning the new data))
    """
    
    # get the data for the variable
    LOG.debug("loading basic data for variable " + variableNameInFile + " from " + fileDescriptionForDisplay)
    variableData = fileObject[variableNameInFile]
    
    # apply the basic filter if there is one
    if dataFilter is not None :
        LOG.debug ("applying filter function to data from " + fileDescriptionForDisplay + " for variable " + variableNameInFile)
        variableData = dataFilter(variableData)
    
    # if we've got another variable to filter on, do that
    if (variableToFilterOn is not None) and (variableBasedFilter is not None) :
        LOG.debug ("filtering data from " + fileDescriptionForDisplay + " for variable " + variableNameInFile
                   + " based on additional data from variable " + variableToFilterOn)
        dataToFilterOn = fileObject[variableToFilterOn]
        variableData = variableBasedFilter(variableData, dataToFilterOn)
    
    return variableData

def _uri_needs_rsync(uri_to_check) :
    """
    check if the uri requires an rsync in order to access the data
    this will return some false positives if you phrase local uri's with the machine name
    for ex. you are on the machine "lotus" and you use the path "rsync:://lotus/data/"
    """
    return not os.path.exists(uri_to_check)

def _get_UV_info_from_magnitude_direction_info(fileObject, magnitudeName, directionName, invalidMask=None) :
    """
    If there are magnitude and direction names, load that information and calculate the u and v that correspond to it
    """
    
    # if we don't have magnitude and direction, we can't calculate the U and V values
    if (magnitudeName is None) or (directionName is None) :
        return None, None
    
    # load the magnitude and direction data sets
    magnitude = _load_variable_data(fileObject, magnitudeName)
    direction = _load_variable_data(fileObject, directionName)
    
    # convert the magnitude and direction data into u and v vectors
    uData, vData = delta.convert_mag_dir_to_U_V_vector(magnitude, direction, invalidMask=invalidMask)
    
    return uData, vData

def rsync_or_copy_files (list_of_files, target_directory='.', additionalFileNameSuffix='') :
    """
    If the files in the list are remote, rsync them, otherwise, just copy
    them to the target directory
    """
    newPaths = [ ]
    
    for file_uri in list_of_files :
        fileName = os.path.split(file_uri)[1]
        baseFile, ext = os.path.splitext(fileName)
        newPath = os.path.join(target_directory, baseFile + additionalFileNameSuffix + ext)
        newPaths.append(newPath)
        
        if _uri_needs_rsync(file_uri) :
            cmd = ['rsync', '-Cuav', file_uri, newPath]
        else :
            cmd = ['cp', os.path.abspath(file_uri), newPath]
        LOG.debug('running ' + ' '.join(cmd)) 
        sh(cmd)
    
    return newPaths

def colocateToFile_library_call(a_path, b_path, var_list=[ ],
                                options_set={ },
                                # todo, this doesn't yet do anything
                                do_document=False,
                                # todo, the output channel does nothing at the moment
                                output_channel=sys.stdout) :
    """
    this method handles the actual work of the colocateData command line tool
    and can be used as a library routine.
    
    TODO, properly document the options
    """
    
    # load the user settings from either the command line or a user defined config file
    pathsTemp, runInfo, defaultValues, requestedNames, usedConfigFile = _load_config_or_options(a_path, b_path,
                                                                                                options_set,
                                                                                                requestedVars = var_list)
    
    # deal with the input and output files
    if not (os.path.isdir(pathsTemp['out'])) :
        LOG.info("Specified output directory (" + pathsTemp['out'] + ") does not exist.")
        LOG.info("Creating output directory.")
        os.makedirs(pathsTemp['out'])
    
    # make copies of the input files for colocation TODO, fix paths
    [pathsTemp['a'], pathsTemp['b']] = rsync_or_copy_files ([pathsTemp['a'], pathsTemp['b']],
                                                            target_directory=pathsTemp['out'],
                                                            additionalFileNameSuffix='-collocated')
    
    # open the files
    LOG.info("Processing File A:")
    aFile = dataobj.FileInfo(pathsTemp['a'], allowWrite=True)
    if aFile is None:
        LOG.warn("Unable to continue with comparison because file a (" + pathsTemp['a'] + ") could not be opened.")
        sys.exit(1)
    LOG.info("Processing File B:")
    bFile = dataobj.FileInfo(pathsTemp['b'], allowWrite=True)
    if bFile is None:
        LOG.warn("Unable to continue with comparison because file b (" + pathsTemp['b'] + ") could not be opened.")
        sys.exit(1)
    
    # get information about the names the user requested
    finalNames, nameStats = _resolve_names(aFile.file_object, bFile.file_object,
                                           defaultValues,
                                           requestedNames, usedConfigFile)
    
    # return for lon_lat_data variables will be in the form 
    #{"lon": longitude_data,      "lat": latitude_data,      "inv_mask": spaciallyInvalidMaskData}
    # or { } if there is no lon/lat info
    lon_lat_data = { }
    try :
        lon_lat_data, _ = _handle_lon_lat_info (runInfo, aFile, bFile, pathsTemp['out'], should_check_equality=False,
                                                fullDPI=runInfo['detail_DPI'], thumbDPI=runInfo['thumb_DPI'])
    except VariableLoadError, vle :
        LOG.warn("Error while loading longitude or latitude: ")
        LOG.warn(vle.msg)
        exit(1)
    except VariableComparisonError, vce :
        LOG.warn("Error while comparing longitude or latitude: ")
        LOG.warn(vce.msg)
        exit(1)
    
    # handle the longitude and latitude colocation
    LOG.info("Colocating raw longitude and latitude information")
    aColocationInfomation, bColocationInformation, totalNumberOfMatchedPoints = \
                    collocation.create_colocation_mapping_within_epsilon((lon_lat_data['a']['lon'], lon_lat_data['a']['lat']),
                                                                         (lon_lat_data['b']['lon'], lon_lat_data['b']['lat']),
                                                                         runInfo['lon_lat_epsilon'],
                                                                         invalidAMask=lon_lat_data['a']['inv_mask'],
                                                                         invalidBMask=lon_lat_data['b']['inv_mask'])
    (colocatedLongitude, colocatedLatitude, (numMultipleMatchesInA, numMultipleMatchesInB)), \
    (unmatchedALongitude, unmatchedALatitude), \
    (unmatchedBLongitude, unmatchedBLatitude) = \
                collocation.create_colocated_lonlat_with_lon_lat_colocation(aColocationInfomation, bColocationInformation,
                                                                            totalNumberOfMatchedPoints,
                                                                            lon_lat_data['a']['lon'], lon_lat_data['a']['lat'],
                                                                            lon_lat_data['b']['lon'], lon_lat_data['b']['lat'])
    
    # TODO, based on unmatched, issue warnings and record info in the file?
    LOG.debug("colocated shape of the longitude: " + str(colocatedLongitude.shape))
    LOG.debug("colocated shape of the latitude:  " + str(colocatedLatitude.shape))
    LOG.debug(str(numMultipleMatchesInA) + " lon/lat pairs contain A points used for multiple matches.")
    LOG.debug(str(numMultipleMatchesInB) + " lon/lat pairs contain B points used for multiple matches.")
    LOG.debug(str(len(unmatchedALatitude)) + " A lon/lat points could not be matched.")
    LOG.debug(str(len(unmatchedBLatitude)) + " B lon/lat points could not be matched.")
    
    # go through each of the possible variables in our files
    # and do our colocation for whichever ones we can
    for displayName in finalNames:
        
        # pull out the information for this variable analysis run
        varRunInfo = finalNames[displayName].copy()
        
        # get the various names
        technical_name, b_variable_technical_name, \
                explanationName = _get_name_info_for_variable(displayName, varRunInfo)
        
        LOG.info('analyzing: ' + explanationName + ')')
        
        # load the variable data
        aData = _load_variable_data(aFile.file_object, technical_name,
                                    dataFilter = varRunInfo['data_filter_function_a'] if 'data_filter_function_a' in varRunInfo else None,
                                    variableToFilterOn = varRunInfo['variable_to_filter_on_a'] if 'variable_to_filter_on_a' in varRunInfo else None,
                                    variableBasedFilter = varRunInfo['variable_based_filter_a'] if 'variable_based_filter_a' in varRunInfo else None,
                                    fileDescriptionForDisplay = "file A")
        bData = _load_variable_data(bFile.file_object, b_variable_technical_name,
                                    dataFilter = varRunInfo['data_filter_function_b'] if 'data_filter_function_b' in varRunInfo else None,
                                    variableToFilterOn = varRunInfo['variable_to_filter_on_b'] if 'variable_to_filter_on_b' in varRunInfo else None,
                                    variableBasedFilter = varRunInfo['variable_based_filter_b'] if 'variable_based_filter_b' in varRunInfo else None,
                                    fileDescriptionForDisplay = "file B")
        
        # colocate the data for this variable if we have longitude/latitude data
        if (len(lon_lat_data.keys()) > 0) and runInfo['doColocate'] :
            
            # figure out the invalid masks
            invalidA = lon_lat_data['a']['inv_mask'] | (aData == varRunInfo['missing_value'])
            invalidB = lon_lat_data['b']['inv_mask'] | (bData == varRunInfo['missing_value_alt_in_b'])
            
            # match up our points in A and B
            (aData, bData, (numberOfMultipleMatchesInA, numberOfMultipleMatchesInB)), \
            (aUnmatchedData,             unmatchedALongitude, unmatchedALatitude), \
            (bUnmatchedData,             unmatchedBLongitude, unmatchedBLatitude) = \
                    collocation.create_colocated_data_with_lon_lat_colocation(aColocationInfomation, bColocationInformation,
                                                                              colocatedLongitude, colocatedLatitude,
                                                                              aData, bData,
                                                                              missingData=varRunInfo['missing_value'],
                                                                              altMissingDataInB=varRunInfo['missing_value_alt_in_b'],
                                                                              # TODO, should missing data be considered?
                                                                              invalidAMask=invalidA,
                                                                              invalidBMask=invalidB)
            
            LOG.debug(str(numberOfMultipleMatchesInA) + " data pairs contain A data points used for multiple matches.")
            LOG.debug(str(numberOfMultipleMatchesInB) + " data pairs contain B data points used for multiple matches.")
            LOG.debug(str(len(aUnmatchedData)) + " A data points could not be matched.")
            LOG.debug(str(len(bUnmatchedData)) + " B data points could not be matched.")
            
            # save the colocated data information in the output files
            
            # all the a file information
            aFile.file_object.create_new_variable(technical_name + '-colocated', # TODO, how should this suffix be handled?
                                      missingvalue = varRunInfo['missing'] if 'missing' in varRunInfo else None,
                                      data = aData,
                                      variabletocopyattributesfrom = technical_name)
            aFile.file_object.add_attribute_data_to_variable(technical_name + '-colocated', 'number of multiple matches', numberOfMultipleMatchesInA)
            aFile.file_object.add_attribute_data_to_variable(technical_name + '-colocated', 'number of unmatched points', len(aUnmatchedData))
            
            # all the b file information
            bFile.file_object.create_new_variable(b_variable_technical_name + '-colocated', # TODO, how should this suffix be handled?
                                      missingvalue = varRunInfo['missing_value_alt_in_b'] if 'missing_value_alt_in_b' in varRunInfo else None,
                                      data = bData,
                                      variabletocopyattributesfrom = b_variable_technical_name)
            bFile.file_object.add_attribute_data_to_variable(b_variable_technical_name + '-colocated', 'number of multiple matches', numberOfMultipleMatchesInB)
            bFile.file_object.add_attribute_data_to_variable(b_variable_technical_name + '-colocated', 'number of unmatched points', len(bUnmatchedData))
            
            # TODO, any additional statistics
            
        else :
            LOG.debug(explanationName + " was not selected for colocation and will be ignored.")
        
    # the end of the loop to examine all the variables
    
    # we're done with the files, so close them up
    aFile.file_object.close()
    bFile.file_object.close()
    
    return

def reportGen_raw_data_simple_call (aData, bData, variableDisplayName,
                                    epsilon=0.0, missingValue=None,
                                    useThreads=True, includeImages=True,
                                    outputDirectory="./") :
    """
    Generate a report for a single variable given raw data and
    some minimal control settings. This method will also generate
    images for the report if includeImages is True.
    """
    
    LOG.info("Setting up basic information")
    
    aData = array(aData)
    bData = array(bData)
    
    # set up the run info
    runInfo = { }
    runInfo.update(glance_setting_defaults)
    runInfo['shouldIncludeImages'] = True
    runInfo['shouldIncludeReport'] = True
    runInfo['doFork']              = False
    runInfo['useThreadsToControlMemory'] = useThreads
    
    # set up the variable specific info
    variableSettings = { }
    variableSettings = glance_analysis_defaults.copy()
    variableSettings['epsilon']                = epsilon
    variableSettings['missing_value']          = missingValue
    variableSettings['missing_value_alt_in_b'] = missingValue
    variableSettings['variable_name']          = variableDisplayName
    
    # hang onto identification info
    runInfo.update(_get_run_identification_info( ))
    
    # deal with the output directories
    outputDirectory = _clean_path(outputDirectory)
    _setup_dir_if_needed(outputDirectory, "output")
    
    LOG.info("Analyzing " + variableDisplayName)
    
    # if things are the same shape, analyze them and make our images
    if aData.shape == bData.shape :
        
        # setup some values in the variable settings for use in the report
        variableSettings['variable_dir'] = outputDirectory
        variableSettings['variable_report_path_escaped'] = quote(os.path.join(variableDisplayName, 'index.html'))
        variableSettings['doc_path'] = quote(os.path.join(outputDirectory, './' + 'doc.html')) 
        
        # calculate the variable statistics
        variable_stats = statistics.StatisticalAnalysis.withSimpleData(aData, bData,
                                                                       missingValue, missingValue,
                                                                       None, None,
                                                                       epsilon, None).dictionary_form()
        
        # add a little additional info
        variableSettings['time'] = datetime.datetime.ctime(datetime.datetime.now())
        didPass, epsilon_failed_fraction, \
            non_finite_fail_fraction, r_squared_value = _check_pass_or_fail(variableSettings, variable_stats, glance_setting_defaults)
        variableSettings['did_pass'] = didPass
        
        # to hold the names of any images created
        image_names = {
                        'original': [ ],
                        'compared': [ ]
                        }
        
        # if we need the images, make them now
        if includeImages :
            
            LOG.info("Plotting images for " + variableDisplayName)
            
            plotFunctionGenerationObjects = [ ]
            
            # add the function to make the histogram and scatter plot
            plotFunctionGenerationObjects.append(plotcreate.BasicComparisonPlotsFunctionFactory())
            
            # add the function to do basic imshow images
            plotFunctionGenerationObjects.append(plotcreate.IMShowPlotFunctionFactory())
            
            # plot our lon/lat related info
            image_names['original'], image_names['compared'] = \
                plot.plot_and_save_comparison_figures \
                        (aData, bData,
                         plotFunctionGenerationObjects,
                         outputDirectory,
                         variableDisplayName,
                         epsilon,
                         missingValue,
                         lonLatDataDict=None,
                         makeSmall=True,
                         doFork=False,
                         shouldClearMemoryWithThreads=useThreads,
                         shouldUseSharedRangeForOriginal=True)
            
            LOG.info("\tfinished creating figures for: " + variableDisplayName)
        
        # create a temporary files object
        files = {'file A': {'path': "raw data input",
                            'lastModifiedTime': "unknown",
                            'md5sum': "n/a"
                            },
                 'file B': {'path': "raw data input",
                            'lastModifiedTime': "unknown",
                            'md5sum': "n/a"
                            }
                    }
        
        # create our report 
        LOG.info ('Generating report for: ' + variableDisplayName) 
        report.generate_and_save_variable_report(files,
                                                 variableSettings, runInfo,
                                                 variable_stats,
                                                 { },
                                                 image_names,
                                                 outputDirectory, "index.html")
        
        # make the glossary page
        LOG.info ('Generating glossary page')
        report.generate_and_save_doc_page(statistics.StatisticalAnalysis.doc_strings(), outputDirectory)
        
    else :
        message = (variableDisplayName + ' ' + 
                'could not be compared. This may be because the data for this variable does not match in shape ' +
                'between the two files (file A data shape: ' + str(aData.shape) + '; file B data shape: '
                + str(bData.shape) + ').')
        LOG.warn(message)

def _setup_dir_if_needed(dirPath, descriptionName) :
    """
    create the directory if that is needed, if not don't
    """
    if not (os.path.isdir(dirPath)) :
        LOG.info("Specified " + descriptionName + " directory (" + dirPath + ") does not exist.")
        LOG.info("Creating " + descriptionName + " directory.")
        os.makedirs(dirPath)

def inspect_library_call (a_path, var_list=[ ],
                          options_set={ },
                          # todo, this doesn't yet do anything
                          do_document=False,
                          # todo, the output channel does nothing at the moment
                          output_channel=sys.stdout) :
    """
    this method handles the actual work of the inspectReport command line tool
    and can also be used as a library routine, pass in the slightly parsed
    command line input, or call it as a library function... be sure to fill
    out the options
    
    TODO at the moment the options are very brittle and need to be fully filled
    or this method will fail badly (note: the addition of some glance defaults
    has minimized the problem, but you still need to be careful when dealing with
    optional boolean values. this needs more work.)
    """
    
    # load the user settings from either the command line or a user defined config file
    pathsTemp, runInfo, defaultValues, requestedNames, usedConfigFile = _load_config_or_options(a_path, None, # there is no B path
                                                                                                options_set,
                                                                                                requestedVars = var_list)
    
    # information for debugging purposes
    LOG.debug('paths: ' +           str(pathsTemp))
    LOG.debug('defaults: ' +        str(defaultValues))
    LOG.debug('run information: ' + str(runInfo))
    
    # if we wouldn't generate anything, just stop now
    if (not runInfo['shouldIncludeImages']) and (not runInfo['shouldIncludeReport']) :
        LOG.warn("User selection of no image generation and no report generation will result in no " +
                 "content being generated. Aborting generation function.")
        return
    
    # hang onto info to identify who/what/when/where/etc. the report is being run by/for 
    runInfo.update(_get_run_identification_info( ))
    
    # deal with the input and output files
    if not (os.path.isdir(pathsTemp['out'])) :
        LOG.info("Specified output directory (" + pathsTemp['out'] + ") does not exist.")
        LOG.info("Creating output directory.")
        os.makedirs(pathsTemp['out'])
    # open the file
    files = {}
    LOG.info("Processing File A:")
    aFile = dataobj.FileInfo(pathsTemp['a'])
    files['file A'] = aFile.get_old_info_dictionary() # FUTURE move to actually using the file object to generate the report
    if aFile.file_object is None:
        LOG.warn("Unable to continue with examination because file (" + pathsTemp['a'] + ") could not be opened.")
        sys.exit(1)
    
    # get information about the names the user requested
    nameStats = {}
    finalNames, nameStats["possibleNames"] = _resolve_names_one_file(aFile.file_object,
                                                                     defaultValues, # TODO, might need a different default set
                                                                     requestedNames,
                                                                     usedConfigFile)
    
    LOG.debug("output dir: " + str(pathsTemp['out']))
    
    # return for lon_lat_data variables will be in the form 
    #{"lon": longitude_data,      "lat": latitude_data,      "inv_mask": spaciallyInvalidMaskData}
    # or { } if there is no lon/lat info
    lon_lat_data = { }
    spatialInfo  = { }
    try :
        lon_lat_data, spatialInfo = _handle_lon_lat_info_for_one_file (runInfo, aFile)
    except VariableLoadError, vle :
        LOG.warn("Error while loading longitude or latitude: ")
        LOG.warn(vle.msg)
        exit(1)
    
    # if there is an approved lon/lat shape, hang on to that for future variable data shape checks
    good_shape_from_lon_lat = None
    if len(lon_lat_data.keys()) > 0:
        good_shape_from_lon_lat = lon_lat_data['lon'].shape
    
    # go through each of the possible variables in our files
    # and make a report section with images for whichever ones we can
    variableInspections = { }
    for displayName in finalNames:
        
        # pull out the information for this variable analysis run
        varRunInfo = finalNames[displayName].copy()
        
        # get the various names
        technical_name, _, explanationName = _get_name_info_for_variable(displayName, varRunInfo)
        
        LOG.info('analyzing: ' + explanationName)
        
        # load the variable data
        aData = _load_variable_data(aFile.file_object, technical_name,
                                    dataFilter = varRunInfo['data_filter_function_a'] if 'data_filter_function_a' in varRunInfo else None,
                                    variableToFilterOn = varRunInfo['variable_to_filter_on_a'] if 'variable_to_filter_on_a' in varRunInfo else None,
                                    variableBasedFilter = varRunInfo['variable_based_filter_a'] if 'variable_based_filter_a' in varRunInfo else None,
                                    fileDescriptionForDisplay = "file A")
        
        # pre-check if this data should be plotted and if it should be compared to the longitude and latitude
        include_images_for_this_variable = ((not('shouldIncludeImages' in runInfo)) or (runInfo['shouldIncludeImages']))
        if 'shouldIncludeImages' in varRunInfo :
            include_images_for_this_variable = varRunInfo['shouldIncludeImages']
        do_not_test_with_lon_lat = (not include_images_for_this_variable) or (len(lon_lat_data.keys()) <= 0)
        
        # handle vector data
        isVectorData = ('magnitudeName' in varRunInfo)  and ('directionName'  in varRunInfo)
        
        # check if this data can be examined 
        # (don't compare lon/lat sizes if we won't be plotting)
        if ( do_not_test_with_lon_lat or (aData.shape == good_shape_from_lon_lat) ) :
            
            # check to see if there is a directory to put information about this variable in,
            # if not then create it
            variableDir = os.path.join(pathsTemp['out'], './' + displayName)
            varRunInfo['variable_dir'] = variableDir
            varRunInfo['variable_report_path_escaped'] = quote(os.path.join(displayName, 'index.html'))
            LOG.debug ("Directory selected for variable information: " + varRunInfo['variable_report_path_escaped'])
            if not (os.path.isdir(variableDir)) :
                LOG.debug("Variable directory (" + variableDir + ") does not exist.")
                LOG.debug("Creating variable directory.")
                os.makedirs(variableDir)
            
            # form the doc and config paths relative to where the variable is
            upwardPath = './'
            for number in range(len(displayName.split('/'))) : # TODO this is not general to windows
                upwardPath = os.path.join(upwardPath, '../')
            varRunInfo['doc_path'] = quote(os.path.join(upwardPath, 'doc.html'))
            if 'config_file_name' in runInfo :
                varRunInfo['config_file_path'] = quote(os.path.join(upwardPath, runInfo['config_file_name']))
            
            #print ("*** lon lat data temp ***")
            #print (str(lon_lat_data))
            
            # figure out the masks we want, and then do our statistical analysis
            mask_a_to_use = None
            if not do_not_test_with_lon_lat :
                mask_a_to_use = lon_lat_data['inv_mask']
            
            variable_stats = statistics.StatisticalInspectionAnalysis.withSimpleData(aData,
                                                                                     missingValue=varRunInfo['missing_value'],
                                                                                     ignoreMask=mask_a_to_use).dictionary_form()
            
            # add a little additional info to our variable run info before we squirrel it away
            varRunInfo['time'] = datetime.datetime.ctime(datetime.datetime.now())  # todo is this needed?
            
            # to hold the names of any images created
            image_names = {
                            'original': [ ],
                            'compared': [ ]
                            }
            
            # create the images for this variable
            if (include_images_for_this_variable) :
                
                plotFunctionGenerationObjects = [ ]
                
                # we are always going to want to draw a basic histogram of the data values to tell which
                # occur most frequently
                plotFunctionGenerationObjects.append(plotcreate.DataHistogramPlotFunctionFactory())
                
                # if it's vector data with longitude and latitude, quiver plot it on the Earth
                if isVectorData and (not do_not_test_with_lon_lat) :
                    # TODO replace this at some point
                    #plotFunctionGenerationObjects.append(plotcreate.MappedQuiverPlotFunctionFactory())
                    pass
                
                # if the data is one dimensional we can plot it as lines
                elif   (len(aData.shape) is 1) : 
                    plotFunctionGenerationObjects.append(plotcreate.InspectLinePlotsFunctionFactory())
                
                # if the data is 2D we have some options based on the type of data
                elif (len(aData.shape) is 2) :
                    
                    # if the data is not mapped to a longitude and latitude, just show it as an image
                    if (do_not_test_with_lon_lat) :
                        plotFunctionGenerationObjects.append(plotcreate.InspectIMShowPlotFunctionFactory())
                    
                    # if it's 2D and mapped to the Earth, contour plot it on the earth
                    else :
                        plotFunctionGenerationObjects.append(plotcreate.InspectMappedContourPlotFunctionFactory())
                
                # if there's magnitude and direction data, figure out the u and v, otherwise these will be None
                aUData, aVData = _get_UV_info_from_magnitude_direction_info (aFile.file_object,
                                                                             varRunInfo['magnitudeName'] if ('magnitudeName') in varRunInfo else None,
                                                                             varRunInfo['directionName'] if ('directionName') in varRunInfo else None,
                                                                             lon_lat_data['inv_mask']
                                                                             if ('inv_mask' in lon_lat_data) else None)
                
                # plot our images
                image_names['original'], image_names['compared'] = \
                    plot.plot_and_save_comparison_figures \
                            (aData, None, # there is no b data
                             plotFunctionGenerationObjects,
                             varRunInfo['variable_dir'],
                             displayName,
                             None, # there is no epsilon
                             varRunInfo['missing_value'],
                             lonLatDataDict=lon_lat_data,
                             dataRanges     = varRunInfo['display_ranges']      if 'display_ranges'      in varRunInfo else None,
                             dataRangeNames = varRunInfo['display_range_names'] if 'display_range_names' in varRunInfo else None,
                             dataColors     = varRunInfo['display_colors']      if 'display_colors'      in varRunInfo else None,
                             makeSmall=True,
                             doFork=runInfo['doFork'],
                             shouldClearMemoryWithThreads=runInfo['useThreadsToControlMemory'],
                             shouldUseSharedRangeForOriginal=runInfo['useSharedRangeForOriginal'],
                             doPlotSettingsDict = varRunInfo,
                             aUData=aUData, aVData=aVData,
                             fullDPI=       runInfo['detail_DPI'],
                             thumbDPI=      runInfo['thumb_DPI'],
                             units_a=       varRunInfo['units_a']               if 'units_a'             in varRunInfo else None,
                             useBData=False)
                
                LOG.info("\tfinished creating figures for: " + explanationName)
            
            # create the report page for this variable
            if (runInfo['shouldIncludeReport']) :
                
                # hang on to some info on our variable
                variableInspections[displayName] = {
                                                    'variable_run_info': varRunInfo
                                                    }
                
                LOG.info ('\tgenerating report for: ' + explanationName) 
                report.generate_and_save_inspect_variable_report(files, varRunInfo, runInfo,
                                                                 variable_stats, spatialInfo, image_names,
                                                                 varRunInfo['variable_dir'], "index.html")
        
        # if we can't do anything with the variable, we should tell the user 
        else :
            message = (explanationName + ' could not be examined. '
                     + 'This may be because the data for this variable (data shape: '
                     + str(aData.shape) + ') does not match the shape of the selected '
                     + 'longitude ' + str(good_shape_from_lon_lat) + ' and '
                     + 'latitude '  + str(good_shape_from_lon_lat) + ' variables.')
            LOG.warn(message)
        
    # the end of the loop to examine all the variables
    
    # generate our general report pages once we've analyzed all the variables
    if (runInfo['shouldIncludeReport']) :
        
        # get the current time
        runInfo['time'] = datetime.datetime.ctime(datetime.datetime.now())
        
        print("*** name stats: " + str (nameStats))
        
        # TODO, create a new report generation function here
        # make the main summary report
        LOG.info ('generating summary report')
        report.generate_and_save_inspection_summary_report (files,
                                                            pathsTemp['out'], 'index.html',
                                                            runInfo,
                                                            variableInspections,
                                                            spatialInfo,
                                                            nameStats)
        
        # make the glossary
        LOG.info ('generating glossary')
        report.generate_and_save_doc_page(statistics.StatisticalInspectionAnalysis.doc_strings(), pathsTemp['out'])
    
    return 0

def reportGen_library_call (a_path, b_path, var_list=[ ],
                            options_set={ },
                            # todo, this doesn't yet do anything
                            do_document=False,
                            # todo, the output channel does nothing at the moment
                            output_channel=sys.stdout) :
    """
    this method handles the actual work of the reportGen command line tool
    and can also be used as a library routine, pass in the slightly parsed
    command line input, or call it as a library function... be sure to fill
    out the options
    TODO at the moment the options are very brittle and need to be fully filled
    or this method will fail badly (note: the addition of some glance defaults
    has minimized the problem, but you still need to be careful when dealing with
    optional boolean values. this needs more work.)
    """
    
    # have all the variables passed test criteria set for them?
    # if no criteria were set then this will be true
    didPassAll = True
    do_pass_fail = options_set['usePassFail'] # todo, this is a temporary hack, should be loaded with other options
    
    # load the user settings from either the command line or a user defined config file
    pathsTemp, runInfo, defaultValues, requestedNames, usedConfigFile = _load_config_or_options(a_path, b_path,
                                                                                                options_set,
                                                                                                requestedVars = var_list)
    
    # note some of this information for debugging purposes
    LOG.debug('paths: ' +           str(pathsTemp))
    LOG.debug('defaults: ' +        str(defaultValues))
    LOG.debug('run information: ' + str(runInfo))
    
    # if we wouldn't generate anything, just stop now
    if (not runInfo['shouldIncludeImages']) and (not runInfo['shouldIncludeReport']) :
        LOG.warn("User selection of no image generation and no report generation will result in no " +
                 "content being generated. Aborting generation function.")
        if do_pass_fail :
            return 0 # nothing went wrong, we just had nothing to do!
        else :
            return
    
    # hang onto info to identify who/what/when/where/etc. the report is being run by/for 
    runInfo.update(_get_run_identification_info( ))
    
    # deal with the input and output files
    if not (os.path.isdir(pathsTemp['out'])) :
        LOG.info("Specified output directory (" + pathsTemp['out'] + ") does not exist.")
        LOG.info("Creating output directory.")
        os.makedirs(pathsTemp['out'])
    # open the files
    files = {}
    LOG.info("Processing File A:")
    aFile = dataobj.FileInfo(pathsTemp['a'])
    files['file A'] = aFile.get_old_info_dictionary() # FUTURE move to actually using the file object to generate the report
    if aFile.file_object is None:
        LOG.warn("Unable to continue with comparison because file a (" + pathsTemp['a'] + ") could not be opened.")
        sys.exit(1)
    LOG.info("Processing File B:")
    bFile = dataobj.FileInfo(pathsTemp['b']) 
    files['file B'] = bFile.get_old_info_dictionary() # FUTURE move to actually using the file object to generate the report
    if bFile.file_object is None:
        LOG.warn("Unable to continue with comparison because file b (" + pathsTemp['b'] + ") could not be opened.")
        sys.exit(1)
    
    # get information about the names the user requested
    finalNames, nameStats = _resolve_names(aFile.file_object, bFile.file_object,
                                           defaultValues,
                                           requestedNames, usedConfigFile)
    
    LOG.debug("output dir: " + str(pathsTemp['out']))
    
    # return for lon_lat_data variables will be in the form 
    #{"lon": longitude_data,      "lat": latitude_data,      "inv_mask": spaciallyInvalidMaskData}
    # or { } if there is no lon/lat info
    lon_lat_data = { }
    spatialInfo  = { }
    try :
        lon_lat_data, spatialInfo = _handle_lon_lat_info (runInfo, aFile, bFile, pathsTemp['out'],
                                                          should_make_images = runInfo["shouldIncludeImages"],
                                                          fullDPI=runInfo['detail_DPI'], thumbDPI=runInfo['thumb_DPI'])
    except VariableLoadError, vle :
        LOG.warn("Error while loading longitude or latitude: ")
        LOG.warn(vle.msg)
        exit(1)
    except VariableComparisonError, vce :
        LOG.warn("Error while comparing longitude or latitude: ")
        LOG.warn(vce.msg)
        exit(1)
    
    # if there is an approved lon/lat shape, hang on to that for future checks
    good_shape_from_lon_lat = None
    if len(lon_lat_data.keys()) > 0:
        good_shape_from_lon_lat = lon_lat_data['common']['lon'].shape
    
    # this will hold information for the summary report
    # it will be in the form
    # [displayName] = {"passEpsilonPercent":     percent ok with epsilon,
    #                  "finite_similar_percent": percent with the same finiteness, 
    #                  "epsilon":                epsilon value used}
    variableComparisons = {}
    
    # go through each of the possible variables in our files
    # and make a report section with images for whichever ones we can
    for displayName in finalNames:
        
        # pull out the information for this variable analysis run
        varRunInfo = finalNames[displayName].copy()
        
        # get the various names
        technical_name, b_variable_technical_name, \
                explanationName = _get_name_info_for_variable(displayName, varRunInfo)
        
        LOG.info('analyzing: ' + explanationName)
        
        # load the variable data
        aData = _load_variable_data(aFile.file_object, technical_name,
                                    dataFilter = varRunInfo['data_filter_function_a'] if 'data_filter_function_a' in varRunInfo else None,
                                    variableToFilterOn = varRunInfo['variable_to_filter_on_a'] if 'variable_to_filter_on_a' in varRunInfo else None,
                                    variableBasedFilter = varRunInfo['variable_based_filter_a'] if 'variable_based_filter_a' in varRunInfo else None,
                                    fileDescriptionForDisplay = "file A")
        bData = _load_variable_data(bFile.file_object, b_variable_technical_name,
                                    dataFilter = varRunInfo['data_filter_function_b'] if 'data_filter_function_b' in varRunInfo else None,
                                    variableToFilterOn = varRunInfo['variable_to_filter_on_b'] if 'variable_to_filter_on_b' in varRunInfo else None,
                                    variableBasedFilter = varRunInfo['variable_based_filter_b'] if 'variable_based_filter_b' in varRunInfo else None,
                                    fileDescriptionForDisplay = "file B")
        
        # pre-check if this data should be plotted and if it should be compared to the longitude and latitude
        include_images_for_this_variable = ((not('shouldIncludeImages' in runInfo)) or (runInfo['shouldIncludeImages']))
        if 'shouldIncludeImages' in varRunInfo :
            include_images_for_this_variable = varRunInfo['shouldIncludeImages']
        do_not_test_with_lon_lat = (not include_images_for_this_variable) or (len(lon_lat_data.keys()) <= 0)
        
        # handle vector data
        isVectorData = ( ('magnitudeName' in varRunInfo)  and ('directionName'  in varRunInfo) and
                         ('magnitudeBName' in varRunInfo) and ('directionBName' in varRunInfo) )
        
        # check if this data can be displayed but
        # don't compare lon/lat sizes if we won't be plotting
        if ( (aData.shape == bData.shape) 
             and 
             ( do_not_test_with_lon_lat
              or
              ((aData.shape == good_shape_from_lon_lat) and (bData.shape == good_shape_from_lon_lat)) ) ) :
            
            # check to see if there is a directory to put information about this variable in,
            # if not then create it
            variableDir = os.path.join(pathsTemp['out'], './' + displayName)
            varRunInfo['variable_dir'] = variableDir
            varRunInfo['variable_report_path_escaped'] = quote(os.path.join(displayName, 'index.html'))
            LOG.debug ("Directory selected for variable information: " + varRunInfo['variable_report_path_escaped'])
            if not (os.path.isdir(variableDir)) :
                LOG.debug("Variable directory (" + variableDir + ") does not exist.")
                LOG.debug("Creating variable directory.")
                os.makedirs(variableDir)
            
            # form the doc and config paths relative to where the variable is
            upwardPath = './'
            for number in range(len(displayName.split('/'))) : # TODO this is not general to windows
                upwardPath = os.path.join(upwardPath, '../')
            varRunInfo['doc_path'] = quote(os.path.join(upwardPath, 'doc.html'))
            if 'config_file_name' in runInfo :
                varRunInfo['config_file_path'] = quote(os.path.join(upwardPath, runInfo['config_file_name']))
            
            # figure out the masks we want, and then do our statistical analysis
            mask_a_to_use = None
            mask_b_to_use = None
            if not do_not_test_with_lon_lat :
                mask_a_to_use = lon_lat_data['a']['inv_mask']
                mask_b_to_use = lon_lat_data['b']['inv_mask']
            variable_stats = statistics.StatisticalAnalysis.withSimpleData(aData, bData,
                                                                           varRunInfo['missing_value'], varRunInfo['missing_value_alt_in_b'],
                                                                           mask_a_to_use, mask_b_to_use,
                                                                           varRunInfo['epsilon'], varRunInfo['epsilon_percent']).dictionary_form()
            
            # add a little additional info to our variable run info before we squirrel it away
            varRunInfo['time'] = datetime.datetime.ctime(datetime.datetime.now())  # todo is this needed?
            didPass, epsilon_failed_fraction, \
                     non_finite_fail_fraction, \
                     r_squared_value = _check_pass_or_fail(varRunInfo, variable_stats, defaultValues)
            varRunInfo['did_pass'] = didPass
            # update the overall pass status
            if didPass is not None :
                didPassAll = didPassAll & didPass
            
            # based on the settings and whether the variable passsed or failed,
            # should we include images for this variable?
            if ('only_plot_on_fail' in varRunInfo) and varRunInfo['only_plot_on_fail'] :
                include_images_for_this_variable = include_images_for_this_variable and (not didPass)
                varRunInfo['shouldIncludeImages'] = include_images_for_this_variable
            
            # to hold the names of any images created
            image_names = {
                            'original': [ ],
                            'compared': [ ]
                            }
            
            # create the images for this variable
            if (include_images_for_this_variable) :
                
                plotFunctionGenerationObjects = [ ]
                
                # if the data is the same size, we can always make our basic statistical comparison plots
                if (aData.shape == bData.shape) :
                    plotFunctionGenerationObjects.append(plotcreate.BasicComparisonPlotsFunctionFactory())
                
                # if the bin and tuple are defined, try to analyze the data as complex
                # multidimentional information requiring careful sampling
                if ('binIndex' in varRunInfo) and ('tupleIndex' in varRunInfo) :
                    plotFunctionGenerationObjects.append(plotcreate.BinTupleAnalysisFunctionFactory())
                    
                else : # if it's not bin/tuple, there are lots of other posibilities
                    
                    # if it's vector data with longitude and latitude, quiver plot it on the Earth
                    if isVectorData and (not do_not_test_with_lon_lat) :
                        plotFunctionGenerationObjects.append(plotcreate.MappedQuiverPlotFunctionFactory())
                    
                    # if the data is one dimensional we can plot it as lines
                    elif   (len(aData.shape) is 1) : 
                        plotFunctionGenerationObjects.append(plotcreate.LinePlotsFunctionFactory())
                    
                    # if the data is 2D we have some options based on the type of data
                    elif (len(aData.shape) is 2) :
                        
                        # if the data is not mapped to a longitude and latitude, just show it as an image
                        if (do_not_test_with_lon_lat) :
                            plotFunctionGenerationObjects.append(plotcreate.IMShowPlotFunctionFactory())
                        
                        # if it's 2D and mapped to the Earth, contour plot it on the earth
                        else :
                            plotFunctionGenerationObjects.append(plotcreate.MappedContourPlotFunctionFactory())
                
                # if there's magnitude and direction data, figure out the u and v, otherwise these will be None
                aUData, aVData = _get_UV_info_from_magnitude_direction_info (aFile.file_object,
                                                                             varRunInfo['magnitudeName'] if ('magnitudeName') in varRunInfo else None,
                                                                             varRunInfo['directionName'] if ('directionName') in varRunInfo else None,
                                                                             lon_lat_data['a']['inv_mask']
                                                                             if ('a' in lon_lat_data) and ('inv_mask' in lon_lat_data['a']) else None)
                bUData, bVData = _get_UV_info_from_magnitude_direction_info (bFile.file_object,
                                                                             varRunInfo['magnitudeBName'] if ('magnitudeBName') in varRunInfo else None,
                                                                             varRunInfo['directionBName'] if ('directionBName') in varRunInfo else None,
                                                                             lon_lat_data['b']['inv_mask']
                                                                             if ('b' in lon_lat_data) and ('inv_mask' in lon_lat_data['b']) else None)
                
                # plot our lon/lat related info
                image_names['original'], image_names['compared'] = \
                    plot.plot_and_save_comparison_figures \
                            (aData, bData,
                             plotFunctionGenerationObjects,
                             varRunInfo['variable_dir'],
                             displayName,
                             varRunInfo['epsilon'],
                             varRunInfo['missing_value'],
                             missingValueAltInB = varRunInfo['missing_value_alt_in_b'] if 'missing_value_alt_in_b' in varRunInfo else None,
                             lonLatDataDict=lon_lat_data,
                             dataRanges     = varRunInfo['display_ranges']      if 'display_ranges'      in varRunInfo else None,
                             dataRangeNames = varRunInfo['display_range_names'] if 'display_range_names' in varRunInfo else None,
                             dataColors     = varRunInfo['display_colors']      if 'display_colors'      in varRunInfo else None,
                             makeSmall=True,
                             doFork=runInfo['doFork'],
                             shouldClearMemoryWithThreads=runInfo['useThreadsToControlMemory'],
                             shouldUseSharedRangeForOriginal=runInfo['useSharedRangeForOriginal'],
                             doPlotSettingsDict = varRunInfo,
                             aUData=aUData, aVData=aVData,
                             bUData=bUData, bVData=bVData,
                             binIndex=      varRunInfo['binIndex']        if 'binIndex'        in varRunInfo else None,
                             tupleIndex=    varRunInfo['tupleIndex']      if 'tupleIndex'      in varRunInfo else None,
                             binName=       varRunInfo['binName']         if 'binName'         in varRunInfo else 'bin',
                             tupleName=     varRunInfo['tupleName']       if 'tupleName'       in varRunInfo else 'tuple',
                             epsilonPercent=varRunInfo['epsilon_percent'] if 'epsilon_percent' in varRunInfo else None,
                             fullDPI=       runInfo['detail_DPI'],
                             thumbDPI=      runInfo['thumb_DPI'],
                             units_a=       varRunInfo['units_a']         if 'units_a'         in varRunInfo else None,
                             units_b=       varRunInfo['units_b']         if 'units_b'         in varRunInfo else None)
                
                LOG.info("\tfinished creating figures for: " + explanationName)
            
            # create the report page for this variable
            if (runInfo['shouldIncludeReport']) :
                
                # hang on to our good % and other info to describe our comparison
                epsilonPassedPercent = (1.0 -  epsilon_failed_fraction) * 100.0
                finitePassedPercent  = (1.0 - non_finite_fail_fraction) * 100.0 
                variableComparisons[displayName] = {'pass_epsilon_percent':   epsilonPassedPercent,
                                                    'finite_similar_percent': finitePassedPercent,
                                                    'r_squared_correlation':  r_squared_value,
                                                    'variable_run_info':      varRunInfo
                                                    }
                
                LOG.info ('\tgenerating report for: ' + explanationName) 
                report.generate_and_save_variable_report(files,
                                                         varRunInfo, runInfo,
                                                         variable_stats,
                                                         spatialInfo,
                                                         image_names,
                                                         varRunInfo['variable_dir'], "index.html")
        
        # if we can't compare the variable, we should tell the user 
        else :
            message = (explanationName + ' ' + 
                     'could not be compared. This may be because the data for this variable does not match in shape ' +
                     'between the two files (file A data shape: ' + str(aData.shape) + '; file B data shape: '
                     + str(bData.shape) + ')')
            if do_not_test_with_lon_lat :
                message = message + '.'
            else :
                message = (message + ' or the data may not match the shape of the selected '
                     + 'longitude ' + str(good_shape_from_lon_lat) + ' and '
                     + 'latitude '  + str(good_shape_from_lon_lat) + ' variables.')
            LOG.warn(message)
        
    # the end of the loop to examine all the variables
    
    # generate our general report pages once we've analyzed all the variables
    if (runInfo['shouldIncludeReport']) :
        
        # get the current time
        runInfo['time'] = datetime.datetime.ctime(datetime.datetime.now())
        
        # make the main summary report
        LOG.info ('generating summary report')
        report.generate_and_save_summary_report(files,
                                                pathsTemp['out'], 'index.html',
                                                runInfo,
                                                variableComparisons, 
                                                spatialInfo,
                                                nameStats)
        
        # make the glossary
        LOG.info ('generating glossary')
        report.generate_and_save_doc_page(statistics.StatisticalAnalysis.doc_strings(), pathsTemp['out'])
    
    returnCode = 0 if didPassAll else 2 # return 2 only if some of the variables failed
    
    # if we are reporting the pass / fail, return an appropriate status code
    if do_pass_fail :
        LOG.debug("Pass/Fail return code: " + str(returnCode))
        return returnCode

def stats_library_call(afn, bfn, var_list=[ ],
                       options_set={ },
                       do_document=False,
                       output_channel=sys.stdout): 
    """
    this method handles the actual work of the stats command line tool and
    can also be used as a library routine, simply pass in an output channel
    and/or use the returned dictionary of statistics for your own form of
    display.
    TODO, should this move to a different file?
    """
    # unpack some options
    epsilon_val  = options_set['epsilon']
    missing_val  = options_set['missing']
    do_pass_fail = options_set['usePassFail']
    
    LOG.debug ("file a: " + afn)
    LOG.debug ("file b: " + bfn)
    
    # open the files
    filesInfo = _open_and_process_files([afn, bfn], 2)
    aFile = filesInfo[afn]['fileObject']
    bFile = filesInfo[bfn]['fileObject']
    
    # information for testing pass/fail if needed
    has_failed = False
    defaultVariablePassFailSettings = {
                                        'epsilon_failure_tolerance': 0.0,
                                        'nonfinite_data_tolerance': 0.0
                                       }
    
    # figure out the variable names and their individual settings
    if len(var_list) <= 0 :
        var_list = ['.*']
    names = _parse_varnames( filesInfo['commonVarNames'], var_list, epsilon_val, missing_val )
    LOG.debug(str(names))
    doc_each  = do_document and len(names)==1
    doc_atend = do_document and len(names)!=1
    
    for name, epsilon, missing in names:
        aData = aFile[name]
        bData = bFile[name]
        if missing is None:
            amiss = aFile.missing_value(name)
            bmiss = bFile.missing_value(name)
        else:
            amiss,bmiss = missing,missing
        LOG.debug('comparing %s with epsilon %s and missing %s,%s' % (name,epsilon,amiss,bmiss))
        print >> output_channel, '-'*32
        print >> output_channel, name
        print >> output_channel, ''
        variable_stats = statistics.StatisticalAnalysis.withSimpleData(aData, bData, amiss, bmiss, epsilon=epsilon)
        # if we're doing pass/fail testing, do that now
        if do_pass_fail :
            didPass, _, _, _ =_check_pass_or_fail(defaultVariablePassFailSettings,
                                                  variable_stats.dictionary_form(),
                                                  glance_analysis_defaults)
            has_failed = has_failed or not(didPass)
        lal = list(variable_stats.dictionary_form().items())
        #lal = list(statistics.summarize(aData, bData, epsilon, (amiss,bmiss)).items()) 
        lal.sort()
        for dictionary_title, dict_data in lal:
            print >> output_channel, '%s' %  dictionary_title
            dict_data
            for each_stat in sorted(list(dict_data)):
                print >> output_channel, '  %s: %s' % (each_stat, dict_data[each_stat])
                if doc_each: print >> output_channel, ('    ' + statistics.StatisticalAnalysis.doc_strings()[each_stat])
            print >> output_channel, '' 
    if doc_atend:
        print >> output_channel, ('\n\n' + statistics.STATISTICS_DOC_STR)
    
    # if we are doing pass/fail, we need to return a status code
    if do_pass_fail :
        status_code = 0
        if has_failed :
            status_code = 3
        LOG.debug("stats is returning status code: " + str(status_code))
        return status_code
    # note: if we aren't doing pass/fail, stats will not return anything

def inspect_stats_library_call (afn, var_list=[ ], options_set={ }, do_document=False, output_channel=sys.stdout): 
    """
    this method handles the actual work of the inspect_stats command line tool and
    can also be used as a library routine, simply pass in an output channel
    and/or use the returned dictionary of statistics for your own form of
    display.
    TODO, should this move to a different file?
    """
    # unpack some options
    missing_val  = options_set['missing']
    
    LOG.debug ("file a: " + afn)
    
    # open the file
    filesInfo = _open_and_process_files([afn], 1)
    aFile = filesInfo[afn]['fileObject']
    
    # figure out the variable names and their individual settings
    if len(var_list) <= 0 :
        var_list = ['.*']
    names = _parse_varnames( filesInfo['commonVarNames'], var_list, epsilon=None, missing=missing_val )
    LOG.debug(str(names))
    doc_each  = do_document and len(names)==1
    doc_atend = do_document and len(names)!=1
    
    for name, epsilon, missing in names:
        aData = aFile[name]
        amiss = missing
        if missing is None:
            amiss = aFile.missing_value(name)
        LOG.debug('analyzing %s with missing data value %s' % (name,amiss))
        print >> output_channel, '-'*32
        print >> output_channel, name
        print >> output_channel, ''
        variable_stats = statistics.StatisticalInspectionAnalysis.withSimpleData(aData, amiss)
        lal = list(variable_stats.dictionary_form().items())
        lal.sort()
        for dictionary_title, dict_data in lal:
            print >> output_channel, '%s' %  dictionary_title
            dict_data
            for each_stat in sorted(list(dict_data)):
                print >> output_channel, '  %s: %s' % (each_stat, dict_data[each_stat])
                if doc_each: print >> output_channel, ('    ' + statistics.StatisticalAnalysis.doc_strings()[each_stat])
            print >> output_channel, '' 
    if doc_atend:
        print >> output_channel, ('\n\n' + statistics.STATISTICS_DOC_STR)

def main():
    import optparse
    usage = """
%prog [options] 
run "%prog help" to list commands
examples:

python -m glance.compare info A.hdf
python -m glance.compare stats A.hdf B.hdf '.*_prof_retr_.*:1e-4' 'nwp_._index:0'
python -m glance.compare plotDiffs A.hdf B.hdf
python -m glance.compare reportGen A.hdf B.hdf
python -m glance.compare gui
python -m glance.compare inspectStats A.hdf

"""
    
    # the following represent options available to the user on the command line:
    
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--test', dest="self_test",
                    action="store_true", default=False, help="run internal unit tests")            
    parser.add_option('-q', '--quiet', dest="quiet",
                    action="store_true", default=False, help="only error output")
    parser.add_option('-v', '--verbose', dest="verbose",
                    action="store_true", default=False, help="enable more informational output")   
    parser.add_option('-w', '--debug', dest="debug",
                    action="store_true", default=False, help="enable debug output")   
    parser.add_option('-e', '--epsilon', dest="epsilon", type='float', default=0.0,
                    help="set default epsilon value for comparison threshold")   
    parser.add_option('-m', '--missing', dest="missing", type='float', default=None,
                    help="set default missing-value")
    #report generation related options
    parser.add_option('-p', '--outputpath', dest="outputpath", type='string', default='./',
                    help="set path to output directory")
    parser.add_option('-o', '--longitude', dest="longitudeVar", type='string',
                    help="set name of longitude variable")
    parser.add_option('-a', '--latitude', dest="latitudeVar", type='string',
                    help="set name of latitude variable")
    parser.add_option('-i', '--imagesonly', dest="imagesOnly", 
                      action="store_true", default=False,
                      help="generate only image files (no html report)")
    parser.add_option('-r', '--reportonly', dest="htmlOnly", 
                      action="store_true", default=False,
                      help="generate only html report files (no images)")
    parser.add_option('-c', '--configfile', dest="configFile", type='string', default=None,
                      help="set optional configuration file")
    parser.add_option('-l', '--llepsilon', dest='lonlatepsilon', type='float', default=0.0,
                      help="set default epsilon for longitude and latitude comparsion")
    parser.add_option('-n', '--version', dest='version',
                      action="store_true", default=False, help="view the glance version")
    parser.add_option('-f', '--fork', dest='doFork',
                      action="store_true", default=False, help="start multiple processes to create images in parallel")
    parser.add_option('-d', '--nolonlat', dest='noLonLatVars',
                      action="store_true", default=False, help="do not try to find or analyze logitude and latitude")
    parser.add_option('-x', '--doPassFail', dest='usePassFail',
                      action="store_true", default=False, help="should the comparison test for pass/fail (currently only affects stats)")
    
    # parse the uers options from the command line
    options, args = parser.parse_args()
    if options.self_test:
        import doctest
        doctest.testmod()
        sys.exit(2)
    
    # set up the logging level based on the options the user selected on the command line
    lvl = logging.WARNING
    if options.debug: lvl = logging.DEBUG
    elif options.verbose: lvl = logging.INFO
    elif options.quiet: lvl = logging.ERROR
    logging.basicConfig(level = lvl)
    
    # display the version
    if options.version :
        print (_get_glance_version_string() + '\n')

    commands = {}
    prior = None
    prior = dict(locals())
    
    """
    The following functions represent available menu selections in glance.
    """
    
    def info(*args):
        """list information about a list of files
        List available variables for comparison.
        """
        for fn in args:
            try :
                lal = list(io.open(fn)())
                lal.sort()
                print fn + ': ' + ('\n  ' + ' '*len(fn)).join(lal)
            except KeyError :
                LOG.warn('Unable to open / process file selection: ' + fn)
    
    def stats(*args):
        """create statistics summary of variables
        Summarize difference statistics between listed variables.
        If no variable names are given, summarize all common variables.
        Variable names can be of the form varname:epsilon:missing to use non-default epsilon or missing value.
        Variable names can be regular expressions, e.g. 'image.*' or '.*prof_retr.*::-999'
        Either epsilon or missing can be empty to stay with default.
        If _FillValue is an attribute of a variable, that will be used to find missing values where no value is given.
        Run with -v to get more detailed information on statistics.
        Examples:
         python -m glance.compare stats hdffile1 hdffile2
         python -m glance.compare stats --epsilon=0.00001 A.hdf B.hdf baseline_cmask_seviri_cloud_mask:0.002:
         python -m glance.compare -w stats --epsilon=0.00001 A.hdf A.hdf imager_prof_retr_abi_total_precipitable_water_low::-999
        """ 
        afn, bfn = args[:2]
        do_doc = (options.verbose or options.debug)
        
        tempOptions = { }
        tempOptions['epsilon']       = options.epsilon
        tempOptions['missing']       = options.missing
        tempOptions['usePassFail']   = options.usePassFail
        # add more if needed for stats
        
        status_result = stats_library_call(_clean_path(afn), _clean_path(bfn),
                                           var_list=args[2:],
                                           options_set=tempOptions,
                                           do_document=do_doc)
        
        if status_result is not None :
            return status_result

    def plotDiffs(*args) :
        """generate a set of images comparing two files
        This option creates a set of graphical comparisons of variables in the two given hdf files.
        The images detailing the differences between variables in the two hdf files will be
        generated and saved to disk. 
        Variables to be compared may be specified after the names of the two input files. If no variables
        are specified, all variables that match the shape of the longitude and latitude will be compared.
        Specified variables that do not exist, do not match the correct data shape, or are the longitude/latitude
        variables will be ignored.
        The user may also use the notation variable_name:epsilon:missing_value to specify the acceptible epsilon
        for comparison and the missing_value which indicates missing data. If one or both of these values is absent
        (in the case of variable_name:epsilon: variable_name::missing_value or just variable_name) the default value
        of 0.0 will be used for epsilon and no missing values will be analyzed. 
        The created images will be stored in the provided path, or if no path is provided, they will be stored in
        the current directory.
        The longitude and latitude variables may be specified with --longitude and --latitude
        If no longitude or latitude are specified the pixel_latitude and pixel_longitude variables will be used.
        Examples:
         python -m glance.compare plotDiffs A.hdf B.hdf variable_name_1:epsilon1: variable_name_2 variable_name_3:epsilon3:missing3 variable_name_4::missing4
         python -m glance.compare --outputpath=/path/where/output/will/be/placed/ plotDiffs A.hdf B.hdf
         python -m glance.compare plotDiffs --longitude=lon_variable_name --latitude=lat_variable_name A.hdf B.hdf variable_name
        """
        # set the options so that a report will not be generated
        options.imagesOnly = True
        
        # make the images
        reportGen(*args)
        
        return

    def reportGen(*args) :
        """generate a report comparing two files
        This option creates a report comparing variables in the two given hdf files.
        An html report and images detailing the differences between variables in the two hdf files will be
        generated and saved to disk. The images will be embedded in the report or visible as separate .png files.
        Variables to be compared may be specified after the names of the two input files. If no variables
        are specified, all variables that match the shape of the longitude and latitude will be compared.
        Specified variables that do not exist, do not match the correct data shape, or are the longitude/latitude
        variables will be ignored.
        The user may also use the notation variable_name:epsilon:missing_value to specify the acceptible epsilon
        for comparison and the missing_value which indicates missing data. If one or both of these values is absent
        (in the case of variable_name:epsilon: variable_name::missing_value or just variable_name) the default value
        of 0.0 will be used for epsilon and no missing values will be analyzed. 
        The html report page(s) and any created images will be stored in the provided path, or if no path is provided,
        they will be stored in the current directory.
        If for some reason you would prefer to generate the report without images, use the --reportonly option. This
        option will generate the html report but omit the images. This may be significantly faster, depending on
        your system, but the differences between the files may be quite a bit more difficult to interpret.
        The longitude and latitude variables may be specified with --longitude and --latitude
        If no longitude or latitude are specified the pixel_latitude and pixel_longitude variables will be used.
        Examples:
         python -m glance.compare reportGen A.hdf B.hdf variable_name_1:epsilon1: variable_name_2 variable_name_3:epsilon3:missing3 variable_name_4::missing4
         python -m glance.compare --outputpath=/path/where/output/will/be/placed/ reportGen A.hdf B.hdf
         python -m glance.compare reportGen --longitude=lon_variable_name --latitude=lat_variable_name A.hdf B.hdf variable_name
         python -m glance.compare reportGen --imagesonly A.hdf B.hdf
        """
        
        tempOptions = { }
        tempOptions['outputpath']    = _clean_path(options.outputpath)
        tempOptions['configFile']    = _clean_path(options.configFile)
        tempOptions['imagesOnly']    = options.imagesOnly
        tempOptions['htmlOnly']      = options.htmlOnly
        tempOptions['doFork']        = options.doFork
        tempOptions['noLonLatVars']  = options.noLonLatVars
        tempOptions['latitudeVar']   = options.latitudeVar
        tempOptions['longitudeVar']  = options.longitudeVar
        tempOptions['lonlatepsilon'] = options.lonlatepsilon
        tempOptions['epsilon']       = options.epsilon
        tempOptions['missing']       = options.missing
        tempOptions['usePassFail']   = options.usePassFail
        
        a_path = _clean_path(args[0])
        b_path = _clean_path(args[1])
        
        return reportGen_library_call(a_path, b_path, args[2:], tempOptions)
    
    def inspectStats(*args):
        """create statistics summary of variables from one file
        Summarize data on variables in a file.
        If no variable names are given, summarize all variables.
        Variable names can be of the form varname::missing to use non-default missing value.
        Variable names can be regular expressions, e.g. 'image.*' or '.*prof_retr.*::-999'
        Missing can be empty to stay with default.
        If _FillValue is an attribute of a variable, that will be used to find missing values where no value is given.
        Run with -v to get more detailed information on the statistics provided.
        Examples:
         python -m glance.compare    inspectStats A.hdf
         python -m glance.compare    inspectStats A.hdf baseline_cmask_seviri_cloud_mask
         python -m glance.compare -w inspectStats A.hdf imager_prof_retr_abi_total_precipitable_water_low::-999
        """ 
        afn = args[0]
        do_doc = (options.verbose or options.debug)
        
        tempOptions = { }
        tempOptions['missing']       = options.missing
        # add more if needed for stats
        
        inspect_stats_library_call(_clean_path(afn), var_list=args[1:], options_set=tempOptions, do_document=do_doc)
    
    def inspectReport(*args) :
        """inspect the contents of a file
        This option creates a report and or images examining variables in a file.
        
        An html report and images detailing the variables in the file will be generated and saved to disk.
        The images will be embedded in the report or visible as separate .png files.
        
        Variables to be compared may be specified after the name of the input file. If no variables
        are specified, all variables that match the shape of the longitude and latitude will be compared.
        Specified variables that do not exist or do not match the lon/lat data shape will be ignored.
        
        The user may also use the notation variable_name::missing_value to specify the missing_value which indicates
        fill data. If this value is absent (in the case of variable_name:: or just variable_name) glance with attempt
        to load the missing value from the file (if none is found, the inspection may fail).
        
        The html report and any created images will be stored in the provided path, or if no path is provided,
        they will be stored in the current directory.
        
        If for some reason you would prefer to generate the report without images, use the --reportonly option. This
        option will generate the html report but omit the images. This may be significantly faster, depending on
        your system, however the information may be more difficult to interpret.
        
        The longitude and latitude variables may be specified with --longitude and --latitude
        If no longitude or latitude are specified the pixel_latitude and pixel_longitude variables will be used by default.
        If no longitude or latitude mappings are desired, the --nolonlat option will disable spatial mapping.
        
        Examples:
         python -m glance.compare inspect_report A.hdf variable_name_1:: variable_name_2 variable_name_3::missing3 variable_name_4::missing4
         python -m glance.compare --outputpath=/path/where/output/will/be/placed/ inspect_report A.hdf 
         python -m glance.compare inspect_report --longitude=lon_variable_name --latitude=lat_variable_name A.hdf variable_name
         python -m glance.compare inspect_report --reportonly A.hdf
        """
        
        tempOptions = { }
        tempOptions['outputpath']    = _clean_path(options.outputpath)
        tempOptions['configFile']    = _clean_path(options.configFile)
        tempOptions['imagesOnly']    = options.imagesOnly
        tempOptions['htmlOnly']      = options.htmlOnly
        tempOptions['doFork']        = False
        tempOptions['noLonLatVars']  = options.noLonLatVars
        tempOptions['latitudeVar']   = options.latitudeVar
        tempOptions['longitudeVar']  = options.longitudeVar
        tempOptions['missing']       = options.missing
        
        # args[0] is the path of the file to be analyzed, an other args should be variable names
        return inspect_library_call(_clean_path(args[0]), args[1:], tempOptions)
    
    def colocateData(*args) :
        """colocate data in two files
        
        This option colocates data in the two given input files and saves it to separate output files.
        Data will be colocated based on its corresponding longitude and latitude. Multiple matches may be
        made between a data point in file A and those in file B if they are within the longitude/latitude epsilon.
        Points from each file that could not be matched and the number of duplicate matches will also be
        recorded in the output file.
        
        The user may also use the notation variable_name::missing_value to specify the missing_value which indicates
        missing data. If no missing value is given, glance will attempt to load a missing value from the input file.
        If there is no missing value defined for that variable in the file, no missing value will be analyzed.
        Missing value data points will not be considered for colocation.
        
        Data which corresponds to longitude or latitude values which fall outside the earth (outside the normally
        accepted valid ranges) will also be considered invalid and will not be considered for colocation.
        
        The longitude and latitude variables may be specified with --longitude and --latitude
        If no longitude or latitude are specified the pixel_latitude and pixel_longitude variables will be used.
        The longitude and latitude epsilon may be specified with --llepsilon
        If no longitude/latitude epsilon is given the value of 0.0 (degrees) will be used
        
        The output data files generated by this option will appear in the selected output directory, or the current
        directory if no out put directory is selected. The output files will be named originalFileName-colocation.nc
        (replacing "originalFileName" with the names of your input files).
        
        Examples:
         python -m glance.compare colocateData A.hdf B.hdf variable_name_1 variable_name_2 variable_name_3::missing3 
         python -m glance.compare colocateData --outputpath=/path/where/output/will/be/placed/ A.nc B.nc
         python -m glance.compare colocateData --longitude=lon_variable_name --latitude=lat_variable_name A.hdf B.hdf variable_name
         python -m glance.compare colocateData --llepsilon=0.0001 A.nc B.hdf
        """
        
        tempOptions = { }
        tempOptions['outputpath']    = _clean_path(options.outputpath)
        tempOptions['configFile']    = _clean_path(options.configFile)
        tempOptions['noLonLatVars']  = options.noLonLatVars
        tempOptions['latitudeVar']   = options.latitudeVar
        tempOptions['longitudeVar']  = options.longitudeVar
        tempOptions['lonlatepsilon'] = options.lonlatepsilon
        tempOptions['epsilon']       = options.epsilon
        tempOptions['missing']       = options.missing
        
        # TODO, remove these eventually
        tempOptions['htmlOnly']      = False
        tempOptions['imagesOnly']    = False
        tempOptions['doFork']        = False
        
        tempOptions['doColocate']    = True
        
        a_path = _clean_path(args[0])
        b_path = _clean_path(args[1])
        
        colocateToFile_library_call(a_path, b_path, args[2:], tempOptions)
    
    # Note: the figure plotting in the GUI is dependant on having selected an interactive renderer in the first "use"
    # statement at the beginning of this module. (It had to be moved into this module to pre-empt other use statempents
    # from imports of other glance modules.)
    def gui (*args) :
        """start the glance graphical user interface
        
        This option launches the graphical user interface for glance. This interface includes only some of the basic
        functionality of glance and may be expanded in the future.
        
        No arguments are required when using this option.
        The various output related arguments (quiet, verbose, debug, etc.) may be used if desired.
        
        Examples:
         python -m glance.compare gui
        """
        
        LOG.debug("Launching Glance GUI")
        temp_controller = gui_control.GlanceGUIController(_get_glance_version_string())
        temp_controller.launch_gui()
    
    # def build(*args):
    #     """build summary
    #     build extended info
    #     """
    #     LOG.info("building database tables")
    #     
    # def grant(*args):
    #     """grant summary
    #     grant extended info
    #     """
    #     LOG.info("granting permissions for tables")
    #     
    # def index(*args):
    #     """index summary
    #     index extended info
    #     """
    #     LOG.info("creating indices for tables")
        
    def help(command=None):
        """print help for a specific command or list of commands
        e.g. help stats
        """
        if command is None: 
            # print first line of docstring
            for cmd in commands:
                ds = commands[cmd].__doc__.split('\n')[0]
                print "%-16s %s" % (cmd,ds)
        else:
            print commands[command].__doc__
            
    # def test():
    #     "run tests"
    #     test1()
    #
    
    # all the local public functions are considered part of glance, collect them up
    commands.update(dict(x for x in locals().items() if x[0] not in prior))    
    
    # if what the user asked for is not one of our existing functions, print the help
    if (not args) or (args[0] not in commands): 
        parser.print_help()
        help()
        return 9
    else:
        # call the function the user named, given the arguments from the command line  
        locals()[args[0]](*args[1:])

    return 0


if __name__=='__main__':
    sys.exit(main())