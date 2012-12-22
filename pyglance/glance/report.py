#!/usr/bin/env python
# encoding: utf-8
"""
PDF/HTML report generation routines

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import sys, logging

from pkg_resources import resource_string, resource_filename #, resource_stream
from mako.template import Template
from mako.lookup   import TemplateLookup
import types as types
import numpy as np
import shutil as shutil

from glance.constants import *

LOG = logging.getLogger(__name__)

# TODO, this should be overridable in the config file when there is one
floatFormat = '%#0.4g'
formattingSettings = {
                      types.FloatType: floatFormat,
                      np.float32: floatFormat,
                      np.float64: floatFormat
                      }

# make and save an html page using a mako template, put all the data you need
# in the template into the kwargs
def _make_and_save_page (fullFilePath, templateFileNameToUse, **kwargs) :
    
    tempFileName = resource_filename(__name__, ".")
    tempLookup = TemplateLookup(directories=[tempFileName])
    
    fileToWrite = open(fullFilePath, 'w')
    tempTemplate = Template(resource_string(__name__, templateFileNameToUse), lookup=tempLookup)

    fileToWrite.write(tempTemplate.render(**kwargs))
    fileToWrite.close()
    
    return

def make_formatted_display_string(displayData, customDisplayFormat=None) :
    '''
    given a piece of data return a display string
    '''
    displayString = ''
    formatStr = customDisplayFormat
    
    # check to see if there is a format string to use
    if type(displayData) in formattingSettings :
        if formatStr is None :
            formatStr = formattingSettings[type(displayData)]
        displayString = formatStr % displayData
    else :
        displayString = str(displayData)
    
    return displayString

def generate_and_save_summary_report(files,
                                     outputPath, reportFileName,
                                     runInfo,
                                     variables,
                                     spatial={},
                                     varNames={}) :
    """
    given two files, and information about them, save a summary of their comparison
    The summary report, in html format will be saved to the given outputPath/outputFile
    
    Variables should be a dictionary keyed on the name of each compared variable and containing the
    % of data values that were judged to be "similar enough" between file A and B (according
    to the epsilon originally inputed for the comparison) and the epsilon used for the comparison
        variables[var_name] = {"pass_epsilon_percent": percent "similar enough" according to epsilon,
                               "variable_run_info": variable Run Info <- as defined for a variable report
                               }
    more keys can be added in the future if more detailed data reporting is desired on the main report page
    
    runInfo should be information about the run in the form
        runInfo = {
                   MACHINE_INFO_KEY:        currentMachine,
                   USER_INFO_KEY:           currentUser,
                   TIME_INFO_KEY:           currentTime,
                   GLANCE_VERSION_INFO_KEY: a version string describing glance,
                   LATITUDE_NAME_KEY:       latitudeName,
                   LONGITUDE_NAME_KEY:      longitudeName,
                   LAT_ALT_NAME_IN_B_KEY:   latitudeNameInB,       # optional, if not defined, B's using the normal latitude
                   LON_ALT_NAME_IN_B_KEY:   longitudeNameInB,     # optional, if not defined, B's using the normal longitude
                   DO_MAKE_IMAGES_KEY:      shouldIncludeImages,      # this key/value is optional, defaults to True
                   SHORT_CIRCUIT_DIFFS_KEY: if the diff related information should be shown
                   }
                   
    files is a dictionary in the form
        files = {A_FILE_TITLE_KEY:
                           {
                            PATH_KEY:          fileAName,
                            #DISPLAY_NAME_KEY:  fileADisplayName, # TODO
                            LAST_MODIFIED_KEY: lastModifiedTimeA,
                            MD5SUM_KEY:        aMD5SUM
                           },
                 B_FILE_TITLE_KEY:
                           {
                            PATH_KEY:          fileBName,
                            #DISPLAY_NAME_KEY:  fileADisplayName, # TODO
                            LAST_MODIFIED_KEY: lastModifiedTimeB,
                            MD5SUM_KEY:        bMD5SUM
                           }
                    }
    
    spatial is a dictionary in the form
        spatial = {
                     A_FILE_TITLE_KEY:
                                {
                                NUMBER_INVALID_PTS_KEY:  number of spatially invalid points only in A (and not corrspondingly in B),
                                PERCENT_INVALID_PTS_KEY: percent of spatially invalid points in A (out of all pts in A)
                                },
                     B_FILE_TITLE_KEY:
                                {
                                NUMBER_INVALID_PTS_KEY:  number of spatially invalid points only in B (and not corrspondingly in A),
                                PERCENT_INVALID_PTS_KEY: percent of spatially invalid points in B (out of all pts in B)
                                },
                     PERCENT_INV_PTS_SHARED_KEY: the percent of points that are spatially invalid in at least one of the two files,
        }
    any of the first level of keys in spatial are optional,
    although it is assumed that if you include an entry for A_FILE_TITLE_KEY or B_FILE_TITLE_KEY it
    will have both of the expected keys in its dictionary
    (ie. both NUMBER_INVALID_PTS_KEY and PERCENT_INVALID_PTS_KEY)
    
    varNames should be information about the variables present in the files, in the form
        varNames = {
                    VAR_NAMES_UNIQUE_TO_A_KEY: uniqueToAVars,
                    VAR_NAMES_UNIQUE_TO_B_KEY: uniqueToBVars,
                    SHARED_VARIABLE_NAMES_KEY: sharedVars
                   }
    all entries in the varNames dictionary are optional.
    """
    
    # pack up all the data needed to build the summary report
    
    
    varNamesToUse = {
                        VAR_NAMES_UNIQUE_TO_A_KEY: [ ],
                        VAR_NAMES_UNIQUE_TO_B_KEY: [ ],
                        SHARED_VARIABLE_NAMES_KEY: [ ]
                    }
    varNamesToUse.update(varNames)
    
    # build the full kwargs with all the info
    kwargs = { 
               RUN_INFO_DICT_KEY:          runInfo,
               FILES_INFO_DICT_KEY:        files,
               SPATIAL_INFO_DICT_KEY:      spatial,
               VARIABLE_NAMES_DICT_KEY:    varNamesToUse,
               VARIABLE_RUN_INFO_DICT_KEY: variables 
               }
              
    _make_and_save_page((outputPath + "/" + reportFileName), 'mainreport.txt', **kwargs)
    
    # copy the pass/fail images, TODO should I move this to a list input in the parameters?
    passFile = resource_filename(__name__, 'pass.gif') # TODO, how can this be done without unzipping the egg?
    failFile = resource_filename(__name__, 'fail.gif') # TODO, how can this be done without unzipping the egg?
    shutil.copy(passFile, outputPath)
    shutil.copy(failFile, outputPath)
    
    # copy the original configuration file, TODO should I move this to a list input in the parameters?
    if (CONFIG_FILE_PATH_KEY in runInfo) :
        originalConfigFile = runInfo[CONFIG_FILE_PATH_KEY]
        shutil.copy(originalConfigFile, outputPath)
    
    return

def generate_and_save_doc_page(definitions, outputPath) :
    """
    generate a page with all the statistics definitions for reference from the reports
    """
    kwargs = {
                DEFINITIONS_INFO_KEY: definitions
             }
    _make_and_save_page(outputPath + "/doc.html", 'doc.txt', **kwargs)
    
    return

def generate_and_save_variable_report(files,
                                      variableRunInfo, # contains variable specific run information
                                      generalRunInfo,  # contains run information not related to the variable
                                      statGroups,
                                      spatial,
                                      imageNames,
                                      outputPath, reportFileName
                                      ) :
    """
    given two files and information about the comparison of one of their variables,
    generate an html report about that variable and store it the outputPath/reportFileName
    provided
    
    statGroups is a dictionary in the form
         statGroups['stat group display name'] = {a dictionary of stats/values to show}
    ie. there may be many different groups of stats that should each be displayed
    
    generalRunInfo is a dictionary in the form
        generalRunInfo = {
                            MACHINE_INFO_KEY:        currentMachine,
                            USER_INFO_KEY:           currentUser,
                            GLANCE_VERSION_INFO_KEY: a version string describing glance,
                            LATITUDE_NAME_KEY:       latitudeName,
                            LONGITUDE_NAME_KEY:      longitudeName,
                            LAT_ALT_NAME_IN_B_KEY:   latitudeNameInB,       # optional, if not defined, B's using the normal latitude
                            LON_ALT_NAME_IN_B_KEY:   longitudeNameInB,     # optional, if not defined, B's using the normal longitude
                            DO_MAKE_IMAGES_KEY:      shouldIncludeImages,
                            SHORT_CIRCUIT_DIFFS_KEY: if the diff related information should be shown
                         }
    
    variableRunInfo is a dictionary in the form
        variableRunInfo = {
                            VARIABLE_TECH_NAME_KEY: variableName,
                            EPSILON_KEY:            epsilon,
                            FILL_VALUE_KEY:         missingDataValue,
                            DISPLAY_NAME_KEY:       displayName
                            DID_VARIABLE_PASS_KEY:      boolean value or None, # optional, boolean means it did or did not pass, None means it was
                                                                           # not qualitatively tested against a set of tolerances
                            'time': currentTime
                          }
                            
    files is a dictionary in the form
        files = {A_FILE_TITLE_KEY:
                        {
                            PATH_KEY:          fileAName,
                            #DISPLAY_NAME_KEY:  fileADisplayName, # TODO
                            LAST_MODIFIED_KEY: lastModifiedTimeA,
                            MD5SUM_KEY:        aMD5SUM
                        },
                 B_FILE_TITLE_KEY:
                        {
                            PATH_KEY:          fileBName,
                            #DISPLAY_NAME_KEY:  fileADisplayName, # TODO
                            LAST_MODIFIED_KEY: lastModifiedTimeB,
                            MD5SUM_KEY:        bMD5SUM
                            }
                    }
                    
    spatial is a dictionary in the form
        spatial = {
                     A_FILE_TITLE_KEY:
                            {
                                NUMBER_INVALID_PTS_KEY:  number of spatially invalid points only in A (and not corrspondingly in B),
                                PERCENT_INVALID_PTS_KEY: percent of spatially invalid points in A (out of all pts in A)
                            },
                     B_FILE_TITLE_KEY:
                            {
                                NUMBER_INVALID_PTS_KEY:  number of spatially invalid points only in B (and not corrspondingly in A),
                                PERCENT_INVALID_PTS_KEY: percent of spatially invalid points in B (out of all pts in B)
                            },
                     PERCENT_INV_PTS_SHARED_KEY: the percent of points that are spatially invalid in at least one of the two files,
        }
    
    imageNames is a dictionary in the form
        imageNames = {
                      ORIGINAL_IMAGES_KEY : [list, of, file-names, of, original, images],
                      COMPARED_IMAGES_KEY : [list, of, file-names, of, compared, images]
                     }
    note: the assumption will be made that smaller versions of these images exist
    in the form small.filename 
    
    any of the first level of keys in spatial are optional,
    although it is assumed that if you include an entry for A_FILE_TITLE_KEY or
    B_FILE_TITLE_KEY it will have both of the expected keys in its dictionary
    (ie. both NUMBER_INVALID_PTS_KEY and PERCENT_INVALID_PTS_KEY)
    
    """
    
    # pack up all the data for a report on a particular variable
    
    # information about the run in general
    runInfo = generalRunInfo.copy()
    runInfo.update(variableRunInfo)
    
    # put all the info together in the kwargs
    kwargs = {
               RUN_INFO_DICT_KEY:        runInfo,
               FILES_INFO_DICT_KEY :     files,
               STATS_INFO_DICT_KEY:      statGroups,
               SPATIAL_INFO_DICT_KEY:    spatial,
               IMAGE_NAME_INFO_DICT_KEY: imageNames
               }
    
    _make_and_save_page((outputPath + "/" + reportFileName), 'variablereport.txt', **kwargs)
    
    return

def generate_and_save_inspect_variable_report(files,
                                              variableRunInfo, # contains variable specific run information
                                              generalRunInfo,  # contains run information not related to the variable
                                              statGroups,
                                              spatial,
                                              imageNames,
                                              outputPath, reportFileName
                                              ) :
    """
    given a file and information about one of the variables in that file,
    generate an html report about that variable and store it the outputPath/reportFileName provided
    
    statGroups is a dictionary in the form
         statGroups['stat group display name'] = {a dictionary of stats/values to show}
    ie. there may be many different groups of stats that should each be displayed
    
    generalRunInfo is a dictionary in the form
        generalRunInfo = {
                            MACHINE_INFO_KEY:        currentMachine,
                            USER_INFO_KEY:           currentUser,
                            GLANCE_VERSION_INFO_KEY: a version string describing glance,
                            LATITUDE_NAME_KEY:       latitudeName,
                            LONGITUDE_NAME_KEY:      longitudeName,
                            DO_MAKE_IMAGES_KEY:      shouldIncludeImages,
                         }
    
    variableRunInfo is a dictionary in the form
        variableRunInfo = { VARIABLE_TECH_NAME_KEY: variableName,
                            FILL_VALUE_KEY:         missingDataValue,
                            DISPLAY_NAME_KEY:       displayName,
                            TIME_INFO_KEY:          currentTime
                            }
                            
    files is a dictionary in the form
        files = {
                    A_FILE_TITLE_KEY:
                        {
                            PATH_KEY:          fileAName,
                            #DISPLAY_NAME_KEY:  fileADisplayName, # TODO
                            LAST_MODIFIED_KEY: lastModifiedTimeA,
                            MD5SUM_KEY:        aMD5SUM
                        },
                }
                    
    spatial is a dictionary in the form
        spatial = {
                    NUMBER_INVALID_PTS_KEY:  number of spatially invalid points only in A (and not corrspondingly in B),
                    PERCENT_INVALID_PTS_KEY: percent of spatially invalid points in A (out of all pts in A)
                  }
    
    imageNames is a dictionary in the form
        imageNames = {
                        ORIGINAL_IMAGES_KEY : [list, of, file-names, of, original, images],
                        COMPARED_IMAGES_KEY : [list, of, file-names, of, analyzed, images]
                     }
    note: the assumption will be made that smaller versions of these images exist
    in the form small.filename 
    
    """
    
    # pack up all the data for a report on a particular variable
    
    # information about the run in general
    runInfo = generalRunInfo.copy()
    runInfo.update(variableRunInfo)
    
    # put all the info together in the kwargs
    kwargs = { 
               RUN_INFO_DICT_KEY:        runInfo,
               FILES_INFO_DICT_KEY :     files,
               STATS_INFO_DICT_KEY:      statGroups,
               SPATIAL_INFO_DICT_KEY:    spatial,
               IMAGE_NAME_INFO_DICT_KEY: imageNames
               }
    
    _make_and_save_page((outputPath + "/" + reportFileName), 'inspectvariablereport.txt', **kwargs)
    
    return

def generate_and_save_inspection_summary_report(files,
                                                outputPath, reportFileName,
                                                runInfo,
                                                variables,
                                                spatial={},
                                                varNames={}) :
    """
    given two files, and information about them, save a summary of their comparison
    The summary report, in html format will be saved to the given outputPath/outputFile
    
    Variables should be a dictionary keyed on the name of each compared variable and containing the
    % of data values that were judged to be "similar enough" between file A and B (according
    to the epsilon originally inputed for the comparison) and the epsilon used for the comparison
        variables[var_name] = {"pass_epsilon_percent": percent "similar enough" according to epsilon,
                               "variable_run_info": variable Run Info <- as defined for a variable report
                               }
    more keys can be added in the future if more detailed data reporting is desired on the main report page
    
    runInfo should be information about the run in the form
        runInfo =   {
                        MACHINE_INFO_KEY:        currentMachine,
                        USER_INFO_KEY:           currentUser,
                        TIME_INFO_KEY:           currentTime,
                        GLANCE_VERSION_INFO_KEY: a version string describing glance,
                        LATITUDE_NAME_KEY:       latitudeName,
                        LONGITUDE_NAME_KEY:      longitudeName,
                        DO_MAKE_IMAGES_KEY:      shouldIncludeImages,      # this key/value is optional, defaults to True
                    }
                   
    files is a dictionary in the form
        files =     {
                        A_FILE_TITLE_KEY:
                            {
                                PATH_KEY:          fileAName,
                                #DISPLAY_NAME_KEY:  fileADisplayName, # TODO
                                LAST_MODIFIED_KEY: lastModifiedTimeA,
                                MD5SUM_KEY:        aMD5SUM
                            },
                    }
    
    spatial is a dictionary in the form
        spatial = {
                    NUMBER_INVALID_PTS_KEY:  number of spatially invalid points only in A (and not corrspondingly in B),
                    PERCENT_INVALID_PTS_KEY: percent of spatially invalid points in A (out of all pts in A)
                  }
    
    varNames should be information about the variables present in the file, in the form
    more information may be added to this dictionary in the future
        varNames = {
                        POSSIBLE_NAMES_KEY: listOfAllVarsInFile,
                    }
    """
    
    # pack up all the data needed to build the summary report
    
    # TODO, more automated defaults
    varNameDefaults = { POSSIBLE_NAMES_KEY: [ ],}
    varNamesToUse = varNameDefaults
    varNamesToUse.update(varNames)
    
    # build the full kwargs with all the info
    kwargs = {
                RUN_INFO_DICT_KEY: runInfo,
                FILES_INFO_DICT_KEY: files,
                SPATIAL_INFO_DICT_KEY: spatial,
                VARIABLE_NAMES_DICT_KEY: varNamesToUse,
                VARIABLE_RUN_INFO_DICT_KEY: variables 
               }
              
    _make_and_save_page((outputPath + "/" + reportFileName), 'inspectmainreport.txt', **kwargs)
    
    # copy the original configuration file, TODO should I move this to a list input in the parameters?
    if (CONFIG_FILE_PATH_KEY in runInfo) :
        originalConfigFile = runInfo[CONFIG_FILE_PATH_KEY]
        shutil.copy(originalConfigFile, outputPath)
    
    return

if __name__=='__main__':
    sys.exit(main())