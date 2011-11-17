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
    
    tempLookup = TemplateLookup(directories=[resource_filename(__name__, "stuff")])
    
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
        runInfo = {'machine': currentMachine,
                   'user': currentUser,
                   'time': currentTime,
                   'version': a version string describing glance,
                   'latitude': latitudeName,
                   'longitude': longitudeName,
                   'latitude_alt_name_in_b': latitudeNameInB,       # optional, if not defined, B's using the normal latitude
                   'longitude_alt_name_in_b': longitudeNameInB,     # optional, if not defined, B's using the normal longitude
                   'shouldIncludeImages': shouldIncludeImages,      # this key/value is optional, defaults to True
                   'short_circuit_diffs': if the diff related information should be shown
                   }
                   
    files is a dictionary in the form
        files = {'file A': {'path': fileAName,
                            #'displayName': fileADisplayName, # TODO
                            'lastModifiedTime': lastModifiedTimeA,
                            'md5sum': aMD5SUM
                            },
                 'file B': {'path': fileBName,
                            #'displayName': fileADisplayName, # TODO
                            'lastModifiedTime': lastModifiedTimeB,
                            'md5sum': bMD5SUM
                            }
                    }
    
    spatial is a dictionary in the form
        spatial = {
                     'file A': {'numInvPts': number of spatially invalid points only in A (and not corrspondingly in B),
                                'perInvPts': percent of spatially invalid points in A (out of all pts in A)
                                },
                     'file B': {'numInvPts': number of spatially invalid points only in B (and not corrspondingly in A),
                                'perInvPts': percent of spatially invalid points in B (out of all pts in B)
                                },
                     'perInvPtsInBoth': the percent of points that are spatially invalid in at least one of the two files,
                     'num_lon_lat_not_equal_points': number of points that do not match in the sets of lon/lat
                                # if the 'num_lon_lat_not_equal_points' key is defined it means there are mismatching
                                # longitude/latitude pairs in the data that was compared!
        }
    any of the first level of keys in spatial are optional,
    although it is assumed that if you include an entry for 'file A' or 'file B' it
    will have both of the expected keys in its dictionary (ie. both numInvPts and perInvPts)
    
    varNames should be information about the variables present in the files, in the form
        varNames = {'uniqueToAVars': uniqueToAVars,
                    'uniqueToBVars': uniqueToBVars,
                    'sharedVars': sharedVars
                    }
    """
    
    # pack up all the data needed to build the summary report
    
    # TODO, more automated defaults
    varNameDefaults = { 'uniqueToAVars': [],
                        'uniqueToBVars': [],
                        'sharedVars': []}
    varNamesToUse = varNameDefaults
    varNamesToUse.update(varNames)
    
    # build the full kwargs with all the info
    kwargs = { 'runInfo': runInfo,
               'files': files,
               'spatial': spatial,
               'varNames': varNamesToUse,
               'variables': variables 
               }
              
    _make_and_save_page((outputPath + "/" + reportFileName), 'stuff/mainreport.txt', **kwargs)
    
    # copy the pass/fail images, TODO should I move this to a list input in the parameters?
    passFile = resource_filename(__name__, 'stuff/pass.gif') # TODO, how can this be done without unzipping the egg?
    failFile = resource_filename(__name__, 'stuff/fail.gif') # TODO, how can this be done without unzipping the egg?
    shutil.copy(passFile, outputPath)
    shutil.copy(failFile, outputPath)
    
    # copy the original configuration file, TODO should I move this to a list input in the parameters?
    if ('config_file_path' in runInfo) :
        originalConfigFile = runInfo['config_file_path']
        shutil.copy(originalConfigFile, outputPath)
    
    return

def generate_and_save_doc_page(definitions, outputPath) :
    """
    generate a page with all the statistics definitions for reference from the reports
    """
    kwargs = {'definitions': definitions
              }
    _make_and_save_page(outputPath + "/doc.html", 'stuff/doc.txt', **kwargs)
    
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
        generalRunInfo = {  'machine': currentMachine,
                            'user': currentUser,
                            'version': a version string describing glance,
                            'latitude': latitudeName,
                            'longitude': longitudeName,
                            'latitude_alt_name_in_b': latitudeNameInB,       # optional, if not defined, B's using the normal latitude
                            'longitude_alt_name_in_b': longitudeNameInB,     # optional, if not defined, B's using the normal longitude
                            'shouldIncludeImages': shouldIncludeImages,
                            'short_circuit_diffs': if the diff related information should be shown
                            }
    
    variableRunInfo is a dictionary in the form
        variableRunInfo = { 'variable_name': variableName,
                            'epsilon': epsilon,
                            'missing_value': missingDataValue,
                            'display_name': displayName
                            'did_pass': boolean value or None, # optional, boolean means it did or did not pass, None means it was
                                                               # not qualitatively tested against a set of tolerances
                            'time': currentTime
                            }
                            
    files is a dictionary in the form
        files = {'file A': {'path': fileAName,
                            #'displayName': fileADisplayName, # TODO
                            'lastModifiedTime': lastModifiedTimeA,
                            'md5sum': aMD5SUM
                            },
                 'file B': {'path': fileBName,
                            #'displayName': fileADisplayName, # TODO
                            'lastModifiedTime': lastModifiedTimeB,
                            'md5sum': bMD5SUM
                            }
                    }
                    
    spatial is a dictionary in the form
        spatial = {
                     'file A': {'numInvPts': number of spatially invalid points only in A (and not corrspondingly in B),
                                'perInvPts': percent of spatially invalid points in A (out of all pts in A)
                                },
                     'file B': {'numInvPts': number of spatially invalid points only in B (and not corrspondingly in A),
                                'perInvPts': percent of spatially invalid points in B (out of all pts in B)
                                },
                     'perInvPtsInBoth': the percent of points that are spatially invalid in at least one of the two files,
                     'num_lon_lat_not_equal_points': number of points that do not match in the sets of lon/lat
                                # if the 'num_lon_lat_not_equal_points' key is defined it means there are mismatching
                                # longitude/latitude pairs in the data that was compared!
        }
    
    imageNames is a dictionary in the form
        imageNames = {
                      'original' : [list, of, file-names, of, original, images],
                      'compared' : [list, of, file-names, of, compared, images]
        }
    note: the assumption will be made that smaller versions of these images exist
    in the form small.filename 
    
    any of the first level of keys in spatial are optional,
    although it is assumed that if you include an entry for 'file A' or 'file B' it
    will have both of the expected keys in its dictionary (ie. both numInvPts and perInvPts)
    
    """
    
    # pack up all the data for a report on a particular variable
    
    # information about the run in general
    runInfo = generalRunInfo.copy()
    runInfo.update(variableRunInfo)
    
    # put all the info together in the kwargs
    kwargs = { 'runInfo': runInfo,
               'files' : files,
               'statGroups': statGroups,
               'spatial': spatial,
               'imageNames': imageNames
               }
    
    _make_and_save_page((outputPath + "/" + reportFileName), 'stuff/variablereport.txt', **kwargs)
    
    return

if __name__=='__main__':
    sys.exit(main())