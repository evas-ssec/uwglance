#!/usr/bin/env python
# encoding: utf-8
"""
PDF/HTML report generation routines

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging

from pkg_resources import resource_string
from mako.template import Template

LOG = logging.getLogger(__name__)

# make and save an html page using a mako template, put all the data you need
# in the template into the kwargs
def _make_and_save_page (fullFilePath, templateFileNameToUse, **kwargs) :
    
    fileToWrite = open(fullFilePath, 'w')
    tempTemplate = Template(resource_string(__name__, templateFileNameToUse))

    fileToWrite.write(tempTemplate.render(**kwargs))
    fileToWrite.close()
    
    return

def generate_and_save_summary_report(fileAName, fileBName,
                                     aMD5SUM, bMD5SUM,
                                     outputPath, reportFileName,
                                     longitudeName, latitudeName,
                                     variables,
                                     lastModifiedTimeA, lastModifiedTimeB,
                                     currentTime, currentUser,
                                     currentMachine,
                                     percentSpaciallyInvalidPts=None,
                                     uniqueToAVars=[], uniqueToBVars=[], sharedVars=[]) :
    """
    given two files, and information about them, save a summary of their comparison
    The summary report, in html format will be saved to the given outputPath/outputFile
    Variables should be a dictionary containing the name of each compared variable and the
    % of data values that were judged to be "similar enough" between file A and B (according
    to the epsilon originally inputed for the comparison)
    """
    kwargs = {'fileAName': fileAName,
              'fileBName': fileBName,
              'aMD5SUM': aMD5SUM,
              'bMD5SUM': bMD5SUM,
              'longitudeName': longitudeName,
              'latitudeName': latitudeName,
              'percentSpaciallyInvalidPts': percentSpaciallyInvalidPts,
              'uniqueToAVars': uniqueToAVars,
              'uniqueToBVars': uniqueToBVars,
              'sharedVars': sharedVars,
              'comparedVariables': variables,
              'lastModifiedTimeA': lastModifiedTimeA,
              'lastModifiedTimeB': lastModifiedTimeB,
              'currentTime': currentTime,
              'currentUser': currentUser,
              'currentMachine': currentMachine
              }
    _make_and_save_page((outputPath + "/" + reportFileName), './mainreport.txt', **kwargs)
    
    return

def generate_and_save_doc_page(definitions, outputPath) :
    """
    generate a page with all the statistics definitions for reference from the reports
    """
    kwargs = {'definitions': definitions
              }
    _make_and_save_page(outputPath + "/doc.html", './doc.txt', **kwargs)
    
    return

def generate_and_save_variable_report(fileAName, fileBName,
                                      aMD5SUM, bMD5SUM,
                                      variableName,
                                      outputPath, reportFileName,
                                      epsilon, missingDataValue,
                                      statsData, shouldIncludeImages,
                                      lastModifiedTimeA, lastModifiedTimeB,
                                      currentTime, currentUser,
                                      currentMachine) :
    """
    given two files and information about the comparison of one of their variables,
    generate an html report about that variable and store it the outputPath/reportFileName
    provided
    """
    kwargs = {'fileAName': fileAName,
              'fileBName': fileBName,
              'aMD5SUM' : aMD5SUM,
              'bMD5SUM' : bMD5SUM,
              'variableName': variableName,
              'statsData': statsData,
              'epsilon': epsilon,
              'missingDataValue': missingDataValue,
              'shouldIncludeImages': shouldIncludeImages,
              'lastModifiedTimeA': lastModifiedTimeA,
              'lastModifiedTimeB': lastModifiedTimeB,
              'currentTime': currentTime,
              'currentUser': currentUser,
              'currentMachine': currentMachine
              }
    _make_and_save_page((outputPath + "/" + reportFileName), './variablereport.txt', **kwargs)
    
    return

if __name__=='__main__':
    sys.exit(main())