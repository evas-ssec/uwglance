#!/usr/bin/env python
# encoding: utf-8
"""
This is a place for the general util functions that multiple parts of
glance need to use.

Created by evas Dec 2012.
Copyright (c) 2012 University of Wisconsin SSEC. All rights reserved.
"""

import os, logging
import pkg_resources
import numpy
from subprocess import check_call

LOG = logging.getLogger(__name__)

def get_glance_version_string() :
    version_num = pkg_resources.require('uwglance')[0].version
    
    return "glance, version " + str(version_num) 

def get_run_identification_info( ) :
    """
    get info about what user/machine/version of glance is being used
    """
    
    # get info on who's doing the run and where
    machine = os.uname()[1]               # the name of the machine running the report
    user    = os.getenv("LOGNAME")        # the name of the user running the report
    version = get_glance_version_string() # the version number of glance
    
    return machine, user, version

def clean_path(string_path) :
    """
    Return a clean form of the path without any '.', '..', or '~'
    """
    clean_path = None
    if string_path is not None :
        clean_path = os.path.abspath(os.path.expanduser(string_path))
    
    return clean_path

def setup_dir_if_needed(dirPath, descriptionName) :
    """
    create the directory if that is needed, if not don't
    """
    if not (os.path.isdir(dirPath)) :
        LOG.info("Specified " + descriptionName + " directory (" + dirPath + ") does not exist.")
        LOG.info("Creating " + descriptionName + " directory.")
        os.makedirs(dirPath)

def _uri_needs_rsync(uri_to_check) :
    """
    check if the uri requires an rsync in order to access the data
    this will return some false positives if you phrase local uri's with the machine name
    for ex. you are on the machine "lotus" and you use the path "rsync:://lotus/data/"
    """
    return not os.path.exists(uri_to_check)

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
        check_call(cmd)
    
    return newPaths

def get_percentage_from_mask(dataMask) :
    """
    given a mask that marks the elements we want the percentage of as True (and is the size of our original data),
    figure out what percentage of the whole they are
    """
    numMarkedDataPts = numpy.sum(dataMask)
    totalDataPts     = dataMask.size
    
    # avoid dividing by 0
    if totalDataPts is 0 :
        return 0.0, 0
    
    percentage = 100.0 * float(numMarkedDataPts) / float(totalDataPts)
    
    return percentage, numMarkedDataPts

if __name__=='__main__':
    pass
