#!/usr/bin/env python
# encoding: utf-8
"""
Plotting routines for difference values using matplotlib

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys

# these first two lines must stay before the pylab import
import matplotlib
#matplotlib.use('Qt4Agg') # use the Anti-Grain Geometry rendering engine

from pylab import *

import matplotlib.pyplot as plt

from PIL import Image

import os, sys, logging
import numpy as np

import glance.graphics as maps
import glance.delta    as delta
import glance.figures  as figures
import glance.data     as dataobj
import glance.plotcreatefns as plotfns
from glance.constants import *

LOG = logging.getLogger(__name__)

# a constant for the larger size dpi
fullSizeDPI = 150 # 200
# a constant for the thumbnail size dpi
thumbSizeDPI = 50

def _handle_fig_creation_task(child_figure_function, log_message,
                              outputPath, fullFigName,
                              shouldMakeSmall, doFork,
                              fullDPI=fullSizeDPI, thumbDPI=thumbSizeDPI) :
    """
    fork to do something.
    the parent will return the child pid
    the child will do it's work and then exit
    """
    
    pid = 0
    if (doFork) :
        # do the fork
        pid = os.fork()
    
    # figure out if we're the parent or child
    isParent = not (pid is 0)
    if (isParent) :
        return pid
    else :
        plt.ioff()
        figure = child_figure_function() 
        LOG.info(log_message)
        if figure is not None :
            figure.savefig(os.path.join(outputPath, fullFigName), dpi=fullDPI)
            if (shouldMakeSmall) :

                tempImage = Image.open(os.path.join(outputPath, fullFigName))
                scaleFactor = float(thumbDPI) / float(fullDPI)
                originalSize = tempImage.size
                newSize = (int(originalSize[0] * scaleFactor), int(originalSize[1] * scaleFactor))
                tempImage = tempImage.resize(newSize, Image.ANTIALIAS)
                tempImage.save(os.path.join(outputPath, 'small.' + fullFigName))

            # get rid of the figure
            plt.close(figure)
            del(figure)
        else:
            LOG.warn("Unable to create plot.")
    
    # if we've reached this point and we did fork,
    # then we're the child process and we should stop now
    if (doFork) :
        sys.exit(0) # the child is done now
    
    # if we didn't fork, return the 0 pid to indicate that
    return pid

def _log_spawn_and_wait_if_needed (imageDescription, childPids, 
                                   taskFunction, taskOutputPath, taskFigName,
                                   doMakeThumb=True, doFork=False, shouldClearMemoryWithThreads=False,
                                   fullDPI=fullSizeDPI, thumbDPI=thumbSizeDPI) :
    """
    create a figure generation task, spawning a process as needed
    save the childPid to the list of pids if the process will remain outstanding after this method ends
    the name of the figure that was generated will be added to the image list
    """
    LOG.info("creating image of "+ imageDescription)
    
    # start the actual task
    pid = _handle_fig_creation_task(taskFunction,
                                    "saving image of " + imageDescription,
                                    taskOutputPath, taskFigName,
                                    doMakeThumb, doFork or shouldClearMemoryWithThreads,
                                    fullDPI=fullDPI, thumbDPI=thumbDPI)
    
    # wait based on the state of the pid we received and why we would have forked
    childPid = None
    if not (pid is 0) :
        if doFork :
            childPids.append(pid)
            LOG.debug ("Started child process (pid: " + str(pid) + ") to create image of " + imageDescription)
        else :
            os.waitpid(pid, 0)
    
    return

def plot_and_save_spacial_mismatch(longitudeObject, latitudeObject, spacialMismatchMask,
                                  fileNameDiscriminator, title, fileBaseName, outputPath, makeSmall=False,
                                  fullDPI=fullSizeDPI, thumbDPI=thumbSizeDPI, units=None) :
    """
    given information on spatially placed mismatch points in A and B, plot only those points in a very obvious way
    on top of a background plot of a's data shown in grayscale, save this plot to the output path given
    if makeSmall is passed as true a smaller version of the image will also be saved
    """
    spaciallyInvalidMask = longitudeObject.masks.ignore_mask | latitudeObject.masks.ignore_mask
    
    # get the bounding axis and make a basemap
    boundingAxes = plotfns.get_visible_axes(longitudeObject.data, latitudeObject.data, ~spaciallyInvalidMask)
    LOG.debug("Visible axes for lon/lat mismatch figure  are: " + str(boundingAxes))
    baseMapInstance, boundingAxes = maps.create_basemap(longitudeObject.data, latitudeObject.data,
                                                        boundingAxes, plotfns.select_projection(boundingAxes))
    
    # make the figure
    LOG.info("Creating spatial mismatch image")
    spatialMismatchFig = figures.create_mapped_figure(None, latitudeObject.data, longitudeObject.data, baseMapInstance,
                                                      boundingAxes, title, invalidMask=spaciallyInvalidMask,
                                                      tagData=spacialMismatchMask, units=units)
    # save the figure
    LOG.info("Saving spatial mismatch image")
    spatialMismatchFig.savefig(outputPath + "/" + fileBaseName + "." + fileNameDiscriminator + ".png", dpi=fullDPI) 
    
    # we may also save a smaller versions of the figure
    if (makeSmall) :
        
        spatialMismatchFig.savefig(outputPath + "/" + fileBaseName + "." + fileNameDiscriminator + ".small.png", dpi=thumbDPI)
    
    # get rid of the figure
    spatialMismatchFig.clf()
    plt.close(spatialMismatchFig)
    del(spatialMismatchFig)
    
    return

def plot_and_save_comparison_figures (aData, bData,
                                     plottingFunctionFactoryObjects,
                                     outputPath,
                                     variableDisplayName,
                                     epsilon,
                                     missingValue,
                                     missingValueAltInB=None,
                                     lonLatDataDict=None,
                                     dataRanges=None,
                                     dataRangeNames=None,
                                     dataColors=None,
                                     makeSmall=False,
                                     shortCircuitComparisons=False,
                                     doFork=False,
                                     shouldClearMemoryWithThreads=False,
                                     shouldUseSharedRangeForOriginal=False,
                                     doPlotSettingsDict={ },
                                     aUData=None, aVData=None,
                                     bUData=None, bVData=None,
                                     binIndex=None, tupleIndex=None,
                                     binName='bin', tupleName='tuple',
                                     epsilonPercent=None,
                                     fullDPI=None, thumbDPI=None,
                                     units_a=None, units_b=None,
                                     useBData=True,
                                     histRange=None) :
    """
    Plot images for a set of figures based on the data sets and settings
    passed in. The images will be saved to disk according to the settings.
    
    aData -             the A file data
    bData -             the B file Data
    lonLatDataDict -    a dictionary of longitude and latitude info in the form
                        (or it may be empty if there is no lon/lat info or the
                         available lon/lat info should not be used for this data
                         set)
    
        lonLatDataDict = {
                          A_FILE_KEY =
                                {
                                 LON_KEY:          longitudeDataForFileA,
                                 LAT_KEY:          latitudeDataForFileA,
                                 INVALID_MASK_KEY: invalidMaskForFileA
                                },
                          B_FILE_KEY =
                                {
                                 LON_KEY:          longitudeDataForFileB,
                                 LAT_KEY:          latitudeDataForFileB,
                                 INVALID_MASK_KEY: invalidMaskForFileB
                                },
                          COMMON_KEY =
                                {
                                 LON_KEY:          longitudeDataCommonToBothFiles,
                                 LAT_KEY:          latitudeDataCommonToBothFiles,
                                 INVALID_MASK_KEY: invalidMaskCommonToBothFiles
                                }
                          }
    
    required parameters:
    
    outputPath -        the path where the output images will be placed
    variableDisplayName - a descriptive name that will be shown on the
                          images to label the variable being analyzed
    epsilon -           the epsilon that should be used for the data comparison
    missingValue -      the missing value that should be used for the data comparison
    plottingFunctionFactoryObjects - a list of objects that will generate the plotting
                                     functions; each object in the list must be an 
                                     instance of a child of PlottingFunctionFactory
    
    optional parameters:
    
    missingValueAltInB - an alternative missing value in the B file; if None is given
                         then the missingValue will be used for the B file
    makeSmall -         should small "thumbnail" images be made for each large image?
    shortCircuitComparisons - should the comparison plots be disabled?
    doFork -            should the process fork to create new processes to create each
                        image? **
    shouldClearMemoryWithThreads - should the process use fork to control long term
                                   memory bloating? **
    shouldUseSharedRangeForOriginal - should the original images share an all-inclusive
                                      data range?
    doPlotSettingsDict - a dictionary containting settings to turn off individual plots
                         if an entry for a plot is not pressent or set to True,
                         then the plot will be created
    useBData -           should the b data be expected or not? when this is false,
                         b data will not be used and no lon/lat data for b will be
                         expected either
    histRange -          the range that should be used for the histogram, or None
    
    ** May fail due to a known bug on MacOSX systems.
    """
    
    # lists to hold information on the images we make
    original_images = [ ]
    compared_images = [ ]
    
    # figure out what missing values we should be using
    if missingValueAltInB is None :
        missingValueAltInB = missingValue
    
    # figure out if we have spatially invalid masks to consider
    spaciallyInvalidMaskA = None
    spaciallyInvalidMaskB = None
    if (lonLatDataDict is not None) and (len(lonLatDataDict.keys()) > 0):
        if useBData :
            spaciallyInvalidMaskA = lonLatDataDict[A_FILE_KEY][INVALID_MASK_KEY]
            spaciallyInvalidMaskB = lonLatDataDict[B_FILE_KEY][INVALID_MASK_KEY]
        else :
            spaciallyInvalidMaskA = lonLatDataDict[INVALID_MASK_KEY]
    
    # compare the two data sets to get our difference data and mismatch info
    aDataObject = dataobj.DataObject(aData, fillValue=missingValue, ignoreMask=spaciallyInvalidMaskA)
    bDataObject = None
    diffInfo    = None
    if useBData :
        bDataObject = dataobj.DataObject(bData, fillValue=missingValueAltInB, ignoreMask=spaciallyInvalidMaskB)
        diffInfo = dataobj.DiffInfoObject(aDataObject, bDataObject, epsilonValue=epsilon, epsilonPercent=epsilonPercent)
    else :
        aDataObject.self_analysis() # if we aren't going to do a diff, make sure basic analysis is done
    
    # from this point on, we're forking to create child processes to parallelize image and report generation
    # and to help control memory leak creeping that matplotlib causes
    isParent = True 
    childPids = [ ]
    
    plottingFunctions = { }
    
    for factoryObject in plottingFunctionFactoryObjects :
        
        # generate our plotting functions
        moreFunctions = factoryObject.create_plotting_functions (
                                       # the most basic data set needed
                                       aDataObject, bDataObject,
                                       variableDisplayName,
                                       epsilon,
                                       doPlotSettingsDict,
                                       
                                       # where the names of the created figures will be stored
                                       original_images, compared_images,
                                       
                                       # parameters that are only needed for geolocated data
                                       lonLatDataDict=lonLatDataDict,
                                       
                                       # only used if we are plotting a contour
                                       dataRanges=dataRanges, dataRangeNames=dataRangeNames,
                                       dataColors=dataColors,
                                       shouldUseSharedRangeForOriginal=shouldUseSharedRangeForOriginal,
                                       
                                       # a comparison of the data if the data comparison info is needed
                                       differences=diffInfo,
                                       
                                       # only used for plotting quiver data
                                       aUData=aUData, aVData=aVData,
                                       bUData=bUData, bVData=bVData,
                                       
                                       # only used for line plots 
                                       binIndex=binIndex, tupleIndex=tupleIndex,
                                       binName=binName, tupleName=tupleName,
                                       epsilonPercent=epsilonPercent,
                                       
                                       # used for display in several types of plots
                                       units_a=units_a, units_b=units_b,
                                       
                                       # range for a histogram
                                       histRange=histRange
                                       )
        plottingFunctions.update(moreFunctions)
    
    LOG.debug ('plotting function information: ' + str(plottingFunctions))
    
    # for each function in the list, run it to create the figure
    for figDesc in sorted(list(plottingFunctions.keys())) :
        
        LOG.info("plotting " + figDesc + " figure for " + variableDisplayName)
        # returnDictionary['descriptive name'] = (function, title, file_name, list_this_figure_should_go_into)
        figFunction, figLongDesc, figFileName, outputInfoList = plottingFunctions[figDesc]
        
        # only plot the compared images if we aren't short circuiting them
        if (outputInfoList is not compared_images) or (not shortCircuitComparisons) :
            _log_spawn_and_wait_if_needed(figLongDesc, childPids, figFunction, outputPath, figFileName,
                                          makeSmall, doFork, shouldClearMemoryWithThreads, fullDPI=fullDPI, thumbDPI=thumbDPI)
            # if we made an attempt to make the file, hang onto the name
            outputInfoList.append(figFileName)
    
    # now we need to wait for all of our child processes to terminate before returning
    if (isParent) : # just in case
        if len(childPids) > 0 :
            LOG.info ("waiting for completion of " + variableDisplayName + " images...")
        for pid in childPids:
            os.waitpid(pid, 0)
        LOG.info("... creation and saving of images for " + variableDisplayName + " completed")
    
    return original_images, compared_images

if __name__=='__main__':
    import doctest
    doctest.testmod()
