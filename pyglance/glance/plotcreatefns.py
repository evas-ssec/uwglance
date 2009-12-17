#!/usr/bin/env python
# encoding: utf-8
"""
Plotting figure generation functions.

Created by evas Dec 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

# these first two lines must stay before the pylab import
import matplotlib
matplotlib.use('Agg') # use the Anti-Grain Geometry rendering engine

from pylab import *

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.ticker import FormatStrFormatter

from PIL import Image

import os, sys, logging
import numpy as np
from numpy import ma 

import glance.graphics as maps
import glance.delta    as delta
import glance.report   as report
import glance.filters  as filters # TODO, should this be moved?
import glance.figures  as figures

LOG = logging.getLogger(__name__)

# a constant for the larger size dpi
fullSizeDPI = 150 # 200
# a constant for the thumbnail size dpi
thumbSizeDPI = 50

# make a custom medium grayscale color map for putting our bad data on top of
mediumGrayColorMapData = {
    'red'   : ((0.0, 1.00, 1.00),
               (0.5, 0.60, 0.60),
               (1.0, 0.20, 0.20)),
    'green' : ((0.0, 1.00, 1.00),
               (0.5, 0.60, 0.60),
               (1.0, 0.20, 0.20)),
    'blue'  : ((0.0, 1.00, 1.00),
               (0.5, 0.60, 0.60),
               (1.0, 0.20, 0.20))
}
mediumGrayColorMap = colors.LinearSegmentedColormap('mediumGrayColorMap', mediumGrayColorMapData, 256)

# ********************* Section of utility functions ***********************

# a method to stop people from calling my fake abstract methods
def _abstract( ) :
    raise NotImplementedError('Method must be implemented in subclass.')

# figure out the bounding axes for the display given a set of
# longitude and latitude and possible a mask of invalid values
# that we should ignore in them
def get_visible_axes(longitudeData, latitudeData, toIgnoreMask) :
    
    # calculate the bounding range for the display
    # this is in the form [longitude min, longitude max, latitude min, latitude max]
    visibleAxes = [delta.min_with_mask(longitudeData, toIgnoreMask),
                   delta.max_with_mask(longitudeData, toIgnoreMask),
                   delta.min_with_mask(latitudeData,  toIgnoreMask),
                   delta.max_with_mask(latitudeData,  toIgnoreMask)]
    
    return visibleAxes

def select_projection(boundingAxes) :
    """
    chose a map projection based on the bounding axes that will be shown
    """
    
    # TODO at the moment the default (lcc) is cutting off the field of view,
    # I think the same problem will occur with all of the conic projections, because
    # they all allow you to specify the field of view either as corners or as a width/height in
    # meters, but they do not take the distortion that a conic projection causes into account.
    # This means that either the corners or the bottom curve of the data area will be clipped off
    # the edge of the screen. There is also some sort of persistent bug that causes the basemap
    # to ignore the requested corners for some data sets and shift west and north, cutting off
    # a pretty considerable amount of data. I've tried various tactics to control the field of
    # view and can't find any way to get the basemap to show an acceptable area programatically
    # that will match arbitrary data sets.
    # For the moment, I am setting this to use a cylindrical projection rather than a conic.
    # At some point in the future this should be revisited so that glance will be able to handle
    # a wider range of projections.
    projToUse = 'cyl'
    
    # how big is the field of view?
    longitudeRange  = abs(boundingAxes[1] - boundingAxes[0])
    latitudeRange   = abs(boundingAxes[3] - boundingAxes[2])
    # chose the projection based on the range we have to cover
    if (longitudeRange > 180) :
        projToUse = 'cyl' # use a equidistant cylindrical projection to show the whole world
    elif (longitudeRange > 100) or (latitudeRange > 70) :
        projToUse = 'ortho' # use an orthographic projection to show about half the globe
    
    return projToUse

def _make_axis_and_basemap(lonLatDataDict, goodInAMask, goodInBMask, shouldUseSharedRangeForOriginal=False, variableDisplayName=None) :
    """
    Determine the appropriate axes for the given data (in longitude and latitude) and create the appropriate basemap and shared
    range information.
    
    fullAxis, baseMapInstance, sharedRange = _make_axis_and_basemap(lonLatDataDict, goodInAMask, goodInBMask,
                                                                                   shouldUseSharedRangeForOriginal=False)
    """
    
    nameMessage = ''
    if variableDisplayName is not None:
        nameMessage = " (" + variableDisplayName + ")"
    
    # figure out the bounding axis
    aAxis = get_visible_axes(lonLatDataDict['a']['lon'], lonLatDataDict['a']['lat'], ~goodInAMask)
    bAxis = get_visible_axes(lonLatDataDict['b']['lon'], lonLatDataDict['b']['lat'], ~goodInBMask)
    fullAxis = [min(aAxis[0], bAxis[0]), max(aAxis[1], bAxis[1]),
                min(aAxis[2], bAxis[2]), max(aAxis[3], bAxis[3])]
    
    LOG.debug("Visible axes for file A variable data" + nameMessage + " are: " + str(aAxis))
    LOG.debug("Visible axes for file B variable data" + nameMessage + " are: " + str(bAxis))
    LOG.debug("Visible axes shared for both file's variable data" + nameMessage + " are: " + str(fullAxis))
    
    if (fullAxis[0] is None) or (fullAxis[1] is None) or (fullAxis[2] is None) or (fullAxis[3] is None) :
        LOG.warn("Unable to display figures for variable" + nameMessage + " because of inability to identify" +
                 " usable bounding longitude and latitude range on the earth. Bounding range that was identified:" + str(fullAxis))
        return # TODO, the figures need to be disabled from the report and possibly a warning on the report? throw exception instead?
    
    # create our basemap
    LOG.info('\t\tloading base map data')
    baseMapInstance, fullAxis = maps.create_basemap(lonLatDataDict['common']['lon'], lonLatDataDict['common']['lat'],
                                                    fullAxis, select_projection(fullAxis))
    
    # figure out the shared range for A and B's data, by default don't share a range
    sharedRange = None
    if (shouldUseSharedRangeForOriginal) :
        sharedRange = figures._make_range(aData, ~goodInAMask, 50, offset_to_range=figures.offsetToRange,
                                   data_b=bData, invalid_b_mask=~goodInBMask)
    
    return fullAxis, baseMapInstance, sharedRange

# a method to make a list of index numbers for reordering a multi-dimensional array
def _make_new_index_list(numberOfIndexes, firstIndexNumber=0, lastIndexNumber=None) :
    """
    the first and last index numbers represent the dimensions you want to be first and last (respectively)
    when the list is reordered; any other indexes will retain their relative ordering
    
    newIndexList = _make_new_index_list(numIndexes, binIndex, tupleIndex)
    
    TODO, move this to delta?
    """
    
    if lastIndexNumber is None:
        lastIndexNumber = numberOfIndexes - 1
    
    newIndexList = range(numberOfIndexes)
    maxSpecial   = max(firstIndexNumber, lastIndexNumber)
    minSpecial   = min(firstIndexNumber, lastIndexNumber)
    del(newIndexList[maxSpecial])
    del(newIndexList[minSpecial])
    newIndexList = [firstIndexNumber] + newIndexList + [lastIndexNumber]
    
    return newIndexList

# ********************* Section of public classes ***********************

"""
this class is intended to be a parent class for our plotting function
creation classes
it contains a "fake" "abstract" method
"""
class PlottingFunctionFactory :
    
    def create_plotting_functions (
                                   
                                   self,
                                   
                                   # the most basic data set needed
                                   aData, bData,
                                   variableDisplayName,
                                   epsilon,
                                   goodInAMask, goodInBMask,
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # parameters that are only used if the data can be compared
                                   # point by point
                                   absDiffData=None, rawDiffData=None,
                                   goodInBothMask=None,
                                   troubleMask=None, outsideEpsilonMask=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None
                                   
                                   ) : _abstract

"""
This class creates the most basic of comparison plots based on two similarly
sized data sets. (Plots created include histogram and scatter plots.)
"""
class BasicComparisonPlotsFunctionFactory (PlottingFunctionFactory) :
    def create_plotting_functions (
                                   self,
                                   
                                   # the most basic data set needed
                                   aData, bData,
                                   variableDisplayName,
                                   epsilon,
                                   goodInAMask, goodInBMask,
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # parameters that are only used if the data can be compared
                                   # point by point
                                   absDiffData=None, rawDiffData=None,
                                   goodInBothMask=None,
                                   troubleMask=None, outsideEpsilonMask=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None
                                   ) :
        
        functionsToReturn = { }
        
        # make the histogram plot
        if ('do_plot_histogram' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_histogram']) :
            
            assert(goodInBothMask.shape == rawDiffData.shape)
            
            # setup the data bins for the histogram
            numBinsToUse = 50
            valuesForHist = rawDiffData[goodInBothMask]
            functionsToReturn['histogram'] = ((lambda : figures.create_histogram(valuesForHist, numBinsToUse,
                                                                         ("Difference in\n" + variableDisplayName),
                                                                         ('Value of (Data File B - Data File A) at a Data Point'),
                                                                         ('Number of Data Points with a Given Difference'),
                                                                         True)),
                                              "histogram of the amount of difference in " + variableDisplayName,
                                              "Hist.png", compared_fig_list)
        # make the scatter plot
        if ('do_plot_scatter' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_scatter']) :
            
            assert(aData.shape    == bData.shape)
            assert(bData.shape    == goodInBothMask.shape)
            assert(goodInBothMask.shape == outsideEpsilonMask.shape)
            
            functionsToReturn['scatter']   = ((lambda : figures.create_scatter_plot(aData[goodInBothMask], bData[goodInBothMask],
                                                                            "Value in File A vs Value in File B",
                                                                            "File A Value", "File B Value",
                                                                            outsideEpsilonMask[goodInBothMask],
                                                                            epsilon)),
                                              "scatter plot of file a values vs file b values for " + variableDisplayName,
                                              "Scatter.png", compared_fig_list)
        
        return functionsToReturn

"""
This class creates contour plots mapped onto a region of the earth.
"""
class MappedContourPlotFunctionFactory (PlottingFunctionFactory) :
    def create_plotting_functions (
                                   self,
                                   
                                   # the most basic data set needed
                                   aData, bData,
                                   variableDisplayName,
                                   epsilon,
                                   goodInAMask, goodInBMask,
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # parameters that are only used if the data can be compared
                                   # point by point
                                   absDiffData=None, rawDiffData=None,
                                   goodInBothMask=None,
                                   troubleMask=None, outsideEpsilonMask=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None
                                   ) :
        
        # the default for plotting geolocated data
        mappedPlottingFunction = figures.create_mapped_figure
        
        functionsToReturn = { }
        
        assert(lonLatDataDict is not None)
        assert(goodInAMask    is not None)
        assert(goodInBMask    is not None)
        
        fullAxis, baseMapInstance, sharedRange = _make_axis_and_basemap(lonLatDataDict,
                                                                        goodInAMask, goodInBMask,
                                                                        shouldUseSharedRangeForOriginal,
                                                                        variableDisplayName)
        
        # make the plotting functions
        
        # make the original data plots
        if ('do_plot_originals' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_originals']) :
            
            assert('a'   in lonLatDataDict)
            assert('lat' in lonLatDataDict['a'])
            assert('lon' in lonLatDataDict['a'])
            assert(lonLatDataDict['a']['lat'].shape == lonLatDataDict['a']['lon'].shape)
            
            functionsToReturn['originalA'] = ((lambda : mappedPlottingFunction(aData,
                                                                               lonLatDataDict['a']['lat'], 
                                                                               lonLatDataDict['a']['lon'],
                                                                               baseMapInstance, fullAxis,
                                                                               (variableDisplayName + "\nin File A"),
                                                                               invalidMask=(~goodInAMask),
                                                                               dataRanges=dataRanges or sharedRange,
                                                                               dataRangeNames=dataRangeNames,
                                                                               dataRangeColors=dataColors)),
                                              variableDisplayName + " in file a",
                                              "A.png",  original_fig_list)
            
            assert('b'   in lonLatDataDict)
            assert('lat' in lonLatDataDict['b'])
            assert('lon' in lonLatDataDict['b'])
            assert(lonLatDataDict['b']['lat'].shape == lonLatDataDict['b']['lon'].shape)
            
            functionsToReturn['originalB'] = ((lambda : mappedPlottingFunction(bData, 
                                                                               lonLatDataDict['b']['lat'], 
                                                                               lonLatDataDict['b']['lon'],
                                                                               baseMapInstance, fullAxis,
                                                                               (variableDisplayName + "\nin File B"),
                                                                               invalidMask=(~goodInBMask),
                                                                               dataRanges=dataRanges or sharedRange,
                                                                               dataRangeNames=dataRangeNames,
                                                                               dataRangeColors=dataColors)),
                                              variableDisplayName + " in file b",
                                              "B.png",  original_fig_list)
        
        # make the absolute value difference plot
        if ('do_plot_abs_diff' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_abs_diff']) :
            
            assert(absDiffData.shape == goodInBothMask.shape)
            assert('common' in lonLatDataDict)
            assert('lat'    in lonLatDataDict['common'])
            assert('lon'    in lonLatDataDict['common'])
            assert(lonLatDataDict['common']['lat'].shape == lonLatDataDict['common']['lon'].shape)
            assert(lonLatDataDict['common']['lon'].shape == absDiffData.shape)
            
            functionsToReturn['diffAbs']   = ((lambda : mappedPlottingFunction(absDiffData,
                                                                               lonLatDataDict['common']['lat'], 
                                                                               lonLatDataDict['common']['lon'],
                                                                               baseMapInstance, fullAxis,
                                                                               ("Absolute value of difference in\n"
                                                                                + variableDisplayName),
                                                                               invalidMask=(~goodInBothMask))),
                                              "absolute value of difference in " + variableDisplayName,
                                              "AbsDiff.png", compared_fig_list)
        # make the subtractive difference plot
        if ('do_plot_sub_diff' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_sub_diff']) :
            
            assert(rawDiffData.shape == goodInBothMask.shape)
            assert('common' in lonLatDataDict)
            assert('lat'    in lonLatDataDict['common'])
            assert('lon'    in lonLatDataDict['common'])
            assert(lonLatDataDict['common']['lat'].shape == lonLatDataDict['common']['lon'].shape)
            assert(lonLatDataDict['common']['lon'].shape == rawDiffData.shape)
            
            functionsToReturn['diffSub']   = ((lambda : mappedPlottingFunction(rawDiffData, 
                                                                               lonLatDataDict['common']['lat'], 
                                                                               lonLatDataDict['common']['lon'],
                                                                               baseMapInstance, fullAxis,
                                                                               ("Value of (Data File B - Data File A) for\n"
                                                                                + variableDisplayName),
                                                                               invalidMask=(~goodInBothMask))),
                                              "the difference in " + variableDisplayName,
                                              "Diff.png",    compared_fig_list)
        # make the trouble data plot
        if ('do_plot_trouble' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_trouble']) :
            
            assert(aData.shape == bData.shape)
            assert(goodInAMask.shape == goodInBMask.shape)
            assert('common' in lonLatDataDict)
            assert('lat'    in lonLatDataDict['common'])
            assert('lon'    in lonLatDataDict['common'])
            assert(lonLatDataDict['common']['lat'].shape == lonLatDataDict['common']['lon'].shape)
            assert(lonLatDataDict['common']['lon'].shape == aData.shape)
            
            # this is not an optimal solution, but we need to have at least somewhat valid data at any mismatched points so
            # that our plot won't be totally destroyed by missing or non-finite data from B
            bDataCopy = bData[:]
            tempMask = goodInAMask & (~goodInBMask) 
            bDataCopy[tempMask] = aData[tempMask]
            functionsToReturn['trouble']   = ((lambda : mappedPlottingFunction(bDataCopy, 
                                                                               lonLatDataDict['common']['lat'], 
                                                                               lonLatDataDict['common']['lon'],
                                                                               baseMapInstance, fullAxis,
                                                                               ("Areas of trouble data in\n" + variableDisplayName),
                                                                               invalidMask=(~(goodInAMask | goodInBMask)),
                                                                               colorMap=mediumGrayColorMap, tagData=troubleMask,
                                                                               dataRanges=dataRanges,
                                                                               dataRangeNames=dataRangeNames)), # TODO, does this need modification?
                                              "trouble data in " + variableDisplayName,
                                              "Trouble.png", compared_fig_list)
        
        return functionsToReturn

"""
This class creates quiver plots mapped onto a region of the earth.
Note: the plotting function requires u and v data for both data sets, but
the size of the two data sets is not required to be the same. If the size
of the two data sets is the same, additional comparison plots will be
created.
TODO, this class has not been fully tested
"""
class MappedQuiverPlotFunctionFactory (PlottingFunctionFactory) :
    def create_plotting_functions (
                                   self,
                                   
                                   # the most basic data set needed
                                   aData, bData,
                                   variableDisplayName,
                                   epsilon,
                                   goodInAMask, goodInBMask,
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # parameters that are only used if the data can be compared
                                   # point by point
                                   absDiffData=None, rawDiffData=None,
                                   goodInBothMask=None,
                                   troubleMask=None, outsideEpsilonMask=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None
                                   ) :
        
        # the default for plotting geolocated data
        mappedPlottingFunction = figures.create_quiver_mapped_figure
        
        functionsToReturn = { }
        
        assert(aUData is not None)
        assert(aVData is not None)
        assert(bUData is not None)
        assert(bVData is not None)
        
        assert(lonLatDataDict is not None)
        assert(goodInAMask    is not None)
        assert(goodInBMask    is not None)
        
        fullAxis, baseMapInstance, _ = _make_axis_and_basemap(lonLatDataDict, goodInAMask, goodInBMask, variableDisplayName=variableDisplayName)
        
        # make the plotting functions
        
        # make the original data plots
        if ('do_plot_originals' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_originals']) :
            
            assert('a'   in lonLatDataDict)
            assert('lat' in lonLatDataDict['a'])
            assert('lon' in lonLatDataDict['a'])
            assert(lonLatDataDict['a']['lat'].shape == lonLatDataDict['a']['lon'].shape)
            
            functionsToReturn['originalA'] = ((lambda : mappedPlottingFunction(aData,
                                                                               lonLatDataDict['a']['lat'], 
                                                                               lonLatDataDict['a']['lon'],
                                                                               baseMapInstance, fullAxis,
                                                                               (variableDisplayName + "\nin File A"),
                                                                               invalidMask=(~goodInAMask),
                                                                               uData=aUData, vData=aVData)),
                                              variableDisplayName + " in file a",
                                              "A.png",  original_fig_list)
            
            assert('b'   in lonLatDataDict)
            assert('lat' in lonLatDataDict['b'])
            assert('lon' in lonLatDataDict['b'])
            assert(lonLatDataDict['b']['lat'].shape == lonLatDataDict['b']['lon'].shape)
            
            functionsToReturn['originalB'] = ((lambda : mappedPlottingFunction(bData, 
                                                                               lonLatDataDict['b']['lat'], 
                                                                               lonLatDataDict['b']['lon'],
                                                                               baseMapInstance, fullAxis,
                                                                               (variableDisplayName + "\nin File B"),
                                                                               invalidMask=(~ goodInBMask),
                                                                               uData=bUData, vData=bVData)),
                                              variableDisplayName + " in file b",
                                              "B.png",  original_fig_list)
            
            # TODO, any additional figures of the original data?
        
        if (aUData.shape == bUData.shape) :
            LOG.info("creating comparison quiver plots for " + variableDisplayName)
            
            # TODO, is there any complication in taking the diff of vectors this way?
            diffUData = aUData - bUData
            diffVData = aVData - bVData
            
            # make the absolute value difference plot
            if ('do_plot_abs_diff' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_abs_diff']) :
                
                assert(absDiffData.shape == goodInBothMask.shape)
                assert('common' in lonLatDataDict)
                assert('lat'    in lonLatDataDict['common'])
                assert('lon'    in lonLatDataDict['common'])
                assert(lonLatDataDict['common']['lat'].shape == lonLatDataDict['common']['lon'].shape)
                assert(lonLatDataDict['common']['lon'].shape == absDiffData.shape)
                
                functionsToReturn['diffAbs']   = ((lambda : mappedPlottingFunction(absDiffData,
                                                                                   lonLatDataDict['common']['lat'], 
                                                                                   lonLatDataDict['common']['lon'],
                                                                                   baseMapInstance, fullAxis,
                                                                                   ("Absolute value of difference in\n"
                                                                                    + variableDisplayName),
                                                                                   invalidMask=(~ goodInBothMask),
                                                                                   uData=diffUData, vData=diffVData)),
                                                  "absolute value of difference in " + variableDisplayName,
                                                  "AbsDiff.png", compared_fig_list)
            # make the subtractive difference plot
            if ('do_plot_sub_diff' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_sub_diff']) :
                
                assert(rawDiffData.shape == goodInBothMask.shape)
                assert('common' in lonLatDataDict)
                assert('lat'    in lonLatDataDict['common'])
                assert('lon'    in lonLatDataDict['common'])
                assert(lonLatDataDict['common']['lat'].shape == lonLatDataDict['common']['lon'].shape)
                assert(lonLatDataDict['common']['lon'].shape == rawDiffData.shape)
                
                functionsToReturn['diffSub']   = ((lambda : mappedPlottingFunction(rawDiffData, 
                                                                                   lonLatDataDict['common']['lat'], 
                                                                                   lonLatDataDict['common']['lon'],
                                                                                   baseMapInstance, fullAxis,
                                                                                   ("Value of (Data File B - Data File A) for\n"
                                                                                    + variableDisplayName),
                                                                                   invalidMask=(~ goodInBothMask),
                                                                                   uData=diffUData, vData=diffVData)),
                                                  "the difference in " + variableDisplayName,
                                                  "Diff.png",    compared_fig_list)
            # make the trouble data plot
            if ('do_plot_trouble' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_trouble']) :
                
                assert(aData.shape == bData.shape)
                assert(goodInAMask.shape == goodInBMask.shape)
                assert('common' in lonLatDataDict)
                assert('lat'    in lonLatDataDict['common'])
                assert('lon'    in lonLatDataDict['common'])
                assert(lonLatDataDict['common']['lat'].shape == lonLatDataDict['common']['lon'].shape)
                assert(lonLatDataDict['common']['lon'].shape == aData.shape)
                
                # this is not an optimal solution, but we need to have at least somewhat valid data at any mismatched points so
                # that our plot won't be totally destroyed by missing or non-finite data from B
                bDataCopy = bData[:]
                tempMask = goodInAMask & (~goodInBMask) 
                bDataCopy[tempMask] = aData[tempMask]
                functionsToReturn['trouble']   = ((lambda : mappedPlottingFunction(bDataCopy, 
                                                                                   lonLatDataDict['common']['lat'], 
                                                                                   lonLatDataDict['common']['lon'],
                                                                                   baseMapInstance, fullAxis,
                                                                                   ("Areas of trouble data in\n" + variableDisplayName),
                                                                                   invalidMask=(~(goodInAMask | goodInBMask)),
                                                                                   colorMap=mediumGrayColorMap, tagData=troubleMask,
                                                                                   dataRanges=dataRanges,
                                                                                   dataRangeNames=dataRangeNames,
                                                                                   # TODO, does this need modification?
                                                                                   uData=bUData, vData=bVData)), 
                                                  "trouble data in " + variableDisplayName,
                                                  "Trouble.png", compared_fig_list)
        
        return functionsToReturn


"""
This class creates simple line plots based on simple one dimentional data.
"""
class LinePlotsFunctionFactory (PlottingFunctionFactory) :
    def create_plotting_functions (
                                   self,
                                   
                                   # the most basic data set needed
                                   aData, bData,
                                   variableDisplayName,
                                   epsilon,
                                   goodInAMask, goodInBMask,
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # parameters that are only used if the data can be compared
                                   # point by point
                                   absDiffData=None, rawDiffData=None,
                                   goodInBothMask=None,
                                   troubleMask=None, outsideEpsilonMask=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None
                                   ) :
        """
        This method generates line plotting functions for one dimensional data
        and returns them in a dictionary of tupples, where each tupple is in the form:
        
            returnDictionary['descriptive name'] = (function, title, file_name, list_this_figure_should_go_into)
        
        The file name is only the name of the file, not the full path.
        """
        assert(aData is not None)
        assert(bData is not None)
        # TODO, assert about more of the data
        
        if len(aData.shape) > 1 :
            assert(binIndex   is not None)
            assert(tupleIndex is not None)
            assert(binIndex   is not tupleIndex)
            
            numIndexes   = len(aData.shape)
            assert(binIndex   < numIndexes)
            assert(tupleIndex < numIndexes)
        
        # TODO, the rest of this function is not totally refactored but should be functional for now
        
        # TODO, temporary
        aList       = [ ]
        bList       = [ ]
        singleAList = [ ]
        singleBList = [ ]
        absDiffList = [ ]
        subDiffList = [ ]
        troubleSingleList = [ ]
        troubleFullList   = [ ]
        if len(aData.shape) > 1 :
            
            newIndexList = _make_new_index_list(numIndexes, binIndex, tupleIndex)
            aData              =              aData.transpose(newIndexList)
            bData              =              bData.transpose(newIndexList)
            goodInAMask        =        goodInAMask.transpose(newIndexList)
            goodInBMask        =        goodInBMask.transpose(newIndexList)
            absDiffData        =        absDiffData.transpose(newIndexList)
            rawDiffData        =        rawDiffData.transpose(newIndexList)
            goodInBothMask     =     goodInBothMask.transpose(newIndexList)
            troubleMask        =        troubleMask.transpose(newIndexList)
            outsideEpsilonMask = outsideEpsilonMask.transpose(newIndexList)
            
            for firstIndexValue in range(aData.shape[0]) :
                # get the a and b data with it's masks
                aSlice           = filters.collapse_to_index(aData[firstIndexValue], (numIndexes - 1), missing_value=-999.0, ignore_below_exclusive=-990.0)
                bSlice           = filters.collapse_to_index(bData[firstIndexValue], (numIndexes - 1), missing_value=-999.0, ignore_below_exclusive=-990.0)
                goodInAMaskSlice = filters.collapse_to_index(goodInAMask[firstIndexValue], (numIndexes - 1), collapsing_function=(lambda data: any(data)))
                goodInBMaskSlice = filters.collapse_to_index(goodInBMask[firstIndexValue], (numIndexes - 1), collapsing_function=(lambda data: any(data)))
                
                # add to the plain a and b lists
                aList.append((aSlice, ~goodInAMaskSlice, 'r', 'A data, ' + str(firstIndexValue), None))
                bList.append((bSlice, ~goodInBMaskSlice, 'b', 'B data, ' + str(firstIndexValue), None))
                
                # deal with the trouble data and it's list
                troubleMaskSlice = filters.collapse_to_index(troubleMask[firstIndexValue], (numIndexes - 1), collapsing_function=(lambda data: any(data)))
                troubleFullList.append((aSlice, ~goodInAMaskSlice, 'r', 'A data, ' + str(firstIndexValue), troubleMaskSlice))
                troubleFullList.append((bSlice, ~goodInBMaskSlice, 'b', 'B data, ' + str(firstIndexValue), troubleMaskSlice))
                
                # get the diff data with it's masks
                absDiffDataSlice = filters.collapse_to_index(absDiffData[firstIndexValue], (numIndexes - 1), missing_value=-999.0, ignore_below_exclusive=-990.0)
                rawDiffDataSlice = filters.collapse_to_index(rawDiffData[firstIndexValue], (numIndexes - 1), missing_value=-999.0, ignore_below_exclusive=-990.0)
                goodMaskSlice    = filters.collapse_to_index(goodInBothMask[firstIndexValue], (numIndexes - 1), collapsing_function=(lambda data: any(data)))
                
                # add to the diff data lists
                absDiffList.append((absDiffDataSlice, ~goodMaskSlice, '', 'abs. diff. data, ' + str(firstIndexValue), None))
                subDiffList.append((rawDiffDataSlice, ~goodMaskSlice, '', 'sub. diff. data, ' + str(firstIndexValue), None))
                
            # also collapse the masks
            troubleMask        = filters.collapse_to_index(troubleMask,        numIndexes, collapsing_function=(lambda data: any(data)))
            #outsideEpsilonMask = filters.collapse_to_index(outsideEpsilonMask, numIndexes, collapsing_function=(lambda data: any(data)))
            #goodInBothMask     = filters.collapse_to_index(goodInBothMask,     numIndexes, collapsing_function=(lambda data: any(data)))
            
            # make some averages
            aDataAvg           = filters.collapse_to_index(aData,       numIndexes, missing_value=-999.0, ignore_below_exclusive=-990.0)
            bDataAvg           = filters.collapse_to_index(bData,       numIndexes, missing_value=-999.0, ignore_below_exclusive=-990.0)
            goodInAMaskAvg     = filters.collapse_to_index(goodInAMask, numIndexes, collapsing_function=(lambda data: any(data)))
            goodInBMaskAvg     = filters.collapse_to_index(goodInBMask, numIndexes, collapsing_function=(lambda data: any(data)))
            
            # make our averaged sets
            singleAList       = [(aDataAvg, ~goodInAMaskAvg, 'r', 'A data', None)]
            singleBList       = [(bDataAvg, ~goodInBMaskAvg, 'b', 'B data', None)]
            troubleSingleList = [(aDataAvg, ~goodInAMaskAvg, 'r', 'A data', troubleMask),
                                 (bDataAvg, ~goodInBMaskAvg, 'b', 'B data', troubleMask)]
        else :
            aList = [(aData, ~goodInAMask, 'r', 'A data', None)]
            bList = [(bData, ~goodInBMask, 'b', 'B data', None)]
            absDiffList = [(absDiffData, ~goodInBothMask, '', 'abs. diff. data', None)]
            subDiffList = [(rawDiffData, ~goodInBothMask, '', 'sub. diff. data', None)]
            
            singleAList = aList
            singleBList = bList
            
            troubleFullList   = [(aData, ~goodInAMask, 'r', 'A data', troubleMask),
                                 (bData, ~goodInBMask, 'b', 'B data', troubleMask)]
            troubleSingleList = troubleFullList
        
        # make a list of both the A and B values together
        allList       = aList       + bList
        allSingleList = singleAList + singleBList
        
        functionsToReturn = { }
        
        # make the original data plots
        if ('do_plot_originals' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_originals']) :
            
            
            functionsToReturn['original'] = ((lambda: figures.create_line_plot_figure(allList,
                                                                               variableDisplayName + "\nin Both Files")),
                                             variableDisplayName + " in both files",
                                             "AB.png", original_fig_list)
            functionsToReturn['originalA'] = ((lambda: figures.create_line_plot_figure(aList,
                                                                                variableDisplayName + "\nin File A")),
                                              variableDisplayName + " in file a",
                                              "A.png",  original_fig_list)
            functionsToReturn['originalB'] = ((lambda: figures.create_line_plot_figure(bList,
                                                                                variableDisplayName + "\nin File B")),
                                              variableDisplayName + " in file b",
                                              "B.png",  original_fig_list)
            
            if len(allSingleList) < len(allList) :
                functionsToReturn['originalS'] = ((lambda: figures.create_line_plot_figure(allSingleList,
                                                                                    variableDisplayName + "\nAvg. in Both Files")),
                                                  variableDisplayName + " avg. in both files",
                                                  "ABsingle.png", original_fig_list)
                functionsToReturn['originalAS'] = ((lambda: figures.create_line_plot_figure(singleAList,
                                                                                     variableDisplayName + "\nAvg. in File A")),
                                                   variableDisplayName + " avg. in file a",
                                                   "Asingle.png",  original_fig_list)
                functionsToReturn['originalBS'] = ((lambda: figures.create_line_plot_figure(singleBList,
                                                                                     variableDisplayName + "\nAvg. in File B")),
                                                   variableDisplayName + " avg. in file b",
                                                   "Bsingle.png",  original_fig_list)
        
        # make the absolute value difference plot
        if ('do_plot_abs_diff' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_abs_diff']) :
            functionsToReturn['diffAbs']   = ((lambda: figures.create_line_plot_figure(absDiffList,
                                                                               "Absolute value of difference in\n" + variableDisplayName)),
                                              "absolute value of difference in " + variableDisplayName,
                                              "AbsDiff.png", compared_fig_list)
        # make the subtractive difference plot
        if ('do_plot_sub_diff' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_sub_diff']) :
            functionsToReturn['diffSub']   = ((lambda: figures.create_line_plot_figure(subDiffList,
                                                                               "Value of (Data File B - Data File A) for\n" + variableDisplayName)),
                                              "the difference in " + variableDisplayName,
                                              "Diff.png",    compared_fig_list)
        
        # make the trouble data plot
        if ('do_plot_trouble' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_trouble']) :
            functionsToReturn['trouble']   = ((lambda: figures.create_line_plot_figure(troubleFullList,
                                                                               "Areas of trouble data in\n" + variableDisplayName)),
                                              "trouble data in " + variableDisplayName,
                                              "Trouble.png", compared_fig_list)
            if len(troubleSingleList) < len(troubleFullList) :
                functionsToReturn['troubleAvg']   = ((lambda: figures.create_line_plot_figure(troubleSingleList,
                                                                                       "Avg. areas of trouble data in\n" + variableDisplayName)),
                                                     "avg. trouble data in " + variableDisplayName,
                                                     "TroubleAvg.png", compared_fig_list)
        
        return functionsToReturn

"""
This class creates simple imshow plots of 2D data
"""
class IMShowPlotFunctionFactory (PlottingFunctionFactory) :
    def create_plotting_functions (
                                   self,
                                   
                                   # the most basic data set needed
                                   aData, bData,
                                   variableDisplayName,
                                   epsilon,
                                   goodInAMask, goodInBMask,
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # parameters that are only used if the data can be compared
                                   # point by point
                                   absDiffData=None, rawDiffData=None,
                                   goodInBothMask=None,
                                   troubleMask=None, outsideEpsilonMask=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None
                                   ) :
        """
        This method generates imshow plotting functions for two dimensional data
        and returns them in a dictionary of tupples, where each tupple is in the form:
        
            returnDictionary['descriptive name'] = (function, title, file_name, list_this_figure_should_go_into)
        
        The file name is only the name of the file, not the full path.
        """
        
        assert(aData is not None)
        assert(bData is not None)
        
        functionsToReturn = { }
        
        # make the original data plots
        if ('do_plot_originals' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_originals']) :
            
            assert(goodInAMask is not None)
            assert(aData.shape == goodInAMask.shape)
            
            functionsToReturn['originalA'] = ((lambda: figures.create_simple_figure(aData, variableDisplayName + "\nin File A",
                                                                            invalidMask=~goodInAMask)),
                                              variableDisplayName + " in file a",
                                              "A.png",  original_fig_list)
            
            assert(goodInBMask is not None)
            assert(bData.shape == goodInBMask.shape)
            
            functionsToReturn['originalB'] = ((lambda: figures.create_simple_figure(bData, variableDisplayName + "\nin File B",
                                                                            invalidMask=~goodInBMask)),
                                              variableDisplayName + " in file b",
                                              "B.png",  original_fig_list)
        
        # make the absolute value difference plot
        if ('do_plot_abs_diff' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_abs_diff']) :
            
            assert(absDiffData    is not None)
            assert(goodInBothMask is not None)
            assert(goodInBothMask.shape == absDiffData.shape)
            
            functionsToReturn['diffAbs']   = ((lambda: figures.create_simple_figure(absDiffData,
                                                                            "Absolute value of difference in\n" + variableDisplayName,
                                                                            invalidMask=~goodInBothMask)),
                                              "absolute value of difference in " + variableDisplayName,
                                              "AbsDiff.png", compared_fig_list)
        # make the subtractive difference plot
        if ('do_plot_sub_diff' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_sub_diff']) :
            
            assert(rawDiffData    is not None)
            assert(goodInBothMask is not None)
            assert(goodInBothMask.shape == rawDiffData.shape)
            
            functionsToReturn['diffSub']   = ((lambda: figures.create_simple_figure(rawDiffData,
                                                                            "Value of (Data File B - Data File A) for\n" + variableDisplayName,
                                                                            invalidMask=~goodInBothMask)),
                                              "the difference in " + variableDisplayName,
                                              "Diff.png",    compared_fig_list)
        # make the trouble data plot
        if ('do_plot_trouble' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_trouble']) :
            
            assert(goodInAMask is not None)
            assert(goodInBMask is not None)
            assert(goodInAMask.shape == goodInBMask.shape)
            assert(aData.shape       == bData.shape)
            assert(aData.shape       == goodInAMask.shape)
            assert(troubleMask is not None)
            assert(troubleMask.shape == aData.shape)
            
            
            # this is not an optimal solution, but we need to have at least somewhat valid data at any mismatched points so
            # that our plot won't be totally destroyed by missing or non-finite data from B
            bDataCopy = bData[:]
            tempMask = goodInAMask & (~goodInBMask) 
            bDataCopy[tempMask] = aData[tempMask]
            functionsToReturn['trouble']   = ((lambda: figures.create_simple_figure(bDataCopy, "Areas of trouble data in\n" + variableDisplayName,
                                                                            invalidMask=~(goodInAMask | goodInBMask), tagData=troubleMask,
                                                                            colorMap=mediumGrayColorMap)),
                                              "trouble data in " + variableDisplayName,
                                              "Trouble.png", compared_fig_list)
        
        return functionsToReturn


if __name__=='__main__':
    import doctest
    doctest.testmod()
