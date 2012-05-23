#!/usr/bin/env python
# encoding: utf-8
"""
Plotting figure generation functions.

Created by evas Dec 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import matplotlib.colors as colors
import matplotlib.cm     as colormapinfo

import logging
import random as random
import numpy as np
from numpy import ma 

import glance.graphics as maps
import glance.delta    as delta
import glance.figures  as figures

LOG = logging.getLogger(__name__)

# a constant for the larger size dpi
fullSizeDPI = 150 # 200
# a constant for the thumbnail size dpi
thumbSizeDPI = 50

# ********************* Section of utility functions ***********************

# a method to stop people from calling my fake abstract methods
def _abstract( ) :
    raise NotImplementedError('Method must be implemented in subclass.')

# figure out the bounding axes for the display given a set of
# longitude and latitude and possible a mask of good values
def get_visible_axes(longitudeData, latitudeData, goodMask) :
    
    # calculate the bounding range for the display
    # this is in the form [longitude min, longitude max, latitude min, latitude max]
    visibleAxes = [delta.min_with_mask(longitudeData, goodMask),
                   delta.max_with_mask(longitudeData, goodMask),
                   delta.min_with_mask(latitudeData,  goodMask),
                   delta.max_with_mask(latitudeData,  goodMask)]
    
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
    
    # TODO, the cylindrical projections now have some sort of bizarre behavior where they
    # show crazy things in the empty space in soundings data. instead we are moving back to
    # conics for the moment and additional testing has been added to widen the viewing window
    #projToUse = 'lcc'
    
    # how big is the field of view?
    longitudeRange  = abs(boundingAxes[1] - boundingAxes[0])
    latitudeRange   = abs(boundingAxes[3] - boundingAxes[2])
    # chose the projection based on the range we have to cover
    if (longitudeRange > 180) :
        projToUse = 'cyl' # use a equidistant cylindrical projection to show the whole world
    elif (longitudeRange > 100) or (latitudeRange > 70) :
        projToUse = 'ortho' # use an orthographic projection to show about half the globe
    
    return projToUse

def _make_shared_range(aData, goodInAMask, bData, goodInBMask, shouldUseSharedRangeForOriginal=False) :
    """
    If a shared range is desired, figure out what the shared range including all the data in
    both sets is and return it. If it is not desired, return None.
    """
    
    # figure out the shared range for A and B's data, by default don't share a range
    sharedRange = None
    if (shouldUseSharedRangeForOriginal) :
        sharedRange = figures._make_range(aData, goodInAMask, 50, offset_to_range=figures.offsetToRange,
                                   data_b=bData, valid_b_mask=goodInBMask)
    
    return sharedRange

def _make_axis_and_basemap(lonLatDataDict, goodInAMask, goodInBMask=None, variableDisplayName=None) :
    """
    Determine the appropriate axes for the given data (in longitude and latitude) and create the appropriate basemap.
    
    fullAxis, baseMapInstance, sharedRange = _make_axis_and_basemap(lonLatDataDict, goodInAMask, goodInBMask)
    """
    
    nameMessage = ''
    if variableDisplayName is not None:
        nameMessage = " (" + variableDisplayName + ")"
    
    # figure out the bounding axis
    aAxis        = get_visible_axes(lonLatDataDict['a']['lon'], lonLatDataDict['a']['lat'], goodInAMask)
    fullAxis     = aAxis
    bAxis        = None
    if goodInBMask is not None :
        bAxis    = get_visible_axes(lonLatDataDict['b']['lon'], lonLatDataDict['b']['lat'], goodInBMask)
        fullAxis = [min(aAxis[0], bAxis[0]), max(aAxis[1], bAxis[1]),
                    min(aAxis[2], bAxis[2]), max(aAxis[3], bAxis[3])]
    else :
        LOG.debug("No file b valid mask provided, using visible axes boundaries from file a.")
    
    LOG.debug("Visible axes for file A variable data" + nameMessage + " are: " + str(aAxis))
    if goodInBMask is not None : 
        LOG.debug("Visible axes for file B variable data" + nameMessage + " are: " + str(bAxis))
        LOG.debug("Visible axes shared for both file's variable data" + nameMessage + " are: " + str(fullAxis))
    
    if (fullAxis[0] is None) or (fullAxis[1] is None) or (fullAxis[2] is None) or (fullAxis[3] is None) :
        LOG.warn("Unable to display figures for variable" + nameMessage + " because of inability to identify" +
                 " usable bounding longitude and latitude range on the earth. Bounding range that was identified:" + str(fullAxis))
        return # TODO, the figures need to be disabled from the report and possibly a warning on the report? throw exception instead?
    
    # create our basemap
    LOG.info('\t\tloading base map data')
    lonToUse = lonLatDataDict['a']['lon']
    latToUse = lonLatDataDict['a']['lat']
    if goodInBMask is not None :
        lonToUse = lonLatDataDict['common']['lon']
        latToUse = lonLatDataDict['common']['lat']
    baseMapInstance, fullAxis = maps.create_basemap(lonToUse, latToUse,
                                                    fullAxis, select_projection(fullAxis))
    
    return fullAxis, baseMapInstance

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
                                   aData, bData, # these should be data objects from glance.data
                                   variableDisplayName,
                                   epsilon,
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # a comparison of the data if the data comparison info is needed
                                   # this should be a DiffInfoObject from glance.data
                                   differences=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None,
                                   binName=None,  tupleName=None,
                                   
                                   # the optional epsilon for comparison of a percent of A
                                   epsilonPercent=None,
                                   # the optional units for display
                                   units_a=None, units_b=None
                                   
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
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # a comparison of the data if the data comparison info is needed
                                   # this should be a DiffInfoObject from glance.data
                                   differences=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None,
                                   binName=None,  tupleName=None,
                                   
                                   # the optional epsilon for comparison of a percent of A
                                   epsilonPercent=None,
                                   # the optional units for display
                                   units_a=None, units_b=None
                                   
                                   ) :
        
        #TODO, for the moment, unpack these values into local variables; FUTURE use the objects directly and as needed
        goodInAMask        = differences.a_data_object.masks.valid_mask
        goodInBMask        = differences.b_data_object.masks.valid_mask
        rawDiffData        = differences.diff_data_object.data
        absDiffData        = np.abs(rawDiffData) # we also want to show the distance between our two, not just which one's bigger/smaller
        goodInBothMask     = differences.diff_data_object.masks.valid_mask
        mismatchMask       = differences.diff_data_object.masks.mismatch_mask
        outsideEpsilonMask = differences.diff_data_object.masks.outside_epsilon_mask
        aData              = aData.data
        bData              = bData.data
        
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
                                                                         True, units=units_a)),
                                              "histogram of the amount of difference in " + variableDisplayName,
                                              "Hist.png", compared_fig_list)
        # make the scatter plot
        if ('do_plot_scatter' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_scatter']) :
            
            assert(aData.shape    == bData.shape)
            assert(bData.shape    == goodInBothMask.shape)
            assert(goodInBothMask.shape == outsideEpsilonMask.shape)
            
            # TODO, if there's an epsilon percent, how should the epsilon lines be drawn?
            functionsToReturn['scatter']   = ((lambda : figures.create_scatter_plot(aData[goodInBothMask], bData[goodInBothMask],
                                                                            "Value in File A vs Value in File B",
                                                                            "File A Value", "File B Value",
                                                                            outsideEpsilonMask[goodInBothMask],
                                                                            epsilon, units_x=units_a, units_y=units_b)),
                                              "scatter plot of file a values vs file b values for " + variableDisplayName,
                                              "Scatter.png", compared_fig_list)
        
        # make a hexplot, which is like a scatter plot with density
        if ('do_plot_hex' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_hex']) :
            
            assert(aData.shape == bData.shape)
            assert(bData.shape == goodInBothMask.shape)
            
            functionsToReturn['scatterD']  = ((lambda : figures.create_hexbin_plot(aData[goodInBothMask], bData[goodInBothMask],
                                                                                   "Value in File A vs Value in File B",
                                                                                   "File A Value", "File B Value", epsilon,
                                                                                   units_x=units_a, units_y=units_b)),
                                              "density of file a values vs file b values for " + variableDisplayName,
                                              "Hex.png", compared_fig_list)
        
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
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # a comparison of the data if the data comparison info is needed
                                   # this should be a DiffInfoObject from glance.data
                                   differences=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None,
                                   binName=None,  tupleName=None,
                                   
                                   # the optional epsilon for comparison of a percent of A
                                   epsilonPercent=None,
                                   # the optional units for display
                                   units_a=None, units_b=None
                                   
                                   ) :
        
        #TODO, for the moment, unpack these values into local variables; FUTURE use the objects directly and as needed
        goodInAMask        = differences.a_data_object.masks.valid_mask
        goodInBMask        = differences.b_data_object.masks.valid_mask
        rawDiffData        = differences.diff_data_object.data
        absDiffData        = np.abs(rawDiffData) # we also want to show the distance between our two, not just which one's bigger/smaller
        goodInBothMask     = differences.diff_data_object.masks.valid_mask
        mismatchMask       = differences.diff_data_object.masks.mismatch_mask
        outsideEpsilonMask = differences.diff_data_object.masks.outside_epsilon_mask
        aData              = aData.data
        bData              = bData.data
        
        # the default for plotting geolocated data
        mappedPlottingFunction = figures.create_mapped_figure
        
        functionsToReturn = { }
        
        assert(lonLatDataDict is not None)
        assert(goodInAMask    is not None)
        assert(goodInBMask    is not None)
        
        # TODO, do I also need to encorporate the lon/lat invalid masks with the good masks?
        fullAxis, baseMapInstance = _make_axis_and_basemap(lonLatDataDict,
                                                           goodInAMask, goodInBMask,
                                                           variableDisplayName)
        sharedRange = _make_shared_range(aData, goodInAMask,
                                         bData, goodInBMask,
                                         shouldUseSharedRangeForOriginal)
        
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
                                                                               dataRangeColors=dataColors,
                                                                               units=units_a)),
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
                                                                               dataRangeColors=dataColors,
                                                                               units=units_b)),
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
                                                                               invalidMask=(~goodInBothMask),
                                                                               units=units_a)),
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
                                                                               invalidMask=(~goodInBothMask),
                                                                               units=units_a)),
                                              "the difference in " + variableDisplayName,
                                              "Diff.png",    compared_fig_list)
        # make the mismatch data plot
        if ('do_plot_mismatch' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_mismatch']) :
            
            assert(aData.shape == bData.shape)
            assert(goodInAMask.shape == goodInBMask.shape)
            assert('common' in lonLatDataDict)
            assert('lat'    in lonLatDataDict['common'])
            assert('lon'    in lonLatDataDict['common'])
            assert(lonLatDataDict['common']['lat'].shape == lonLatDataDict['common']['lon'].shape)
            assert(lonLatDataDict['common']['lon'].shape == aData.shape)
            
            # this is not an optimal solution, but we need to have at least somewhat valid data at any mismatched points so
            # that our plot won't be totally destroyed by missing or non-finite data from B
            bDataCopy = np.array(bData)
            tempMask = goodInAMask & (~goodInBMask) 
            bDataCopy[tempMask] = aData[tempMask]
            functionsToReturn['mismatch']   = ((lambda : mappedPlottingFunction(bDataCopy, 
                                                                               lonLatDataDict['common']['lat'], 
                                                                               lonLatDataDict['common']['lon'],
                                                                               baseMapInstance, fullAxis,
                                                                               ("Areas of mismatch data in\n" + variableDisplayName),
                                                                               invalidMask=(~(goodInAMask | goodInBMask)),
                                                                               colorMap=figures.MEDIUM_GRAY_COLOR_MAP, tagData=mismatchMask,
                                                                               dataRanges=dataRanges,
                                                                               dataRangeNames=dataRangeNames,
                                                                               units=units_a)), # TODO, does this need modification?
                                              "mismatch data in " + variableDisplayName,
                                              "Mismatch.png", compared_fig_list)
        
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
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # a comparison of the data if the data comparison info is needed
                                   # this should be a DiffInfoObject from glance.data
                                   differences=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None,
                                   binName=None,  tupleName=None,
                                   
                                   # the optional epsilon for comparison of a percent of A
                                   epsilonPercent=None,
                                   # the optional units for display
                                   units_a=None, units_b=None
                                   
                                   ) :
        
        #TODO, for the moment, unpack these values into local variables; FUTURE use the objects directly and as needed
        goodInAMask        = differences.a_data_object.masks.valid_mask
        goodInBMask        = differences.b_data_object.masks.valid_mask
        rawDiffData        = differences.diff_data_object.data
        absDiffData        = np.abs(rawDiffData) # we also want to show the distance between our two, not just which one's bigger/smaller
        goodInBothMask     = differences.diff_data_object.masks.valid_mask
        mismatchMask       = differences.diff_data_object.masks.mismatch_mask
        outsideEpsilonMask = differences.diff_data_object.masks.outside_epsilon_mask
        aData              = aData.data
        bData              = bData.data
        
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
        
        # TODO, do I also need to encorporate the lon/lat invalid masks with the good masks?
        fullAxis, baseMapInstance = _make_axis_and_basemap(lonLatDataDict, goodInAMask, goodInBMask,
                                                           variableDisplayName=variableDisplayName)
        
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
                                                                               uData=aUData, vData=aVData,
                                                                               units=units_a)),
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
                                                                               uData=bUData, vData=bVData,
                                                                               units=units_b)),
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
                                                                                   uData=diffUData, vData=diffVData,
                                                                                   units=units_a)),
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
                                                                                   uData=diffUData, vData=diffVData,
                                                                                   units=units_a)),
                                                  "the difference in " + variableDisplayName,
                                                  "Diff.png",    compared_fig_list)
            # make the mismatch data plot
            if ('do_plot_mismatch' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_mismatch']) :
                
                assert(aData.shape == bData.shape)
                assert(goodInAMask.shape == goodInBMask.shape)
                assert('common' in lonLatDataDict)
                assert('lat'    in lonLatDataDict['common'])
                assert('lon'    in lonLatDataDict['common'])
                assert(lonLatDataDict['common']['lat'].shape == lonLatDataDict['common']['lon'].shape)
                assert(lonLatDataDict['common']['lon'].shape == aData.shape)
                
                # this is not an optimal solution, but we need to have at least somewhat valid data at any mismatched points so
                # that our plot won't be totally destroyed by missing or non-finite data from B
                bDataCopy = np.array(bData)
                tempMask = goodInAMask & (~goodInBMask) 
                bDataCopy[tempMask] = aData[tempMask]
                functionsToReturn['mismatch']   = ((lambda : mappedPlottingFunction(bDataCopy, 
                                                                                   lonLatDataDict['common']['lat'], 
                                                                                   lonLatDataDict['common']['lon'],
                                                                                   baseMapInstance, fullAxis,
                                                                                   ("Areas of mismatch data in\n" + variableDisplayName),
                                                                                   invalidMask=(~(goodInAMask | goodInBMask)),
                                                                                   colorMap=figures.MEDIUM_GRAY_COLOR_MAP, tagData=mismatchMask,
                                                                                   dataRanges=dataRanges,
                                                                                   dataRangeNames=dataRangeNames,
                                                                                   # TODO, does this need modification?
                                                                                   uData=bUData, vData=bVData,
                                                                                   units=units_a)), 
                                                  "mismatch data in " + variableDisplayName,
                                                  "Mismatch.png", compared_fig_list)
        
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
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # a comparison of the data if the data comparison info is needed
                                   # this should be a DiffInfoObject from glance.data
                                   differences=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None,
                                   binName=None,  tupleName=None,
                                   
                                   # the optional epsilon for comparison of a percent of A
                                   epsilonPercent=None,
                                   # the optional units for display
                                   units_a=None, units_b=None
                                   
                                   ) :
        """
        This method generates line plotting functions for one dimensional data
        and returns them in a dictionary of tupples, where each tupple is in the form:
        
            returnDictionary['descriptive name'] = (function, title, file_name, list_this_figure_should_go_into)
        
        The file name is only the name of the file, not the full path.
        """
        
        #TODO, for the moment, unpack these values into local variables; FUTURE use the objects directly and as needed
        goodInAMask        = differences.a_data_object.masks.valid_mask
        goodInBMask        = differences.b_data_object.masks.valid_mask
        rawDiffData        = differences.diff_data_object.data
        absDiffData        = np.abs(rawDiffData) # we also want to show the distance between our two, not just which one's bigger/smaller
        goodInBothMask     = differences.diff_data_object.masks.valid_mask
        mismatchMask       = differences.diff_data_object.masks.mismatch_mask
        outsideEpsilonMask = differences.diff_data_object.masks.outside_epsilon_mask
        aData              = aData.data
        bData              = bData.data
        
        assert(aData is not None)
        assert(bData is not None)
        assert(len(aData.shape) is 1)
        assert(aData.shape == bData.shape)
        
        # make all our data sets for plotting ahead of time for simplicity
        aList = [(aData, ~goodInAMask, 'r', 'A data', None, units_a)]
        bList = [(bData, ~goodInBMask, 'b', 'B data', None, units_b)]
        absDiffList = [(absDiffData, ~goodInBothMask, '', 'abs. diff. data', None, units_a)] # todo, should this be a units?
        subDiffList = [(rawDiffData, ~goodInBothMask, '', 'sub. diff. data', None, units_a)] # todo, should this be a units?
        
        mismatchList   = [(aData, ~goodInAMask, 'r', 'A data', mismatchMask, units_a),
                          (bData, ~goodInBMask, 'b', 'B data', mismatchMask, units_b)]
        
        functionsToReturn = { }
        
        # make the original data plots
        if ('do_plot_originals' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_originals']) :
            
            
            functionsToReturn['original'] = ((lambda: figures.create_line_plot_figure((aList + bList),
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
        
        # make the absolute value difference plot
        if ('do_plot_abs_diff' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_abs_diff']) :
            functionsToReturn['diffAbs']   = ((lambda: figures.create_line_plot_figure(absDiffList,
                                                                               "Absolute value of difference in\n" + variableDisplayName)),
                                              "absolute value of difference in " + variableDisplayName,
                                              "AbsDiff.png", compared_fig_list)
        # make the subtractive difference plot
        if ('do_plot_sub_diff' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_sub_diff']) :
            functionsToReturn['diffSub']   = ((lambda: figures.create_line_plot_figure(subDiffList,
                                                                               "Value of (Data File B - Data File A) for\n"
                                                                               + variableDisplayName)),
                                              "the difference in " + variableDisplayName,
                                              "Diff.png",    compared_fig_list)
        
        # make the mismatch data plot
        if ('do_plot_mismatch' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_mismatch']) :
            functionsToReturn['mismatch']   = ((lambda: figures.create_line_plot_figure(mismatchList,
                                                                               "Areas of mismatch data in\n" + variableDisplayName)),
                                              "mismatch data in " + variableDisplayName,
                                              "Mismatch.png", compared_fig_list)
        
        return functionsToReturn

class BinTupleAnalysisFunctionFactory (PlottingFunctionFactory) :
    def create_plotting_functions (
                                   self,
                                   
                                   # the most basic data set needed
                                   aData, bData,
                                   variableDisplayName,
                                   epsilon,
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # a comparison of the data if the data comparison info is needed
                                   # this should be a DiffInfoObject from glance.data
                                   differences=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None,
                                   binName=None,  tupleName=None,
                                   
                                   # the optional epsilon for comparison of a percent of A
                                   epsilonPercent=None,
                                   # the optional units for display
                                   units_a=None, units_b=None
                                   
                                   ) :
        """
        This method generates histogram and sample line plot functions for complex three dimensional data
        and returns them in a dictionary of tupples, where each tupple is in the form:
        
            returnDictionary['descriptive name'] = (function, title, file_name, list_this_figure_should_go_into)
        
        The file name is only the name of the file, not the full path.
        """
        
        #TODO, for the moment, unpack these values into local variables; FUTURE use the objects directly and as needed
        goodInAMask        = differences.a_data_object.masks.valid_mask
        goodInBMask        = differences.b_data_object.masks.valid_mask
        rawDiffData        = differences.diff_data_object.data
        absDiffData        = np.abs(rawDiffData) # we also want to show the distance between our two, not just which one's bigger/smaller
        goodInBothMask     = differences.diff_data_object.masks.valid_mask
        mismatchMask       = differences.diff_data_object.masks.mismatch_mask
        outsideEpsilonMask = differences.diff_data_object.masks.outside_epsilon_mask
        aData              = aData.data
        bData              = bData.data
        
        # confirm that our a and b data are minimally ok
        assert(aData is not None)
        assert(bData is not None)
        assert(len(aData.shape) >= 2)
        assert(aData.shape == bData.shape)
        
        # confirm that our bin and tuple indexes are valid
        assert(binIndex   is not None)
        assert(tupleIndex is not None)
        assert(binIndex   is not tupleIndex)
        assert(binIndex   < len(aData.shape))
        assert(tupleIndex < len(aData.shape))
        
        # reorder and reshape our data into the [bin][case][tuple] form
        reorderMapObject   = delta.BinTupleMapping(aData.shape,
                                                   binIndexNumber=binIndex,
                                                   tupleIndexNumber=tupleIndex)
        aData              = reorderMapObject.reorder_for_bin_tuple(aData)
        bData              = reorderMapObject.reorder_for_bin_tuple(bData)
        goodInAMask        = reorderMapObject.reorder_for_bin_tuple(goodInAMask)
        goodInBMask        = reorderMapObject.reorder_for_bin_tuple(goodInBMask)
        absDiffData        = reorderMapObject.reorder_for_bin_tuple(absDiffData)
        rawDiffData        = reorderMapObject.reorder_for_bin_tuple(rawDiffData)
        goodInBothMask     = reorderMapObject.reorder_for_bin_tuple(goodInBothMask)
        mismatchMask       = reorderMapObject.reorder_for_bin_tuple(mismatchMask)
        outsideEpsilonMask = reorderMapObject.reorder_for_bin_tuple(outsideEpsilonMask)
        
        # our list of functions that will later create the plots
        functionsToReturn = { }
        
        if (aData.size <= 0) :
            return functionsToReturn
        
        # create the scatter plot with colors for each section
        scatterPlotList = [ ]
        tempColorMap = colormapinfo.get_cmap('jet', rawDiffData.shape[0])
        for binNumber in range(rawDiffData.shape[0]) :
            tempColor = tempColorMap(binNumber)
            if len(tempColor) > 3 :
                tempColor = tempColor[:3]
            scatterPlotList.append(((aData[binNumber][goodInBothMask[binNumber]]).ravel(),
                                    (bData[binNumber][goodInBothMask[binNumber]]).ravel(), None,
                                    colors.rgb2hex(tempColor), None, 'bin ' + str(binNumber + 1), None))
        functionsToReturn['multi-scatter'] = ((lambda : figures.create_complex_scatter_plot(scatterPlotList,
                                                                        "Value in File A vs Value in File B, Colored by Bin",
                                                                        "File A Value", "File B Value",
                                                                        epsilon, units_x=units_a, units_y=units_b)),
                                          "scatter plot of file a values vs file b values for " + variableDisplayName + " by bin",
                                          "MultiScatter.png", compared_fig_list)
        
        # for each of the bins, make the rms histogram data
        numHistogramSections = 7 # TODO at some point make this a user controlled setting
        for binNumber in range(rawDiffData.shape[0]) :
            
            new_list = [ ]
            compared_fig_list.append(new_list)
            
            # figure out all the rms diff values for the various cases
            rmsDiffValues = np.zeros(rawDiffData.shape[1])
            for caseNumber in range(rawDiffData.shape[1]) :
                rmsDiffValues[caseNumber] = delta.calculate_root_mean_square(rawDiffData[binNumber][caseNumber],
                                                                             goodInBothMask[binNumber][caseNumber])
            
            # make the basic histogram for this binNumber
            dataForHistogram = rmsDiffValues[np.isfinite(rmsDiffValues)] # remove any invalid data "nan" values
            if ('do_plot_histogram' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_histogram']) :
                def make_histogram(binNumber=binNumber, dataForHistogram=dataForHistogram):
                    return figures.create_histogram(dataForHistogram, numHistogramSections,
                                                                             ("RMS Diff. in " + variableDisplayName +
                                                                              "\nfor " + binName + " # " + str(binNumber + 1)),
                                                                             ('RMS Difference across ' + tupleName + ' dimension'),
                                                                             ('Number of Cases with a Given RMS Diff.'),
                                                                             True, units=units_a)
                functionsToReturn[str(binNumber + 1) + 'histogram'] = (make_histogram,
                                                  "histogram of rms differences in " + variableDisplayName,
                                                  str(binNumber + 1) + "Hist.png", new_list)
            
            # we will need to be able to mask out the non-finite data
            tempFiniteMap = np.isfinite(rmsDiffValues)
            
            # figure out the min/max rms diff values
            minRMSDiff = np.min(rmsDiffValues[tempFiniteMap])
            maxRMSDiff = np.max(rmsDiffValues[tempFiniteMap])
            
            # sort the cases by their rms diff values
            counts = np.zeros(numHistogramSections)
            histogramSections = { }
            histogramSectionLimits = np.linspace(minRMSDiff, maxRMSDiff, numHistogramSections + 1)
            histogramSectionLimits[0] = histogramSectionLimits[0] - 0.00000001
            for caseNumber in range(rmsDiffValues.size) :
                
                # check each of the sections to see which one it falls in
                for limitIndex in range(histogramSectionLimits.size - 1) :
                    
                    # if it falls in this section, add it's case number index to the list for this section
                    if ( (rmsDiffValues[caseNumber] >  histogramSectionLimits[limitIndex]) and
                         (rmsDiffValues[caseNumber] <= histogramSectionLimits[limitIndex + 1]) ) :
                        
                        if limitIndex not in histogramSections :
                            histogramSections[limitIndex] = [ ]
                        
                        histogramSections[limitIndex].append(caseNumber)
            
            # select example cases for the histogram
            random.seed('test') # TODO, seed with something else?
            for section in sorted(histogramSections.keys()) :
                listOfCases = histogramSections[section]
                caseNumber  = listOfCases[random.randint(0, len(listOfCases) - 1)]
                
                # make lineplot functions for the example cases
                caseIndexes = reorderMapObject.determine_case_indecies(caseNumber)
                caseNumText = ''
                for caseIndex in caseIndexes :
                    caseNumText = caseNumText + str(caseIndex)
                dataList = [(aData[binNumber][caseNumber], ~goodInAMask[binNumber][caseNumber], 'r', 'A case', None, units_a),
                            (bData[binNumber][caseNumber], ~goodInBMask[binNumber][caseNumber], 'b', 'B case', None, units_b)]
                def make_lineplot(data=dataList, binNumber=binNumber, caseNumberText=caseNumText):
                    return figures.create_line_plot_figure(data,
                                                           variableDisplayName + " in both files" + "\n" + "for "
                                                           + binName + " # " + str(binNumber + 1) + " and case # "
                                                           + caseNumberText)
                dataDiff = aData[binNumber][caseNumber] - bData[binNumber][caseNumber]
                maskDiff = ~goodInAMask[binNumber][caseNumber] | ~goodInBMask[binNumber][caseNumber]
                def make_diffplot(data=dataDiff, badMask=maskDiff, binNumber=binNumber, caseNumberText=caseNumText, units=units_a):
                    return figures.create_line_plot_figure([(data, badMask, 'm', 'A - B', None, units)],
                                                           "Value of " + variableDisplayName + " in File A - the value in File B\n" +
                                                           "for " + binName + " # " + str(binNumber + 1) + " and case # " + caseNumberText)
                    
                functionsToReturn[str(binNumber + 1) + 'sample' + str(section + 1) + '.AB.png'] = (make_lineplot,
                                                                                           variableDisplayName + " sample in both files",
                                                                                           str(binNumber + 1) + 'sample' + str(section + 1) + '.AB.png',
                                                                                           new_list)
                functionsToReturn[str(binNumber + 1) + 'sample' + str(section + 1) + 'ABdiff.png] '] = (make_diffplot,
                                                                                           variableDisplayName + " difference between files",
                                                                                           str(binNumber + 1) + 'sample' + str(section + 1) + '.ABdiff.png',
                                                                                           new_list)
        
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
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # a comparison of the data if the data comparison info is needed
                                   # this should be a DiffInfoObject from glance.data
                                   differences=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None,
                                   binName=None,  tupleName=None,
                                   
                                   # the optional epsilon for comparison of a percent of A
                                   epsilonPercent=None,
                                   # the optional units for display
                                   units_a=None, units_b=None
                                   
                                   ) :
        """
        This method generates imshow plotting functions for two dimensional data
        and returns them in a dictionary of tupples, where each tupple is in the form:
        
            returnDictionary['descriptive name'] = (function, title, file_name, list_this_figure_should_go_into)
        
        The file name is only the name of the file, not the full path.
        """
        
        #TODO, for the moment, unpack these values into local variables; FUTURE use the objects directly and as needed
        goodInAMask        = differences.a_data_object.masks.valid_mask
        goodInBMask        = differences.b_data_object.masks.valid_mask
        rawDiffData        = differences.diff_data_object.data
        absDiffData        = np.abs(rawDiffData) # we also want to show the distance between our two, not just which one's bigger/smaller
        goodInBothMask     = differences.diff_data_object.masks.valid_mask
        mismatchMask       = differences.diff_data_object.masks.mismatch_mask
        outsideEpsilonMask = differences.diff_data_object.masks.outside_epsilon_mask
        aData              = aData.data
        bData              = bData.data
        
        assert(aData is not None)
        assert(bData is not None)
        
        functionsToReturn = { }
        
        sharedRange = _make_shared_range(aData, goodInAMask,
                                         bData, goodInBMask,
                                         shouldUseSharedRangeForOriginal)
        
        # make the original data plots
        if ('do_plot_originals' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_originals']) :
            
            assert(goodInAMask is not None)
            assert(aData.shape == goodInAMask.shape)
            
            functionsToReturn['originalA'] = ((lambda: figures.create_simple_figure(aData, variableDisplayName + "\nin File A",
                                                                            invalidMask=~goodInAMask, colorbarLimits=sharedRange, 
                                                                            units=units_a)),
                                              variableDisplayName + " in file a",
                                              "A.png",  original_fig_list)
            
            assert(goodInBMask is not None)
            assert(bData.shape == goodInBMask.shape)
            
            functionsToReturn['originalB'] = ((lambda: figures.create_simple_figure(bData, variableDisplayName + "\nin File B",
                                                                            invalidMask=~goodInBMask, colorbarLimits=sharedRange, 
                                                                            units=units_b)),
                                              variableDisplayName + " in file b",
                                              "B.png",  original_fig_list)
        
        # make the absolute value difference plot
        if ('do_plot_abs_diff' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_abs_diff']) :
            
            assert(absDiffData    is not None)
            assert(goodInBothMask is not None)
            assert(goodInBothMask.shape == absDiffData.shape)
            
            functionsToReturn['diffAbs']   = ((lambda: figures.create_simple_figure(absDiffData,
                                                                            "Absolute value of difference in\n" + variableDisplayName,
                                                                            invalidMask=~goodInBothMask, units=units_a)),
                                              "absolute value of difference in " + variableDisplayName,
                                              "AbsDiff.png", compared_fig_list)
        # make the subtractive difference plot
        if ('do_plot_sub_diff' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_sub_diff']) :
            
            assert(rawDiffData    is not None)
            assert(goodInBothMask is not None)
            assert(goodInBothMask.shape == rawDiffData.shape)
            
            functionsToReturn['diffSub']   = ((lambda: figures.create_simple_figure(rawDiffData,
                                                                            "Value of (Data File B - Data File A) for\n" + variableDisplayName,
                                                                            invalidMask=~goodInBothMask, units=units_a)),
                                              "the difference in " + variableDisplayName,
                                              "Diff.png",    compared_fig_list)
        # make the mismatch data plot
        if ('do_plot_mismatch' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_mismatch']) :
            
            assert(goodInAMask is not None)
            assert(goodInBMask is not None)
            assert(goodInAMask.shape == goodInBMask.shape)
            assert(aData.shape       == bData.shape)
            assert(aData.shape       == goodInAMask.shape)
            assert(mismatchMask is not None)
            assert(mismatchMask.shape == aData.shape)
            
            
            # this is not an optimal solution, but we need to have at least somewhat valid data at any mismatched points so
            # that our plot won't be totally destroyed by missing or non-finite data from B
            bDataCopy = np.array(bData)
            tempMask = goodInAMask & (~goodInBMask) 
            bDataCopy[tempMask] = aData[tempMask]
            functionsToReturn['mismatch']   = ((lambda: figures.create_simple_figure(bDataCopy, "Areas of mismatch data in\n" + variableDisplayName,
                                                                            invalidMask=~(goodInAMask | goodInBMask), tagData=mismatchMask,
                                                                            colorMap=figures.MEDIUM_GRAY_COLOR_MAP, units=units_a)),
                                              "mismatch data in " + variableDisplayName,
                                              "Mismatch.png", compared_fig_list)
        
        return functionsToReturn

# ********** below here are the inspection plotting functions **********
# note: for all of these the epsilons are meaningless, for some the bData has optional effects that won't come up when they're used for inspection reports

"""
This class creates the most basic of histogram plots based on one data set.
"""
class DataHistogramPlotFunctionFactory (PlottingFunctionFactory) :
    def create_plotting_functions (
                                   self,
                                   
                                   # the most basic data set needed
                                   aData, bData, # Note, bData is not used
                                   variableDisplayName,
                                   epsilon,
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # no comparison is needed for this call, should be None
                                   differences=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None,
                                   binName=None,  tupleName=None,
                                   
                                   # the optional epsilon for comparison of a percent of A
                                   epsilonPercent=None,
                                   # the optional units for display
                                   units_a=None, units_b=None
                                   
                                   ) :
        
        functionsToReturn = { }
        
        # right now this function is not intended to handle both types of data; in the FUTURE this may change
        assert (aData is not None)
        assert (bData is     None)
        
        # make the histogram plot; TODO, for now reuse this setting in FUTURE, should a new one be made to control this?
        if ('do_plot_histogram' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_histogram']) :
            
            goodInAMask = aData.masks.valid_mask
            assert(goodInAMask.shape == aData.data.shape)
            
            # setup the data bins for the histogram
            numBinsToUse = 50
            valuesForHist = aData.data[goodInAMask]
            functionsToReturn['histogram'] = ((lambda : figures.create_histogram(valuesForHist, numBinsToUse,
                                                                         ("Values of\n" + variableDisplayName),
                                                                         ('Value of Data Point'),
                                                                         ('Number of Data Points with a Given Value'),
                                                                         True, units=units_a)),
                                              "histogram of the values in " + variableDisplayName,
                                              "HistA.png", original_fig_list)
        
        return functionsToReturn

"""
This class creates a simple imshow plot of 2D data
"""
class InspectIMShowPlotFunctionFactory (PlottingFunctionFactory) :
    def create_plotting_functions (
                                   self,
                                   
                                   # the most basic data set needed
                                   aData, bData, # bData is not expected and will not be used
                                   variableDisplayName,
                                   epsilon,
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # no comparison is needed for this call, should be None
                                   differences=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None,
                                   binName=None,  tupleName=None,
                                   
                                   # the optional epsilon for comparison of a percent of A
                                   epsilonPercent=None,
                                   # the optional units for display
                                   units_a=None, units_b=None
                                   
                                   ) :
        """
        This method generates imshow plotting functions for two dimensional data
        and returns them in a dictionary of tupples, where each tupple is in the form:
        
            returnDictionary['descriptive name'] = (function, title, file_name, list_this_figure_should_go_into)
        
        The file name is only the name of the file, not the full path.
        """
        
        assert(aData      is not None)
        assert(aData.data is not None)
        assert(bData      is     None)
        
        functionsToReturn = { }
        
        # make the original data plots
        if ('do_plot_originals' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_originals']) :
            
            goodInAMask = aData.masks.valid_mask
            
            assert(goodInAMask is not None)
            assert(aData.data.shape == goodInAMask.shape)
            
            functionsToReturn['originalA'] = ((lambda: figures.create_simple_figure(aData.data, variableDisplayName + "\nin File",
                                                                            invalidMask=~goodInAMask,
                                                                            units=units_a)),
                                              variableDisplayName + " in file",
                                              "origA.png",  original_fig_list)
        
        return functionsToReturn

"""
This class creates simple line plots based on simple one dimentional data.
"""
class InspectLinePlotsFunctionFactory (PlottingFunctionFactory) :
    def create_plotting_functions (
                                   self,
                                   
                                   # the most basic data set needed
                                   aData, bData, # bData will not be used
                                   variableDisplayName,
                                   epsilon,
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # no comparison is needed for this call, should be None
                                   differences=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None,
                                   binName=None,  tupleName=None,
                                   
                                   # the optional epsilon for comparison of a percent of A
                                   epsilonPercent=None,
                                   # the optional units for display
                                   units_a=None, units_b=None
                                   
                                   ) :
        """
        This method generates line plotting functions for one dimensional data
        and returns them in a dictionary of tupples, where each tupple is in the form:
        
            returnDictionary['descriptive name'] = (function, title, file_name, list_this_figure_should_go_into)
        
        The file name is only the name of the file, not the full path.
        """
        
        assert(aData      is not None)
        assert(aData.data is not None)
        assert(bData      is     None)
        #assert(len(aData.data.shape) is 1)
        
        functionsToReturn = { }
        
        # make the original data plots
        if ('do_plot_originals' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_originals']) :
            
            aDataTemp =  aData.data
            aMaskTemp = ~aData.masks.valid_mask
            if len(aDataTemp.shape) > 1 :
                aDataTemp = aDataTemp.ravel()
                aMaskTemp = aMaskTemp.ravel()
            aList = [(aDataTemp, aMaskTemp, 'b', 'A data', None, units_a)]
            
            functionsToReturn['originalA'] = ((lambda: figures.create_line_plot_figure(aList,
                                                                                variableDisplayName + "\nin File")),
                                              variableDisplayName + " in file",
                                              "lineA.png",  original_fig_list)
        
        return functionsToReturn

"""
This class creates contour plots mapped onto a region of the earth.
"""
class InspectMappedContourPlotFunctionFactory (PlottingFunctionFactory) :
    def create_plotting_functions (
                                   self,
                                   
                                   # the most basic data set needed
                                   aData, bData, # bData will not be used
                                   variableDisplayName,
                                   epsilon,
                                   doPlotSettingsDict,
                                   
                                   # where the names of the created figures will be stored
                                   original_fig_list, compared_fig_list,
                                   
                                   # parameters that are only needed for geolocated data
                                   lonLatDataDict=None,
                                   
                                   # only used if we are plotting a contour
                                   dataRanges=None, dataRangeNames=None, dataColors=None,
                                   shouldUseSharedRangeForOriginal=True,
                                   
                                   # this isn't needed for this method, should be None
                                   differences=None,
                                   
                                   # only used for plotting quiver data
                                   aUData=None, aVData=None,
                                   bUData=None, bVData=None,
                                   
                                   # only used for line plots
                                   binIndex=None, tupleIndex=None,
                                   binName=None,  tupleName=None,
                                   
                                   # the optional epsilon for comparison of a percent of A
                                   epsilonPercent=None,
                                   # the optional units for display
                                   units_a=None, units_b=None
                                   
                                   ) :
        
        # for simplicity, get the valid mask for a
        goodInAMask = aData.masks.valid_mask
        
        # the default for plotting geolocated data
        mappedPlottingFunction = figures.create_mapped_figure
        
        functionsToReturn = { }
        
        assert(lonLatDataDict is not None)
        assert(goodInAMask    is not None)
        
        print ("lon lat dictionary form: " + str(lonLatDataDict))
        
        # figure out which part of the earth is visible and construct a basemap using that
        fullAxis, baseMapInstance = _make_axis_and_basemap({'a':lonLatDataDict},
                                                           goodInAMask, None, # there is no b mask
                                                           variableDisplayName)
        
        # make the original data plot
        if ('do_plot_originals' not in doPlotSettingsDict) or (doPlotSettingsDict['do_plot_originals']) :
            
            assert('lat' in lonLatDataDict)
            assert('lon' in lonLatDataDict)
            assert(lonLatDataDict['lat'].shape == lonLatDataDict['lon'].shape)
            
            functionsToReturn['originalA'] = ((lambda : mappedPlottingFunction(aData.data,
                                                                               lonLatDataDict['lat'], 
                                                                               lonLatDataDict['lon'],
                                                                               baseMapInstance, fullAxis,
                                                                               (variableDisplayName + "\nin File"),
                                                                               invalidMask=(~goodInAMask),
                                                                               dataRanges=dataRanges,
                                                                               dataRangeNames=dataRangeNames,
                                                                               dataRangeColors=dataColors,
                                                                               units=units_a)),
                                              variableDisplayName + " in file",
                                              "mapA.png",  original_fig_list)
        
        return functionsToReturn

if __name__=='__main__':
    import doctest
    doctest.testmod()
