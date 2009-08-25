#!/usr/bin/env python
# encoding: utf-8
"""
Plotting routines for difference values using matplotlib

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

# these first two lines must stay before the pylab import
import matplotlib
matplotlib.use('Agg') # use the Anti-Grain Geometry rendering engine

from pylab import *

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.ticker import FormatStrFormatter

from PIL import Image

import os, sys, logging
import numpy as np

import glance.graphics as maps
import glance.delta as delta
import glance.report as report

LOG = logging.getLogger(__name__)

# TODO this value is being used to work around a problem with the contourf
# and how it handles range boundaries. Find a better solution if at all possible.
offsetToRange = 0.0000000000000000001

# the value that will denote "bad" longitudes and latitudes
badLonLat = maps.badLonLat

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

# build a scatter plot of the x,y points
def _create_scatter_plot(dataX, dataY, title, xLabel, yLabel, badMask=None, epsilon=None) :
    
    # make the figure
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    # if we have "bad" data to plot, pull it out
    badX = None
    badY = None
    if (badMask != None) :
        badX = dataX[badMask]
        badY = dataY[badMask]
        dataX = dataX[~badMask]
        dataY = dataY[~badMask]
    
    # the scatter plot of the good data 
    axes.plot(dataX, dataY, 'b,', label='within\nepsilon')
    
    # plot the bad data
    numTroublePts = 0
    if (badMask != None) :
        numTroublePts = badX.shape[0]
        LOG.debug('\t\tnumber of trouble points in scatter plot: ' + str(badX.shape[0]))
        if numTroublePts > 0 :
            axes.plot(badX, badY, 'r,', label='outside\nepsilon')
    
    # draw the line for the "perfect fit" 
    xbounds = axes.get_xbound()
    xrange = xbounds[1] - xbounds[0]
    ybounds = axes.get_ybound()
    yrange = ybounds[1] - ybounds[0]
    perfect = [max(xbounds[0], ybounds[0]), min(xbounds[1], ybounds[1])]
    axes.plot(perfect, perfect, 'k--', label='A = B')
    
    # now draw the epsilon bound lines if they are visible and the lines won't be the same as A = B
    if (not (epsilon is None)) and (epsilon > 0.0) and (epsilon < xrange) and (epsilon < yrange):
        # plot the top line
        axes.plot([perfect[0], perfect[1] - epsilon], [perfect[0] + epsilon, perfect[1]], '--', color='#00FF00', label='+/-epsilon')
        # plot the bottom line
        axes.plot([perfect[0] + epsilon, perfect[1]], [perfect[0], perfect[1] - epsilon], '--', color='#00FF00')
    
    # make a key to explain our plot
    # as long as things have been plotted with proper labels they should show up here
    axes.legend(loc=0, markerscale=3.0) # Note: at the moment markerscale doesn't seem to work
    
    # and some informational stuff
    axes.set_title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    
    # format our axes so they display gracefully
    yFormatter = FormatStrFormatter("%4.4g")
    axes.yaxis.set_major_formatter(yFormatter)
    xFormatter = FormatStrFormatter("%4.4g")
    axes.xaxis.set_major_formatter(xFormatter)
    
    return figure

# build a histogram figure of the given data with the given title and number of bins
def _create_histogram(data, bins, title, xLabel, yLabel, displayStats=False) :
    
    # make the figure
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    # the histogram of the data
    n, bins, patches = plt.hist(data, bins)
    
    # format our axes so they display gracefully
    yFormatter = FormatStrFormatter("%3.3g")
    axes.yaxis.set_major_formatter(yFormatter)
    xFormatter = FormatStrFormatter("%.4g")
    axes.xaxis.set_major_formatter(xFormatter)
    
    # and some informational stuff
    axes.set_title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    
    # if stats were passed in, put some of the information on the graph
    # the location is in the form x, y (I think)
    if displayStats :
        # info on the basic stats
        tempMask = ones(data.shape,dtype=bool)
        tempStats = delta.stats(data, tempMask, None)
        medianVal = tempStats['median_diff']
        meanVal = tempStats['mean_diff']
        stdVal = tempStats['std_diff']
        numPts = data.size
        
        # info on the display of our statistics
        xbounds = axes.get_xbound()
        numBinsToUse = len(bins)
        xrange = xbounds[1] - xbounds[0]
        binSize = xrange / float(numBinsToUse)
        
        # build the display string
        statText = ('%1.2e' % numPts) + ' data points'
        statText = statText + '\n' + 'mean: ' + report.make_formatted_display_string(meanVal)
        statText = statText + '\n' + 'median: ' + report.make_formatted_display_string(medianVal)
        statText = statText + '\n' + 'std: ' + report.make_formatted_display_string(stdVal)
        statText = statText + '\n\n' + 'bins: ' + report.make_formatted_display_string(numBinsToUse)
        statText = statText + '\n' + 'bin size ' + report.make_formatted_display_string(binSize)
        
        # figure out where to place the text and put it on the figure
        centerOfDisplay = xbounds[0] + (float(xrange) / 2.0)
        xValToUse = 0.67
        # if most of the values will be on the right, move our text to the left...
        if (medianVal > centerOfDisplay) :
            xValToUse = 0.17
        figtext(xValToUse, 0.60, statText)
    
    return figure

# make a "clean" copy of the longitude or latitude array passed in by replacing all
# values which would fall in the provided "invalid" mask with the default bad-lon-lat
# value, which will keep them from being plotted
def _clean_lon_or_lat_with_mask(lon_or_lat_data, invalid_data_mask):
    clean_data = empty_like(lon_or_lat_data)
    clean_data[~invalid_data_mask] = lon_or_lat_data[~invalid_data_mask]
    clean_data[ invalid_data_mask] = badLonLat
    
    return clean_data

# chose a projection based on the bounding axes that will be shown
def _select_projection(boundingAxes) :
    
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
    projToUse = 'mill'
    
    # how big is the field of view?
    longitudeRange  = abs(boundingAxes[1] - boundingAxes[0])
    latitudeRange   = abs(boundingAxes[3] - boundingAxes[2])
    # chose the projection based on the range we have to cover
    if (longitudeRange > 180) :
        projToUse = 'mill' # use a miller cylindrical projection to show the whole world
    elif (longitudeRange > 100) or (latitudeRange > 70) :
        projToUse = 'ortho' # use an orthographic projection to show about half the globe
    
    return projToUse

# create a figure including our data mapped onto a map at the lon/lat given
# the colorMap parameter can be used to control the colors the figure is drawn in
# if any masks are passed in the tagData list they will be plotted as an overlays
# set on the existing image
def _create_mapped_figure(data, latitude, longitude, baseMapInstance, boundingAxes, title,
                          invalidMask=None, colorMap=None, tagData=None,
                          dataRanges=None, dataRangeNames=None, dataRangeColors=None) :
    
    # make a clean version of our lon/lat
    latitudeClean  = _clean_lon_or_lat_with_mask(latitude,  invalidMask)
    longitudeClean = _clean_lon_or_lat_with_mask(longitude, invalidMask)
    
    # build the plot
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    # build extra info to go to the map plotting function
    kwargs = {}
    
    # figure the range for the color bars
    # this is controllable with the "dataRanges" parameter for discrete data display
    if not (data is None) :
        if dataRanges is None :
            # todo, the use off the offset here is covering a problem with
            # contourf hiding data exactly at the end of the range and should
            # be removed if a better solution can be found
            minVal = delta.min_with_mask(data, invalidMask) - offsetToRange
            maxVal = delta.max_with_mask(data, invalidMask) + offsetToRange
            dataRanges = np.linspace(minVal, maxVal, 50)
        else: # make sure the user range will not discard data TODO, find a better way to handle this
            dataRanges[0] = dataRanges[0] - offsetToRange
            dataRanges[len(dataRanges) - 1] = dataRanges[len(dataRanges) - 1] + offsetToRange
        kwargs['levelsToUse'] = dataRanges
        if dataRangeColors is not None :
            kwargs['colors'] = dataRangeColors # add in the list of colors (may be None)
    
    # if we've got a color map, pass it to the list of things we want to tell the plotting function
    if not (colorMap is None) :
        kwargs['cmap'] = colorMap
    
    # draw our data placed on a map
    #bMap, x, y = maps.mapshow(longitudeClean, latitudeClean, data, boundingAxes, **kwargs)
    maps.draw_basic_features(baseMapInstance, boundingAxes)
    bMap, x, y = maps.show_lon_lat_data(longitudeClean, latitudeClean, baseMapInstance, data, **kwargs)
    
    # and some informational stuff
    axes.set_title(title)
    # show a generic color bar
    doLabelRanges = False
    if not (data is None) :
        cbar = colorbar(format='%.3f')
        # if there are specific requested labels, add them
        if not (dataRangeNames is None) :
            
            # if we don't have exactly the right number of range names to label the ranges
            # then label the tick marks
            if not (len(dataRangeNames) is (len(dataRanges) - 1)) :
                cbar.ax.set_yticklabels(dataRangeNames)
            else : # we will want to label the ranges themselves
                cbar.ax.set_yticklabels(dataRangeNames) # todo, this line is temporary
                doLabelRanges = True
    
    # if there are "tag" masks, plot them over the existing map
    if not (tagData is None) :
        
        # pick out the cooridinates of the points we want to plot
        newX = x[tagData]
        newY = y[tagData]
        
        
        # look at how many trouble points we have
        numTroublePoints = newX.size
        hasTrouble = False
        neededHighlighting = False
        
        if numTroublePoints > 0 :
            hasTrouble = True
            # figure out how many bad points there are
            totalNumPoints = x.size # the number of points
            percentBad = (float(numTroublePoints) / float(totalNumPoints)) * 100.0
            LOG.debug('\t\tnumber of trouble points: ' + str(numTroublePoints))
            LOG.debug('\t\tpercent of trouble points: ' + str(percentBad))
            
            # if there are very few points, make them easier to notice
            # by plotting some colored circles underneath them
            if (percentBad < 0.25) :
                neededHighlighting = True
                p = bMap.plot(newX, newY, 'o', color='#993399', markersize=5)
            elif (percentBad < 1.0) :
                neededHighlighting = True
                p = bMap.plot(newX, newY, 'o', color='#993399', markersize=3)
            
            # plot our point on top of the existing figure
            p = bMap.plot(newX, newY, '.', color='#00FF00', markersize=1)
        
        # display the number of trouble points on the report if we were passed a set of tag data
        # I'm not thrilled with this solution for getting it below the labels drawn by the basemap
        # but I don't think there's a better one at the moment given matplotlib's workings
        troublePtString = '\n\nShowing ' + str(numTroublePoints) + ' Trouble Points'
        # if our plot is more complex, add clarification
        if hasTrouble :
            troublePtString = troublePtString + ' in Green'
            if neededHighlighting :
                troublePtString = troublePtString + '\nwith Purple Circles for Visual Clarity'
        plt.xlabel(troublePtString)
    
    # if we still need to label the ranges, do it now that our fake axis won't mess the trouble points up
    if doLabelRanges :
        """ TODO get this working properly
        fakeAx = plt.axes ([0.77, 0.05, 0.2, 0.9], frameon=False)
        fakeAx.xaxis.set_visible(False)
        fakeAx.yaxis.set_visible(False)
        
        testRect = Rectangle((0, 0), 1, 1, fc="r")
        legendKey = fakeAx.legend([testRect], ["r\n\n\n"], mode="expand", ncol=1, borderaxespad=0.)
        """
    
    return figure

# figure out the bounding axes for the display given a set of
# longitude and latitude and possible a mask of invalid values
# that we should ignore in them
def _get_visible_axes(longitudeData, latitudeData, toIgnoreMask) :
    
    # calculate the bounding range for the display
    # this is in the form [longitude min, longitude max, latitude min, latitude max]
    visibleAxes = [delta.min_with_mask(longitudeData, toIgnoreMask),
                   delta.max_with_mask(longitudeData, toIgnoreMask),
                   delta.min_with_mask(latitudeData,  toIgnoreMask),
                   delta.max_with_mask(latitudeData,  toIgnoreMask)]
    
    return visibleAxes

def spectral_diff_plot( mean_diff, std_diff, max_diff, min_diff, acceptable_diff=None, x=None ):
    """plot spectral difference in current axes, wrapped in std
    >>> x = arange(0,9,0.1)
    >>> y1 = sin(x)
    >>> y2 = sin(x+0.05)
    >>> d = y2-y1
    >>> s = std(d)
    >>> spectral_diff_plot( d, s, array([0.1] * len(d)), x=x )
    """
    cla()
    if x is None: x = range(len(mean_diff))
    if acceptable_diff is not None:
        plot(x, acceptable_diff, 'g.', hold=True, alpha=0.2)
        plot(x, -acceptable_diff, 'g.', hold=True, alpha=0.2)
    plot(x, mean_diff+std_diff, 'r', hold=True, alpha=0.5)
    plot(x,min_diff,'b')
    plot(x,max_diff,'b')
    plot(x, mean_diff-std_diff, 'r', hold=True, alpha=0.5)
    plot(x, mean_diff, 'k', hold=True)

def compare_spectra(actual, desired=None, acceptable=None, x=None):
    """ given an actual[spectrum][wnum], desired[spectrum][wnum], plot comparisons of differences
    """
    delta = actual-desired if (desired is not None) else actual
    d_mean = mean(delta,axis=0)
    d_max = delta.max(axis=0)
    d_min = delta.min(axis=0)
    d_std = std(delta,axis=0)
    des_mean = mean(desired,axis=0)
    subplot(211)
    cla()
    if x is None: x = range(len(des_mean))
    # plot(x,des_mean+d_max,'b')
    # plot(x,des_mean+d_min,'b')
    plot(x,des_mean,'k')
    grid()
    title("mean spectrum")
    subplot(212)
    spectral_diff_plot(d_mean, d_std, d_max, d_min, acceptable, x)
    grid()
    title("difference min-max (blue), mean (black), mean +/- std (red)")
    
def plot_and_save_spacial_trouble(longitude, latitude,
                                  spacialTroubleMask, spaciallyInvalidMask,
                                  fileNameDiscriminator, title, fileBaseName, outputPath, makeSmall=False) :
    """
    given information on spatially placed trouble points in A and B, plot only those points in a very obvious way
    on top of a background plot of a's data shown in grayscale, save this plot to the output path given
    if makeSmall is passed as true a smaller version of the image will also be saved
    """
    
    # get the bounding axis and make a basemap
    boundingAxes = _get_visible_axes(longitude, latitude, spaciallyInvalidMask)
    LOG.debug("Visible axes for lon/lat trouble figure  are: " + str(boundingAxes))
    baseMapInstance, boundingAxes = maps.create_basemap(longitude, latitude, boundingAxes, _select_projection(boundingAxes))
    
    # make the figure
    LOG.info("Creating spatial trouble image")
    spatialTroubleFig = _create_mapped_figure(None, latitude, longitude, baseMapInstance, boundingAxes,
                                                    title, spaciallyInvalidMask, None, spacialTroubleMask)
    # save the figure
    LOG.info("Saving spatial trouble image")
    spatialTroubleFig.savefig(outputPath + "/" + fileBaseName + "." + fileNameDiscriminator + ".png", dpi=fullSizeDPI) 
    
    # we may also save a smaller versions of the figure
    if (makeSmall) :
        
        spatialTroubleFig.savefig(outputPath + "/" + fileBaseName + "." + fileNameDiscriminator + ".small.png", dpi=thumbSizeDPI)
    
    # clear the figure
    spatialTroubleFig.clf()
    plt.close(spatialTroubleFig)
    del(spatialTroubleFig)
    
    return
    
def _handle_fig_creation_task(child_figure_function, log_message,
                              fullFigOutputNamePath, smallFigOutputNamePath,
                              shouldMakeSmall, doFork) :
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
        figure = child_figure_function() 
        LOG.info(log_message)
        figure.savefig(fullFigOutputNamePath, dpi=fullSizeDPI)
        if (shouldMakeSmall) :
            
            tempImage = Image.open(fullFigOutputNamePath)
            scaleFactor = float(thumbSizeDPI) / float(fullSizeDPI)
            originalSize = tempImage.size
            newSize = (int(originalSize[0] * scaleFactor), int(originalSize[1] * scaleFactor))
            tempImage = tempImage.resize(newSize, Image.ANTIALIAS)
            tempImage.save(smallFigOutputNamePath)
            
            #figure.savefig(smallFigOutputNamePath, dpi=thumbSizeDPI)

        # get rid of the figure 
        plt.close(figure)
        del(figure)
    
    # if we've reached this point and we did fork,
    # then we're the child process and we should stop now
    if (doFork) :
        sys.exit(0) # the child is done now
    
    # if we didn't fork, return the 0 pid to indicate that
    return pid

def plot_and_save_figure_comparison(aData, bData,
                                    variableRunInfo, 
                                    fileAName, fileBName,
                                    latitudeAData, longitudeAData,
                                    latitudeBData, longitudeBData,
                                    latitudeCommonData, longitudeCommonData,
                                    spaciallyInvalidMaskA,
                                    spaciallyInvalidMaskB,
                                    outputPath, 
                                    makeSmall=False,
                                    shortCircuitComparisons=False,
                                    doFork=False) : 
    """
    given two files, and information on what to compare, make comparison
    figures and save them to the given output graph.
    Four figures will be generated, one for each file, a comparison image
    that represents the amount of difference in that variable between the
    two files, and an image highlighting the trouble spots where the
    difference exceeds epsilon or there are missing or nan values in one
    or both of the files
    
    variableRunInfo is a dictionary in the form
        variableRunInfo = { 'variable_name': variableName,
                            'epsilon': epsilon,
                            'missing_value': missingDataValue,
                            'display_name': displayName   # this entry is optional; the variable_name should be used
                                                          # for descriptive purposes if it is not defined
                            }
    """
    # if we weren't given a variable display name,
    # just use the standard variable name instead
    variableName = variableRunInfo['variable_name']
    variableDisplayName = variableName
    if 'display_name' in variableRunInfo :
        variableDisplayName = variableRunInfo['display_name']
    
    # figure out what missing values we should be using
    missing_value = variableRunInfo['missing_value']
    missing_value_b = missing_value
    if ('missing_value_alt_in_b' in variableRunInfo) :
        missing_value_b = variableRunInfo['missing_value_alt_in_b']
    
    # compare the two data sets to get our difference data and trouble info
    rawDiffData, goodMask, (goodInAMask, goodInBMask), troubleMask, outsideEpsilonMask, \
    (aNotFiniteMask, bNotFiniteMask), (aMissingMask, bMissingMask), \
    (spaciallyInvalidMaskA, spaciallyInvalidMaskB) = delta.diff(aData, bData, variableRunInfo['epsilon'],
                                                                (missing_value, missing_value_b),
                                                                (spaciallyInvalidMaskA, spaciallyInvalidMaskB))
    absDiffData = np.abs(rawDiffData) # we also want to show the distance between our two, not just which one's bigger/smaller
    
    # some more display info, pull it out for convenience
    dataRanges = None
    if ('display_ranges' in variableRunInfo) :
        dataRanges = variableRunInfo['display_ranges']
    dataRangeNames = None
    if ('display_range_names' in variableRunInfo) : 
        dataRangeNames = variableRunInfo['display_range_names']
    dataColors = None
    if ('display_colors' in variableRunInfo) :
        dataColors = variableRunInfo['display_colors'] 
    
    # figure out the bounding axis
    aAxis = _get_visible_axes(longitudeAData, latitudeAData, ~goodInAMask)
    bAxis = _get_visible_axes(longitudeBData, latitudeBData, ~goodInBMask)
    fullAxis = [min(aAxis[0], bAxis[0]), max(aAxis[1], bAxis[1]),
                min(aAxis[2], bAxis[2]), max(aAxis[3], bAxis[3])]
    LOG.debug("Visible axes for file A variable data (" + variableDisplayName + ") are: " + str(aAxis))
    LOG.debug("Visible axes for file B variable data (" + variableDisplayName + ") are: " + str(bAxis))
    LOG.debug("Visible axes shared for both file's variable data (" + variableDisplayName + ") are: " + str(fullAxis))
    
    if (fullAxis[0] is None) or (fullAxis[1] is None) or (fullAxis[2] is None) or (fullAxis[3] is None) :
        LOG.warn("Unable to display figures for variable (" + variableDisplayName + ") because of inability to identify" +
                 " usable bounding longitude and latitude range on the earth. Bounding range that was identified:" + str(fullAxis))
        return # TODO, the figures need to be disabled from the report and possibly a warning on the report?
    
    # create our basemap
    LOG.info('\t\tloading base map data')
    baseMapInstance, fullAxis = maps.create_basemap(longitudeCommonData, latitudeCommonData, fullAxis, _select_projection(fullAxis))
    
    # from this point on, we will be forking to create child processes so we can parallelize our image and
    # report generation
    isParent = True 
    childPids = []
    
    # the original data figures
    
    # the original A data
    LOG.info("\t\tcreating image of " + variableDisplayName + " in file a")
    pid = _handle_fig_creation_task((lambda : _create_mapped_figure(aData, latitudeAData, longitudeAData,
                                                                    baseMapInstance, fullAxis,
                                                                    (variableDisplayName + "\nin File A"),
                                                                    invalidMask=(~goodInAMask),
                                                                    dataRanges=dataRanges,
                                                                    dataRangeNames=dataRangeNames,
                                                                    dataRangeColors=dataColors)),
                                    "\t\tsaving image of " + variableDisplayName + " for file a",
                                    outputPath + "/" + variableName + ".A.png",
                                    outputPath + "/" + variableName + ".A.small.png",
                                    makeSmall, doFork)
    if not (pid is 0) :
        childPids.append(pid)
        LOG.debug ("Started child process (pid: " + str(pid) + ") to create file a image for " + variableDisplayName)
    
    # the original B data
    LOG.info("\t\tcreating image of " + variableDisplayName + " in file b")
    pid = _handle_fig_creation_task((lambda : _create_mapped_figure(bData, latitudeBData, longitudeBData,
                                                                    baseMapInstance, fullAxis,
                                                                    (variableDisplayName + "\nin File B"),
                                                                    invalidMask=(~ goodInBMask),
                                                                    dataRanges=dataRanges,
                                                                    dataRangeNames=dataRangeNames,
                                                                    dataRangeColors=dataColors)),
                                    "\t\tsaving image of " + variableDisplayName + " for file b",
                                    outputPath + "/" + variableName + ".B.png",
                                    outputPath + "/" + variableName + ".B.small.png",
                                    makeSmall, doFork)
    if not (pid is 0) :
        childPids.append(pid)
        LOG.debug ("Started child process (pid: " + str(pid) + ") to create file b image for " + variableDisplayName)
    
    # make the data comparison figures
    if not shortCircuitComparisons :
        
        # the distance between the two data sets
        LOG.info("\t\tcreating image of the absolute value of difference in " + variableDisplayName)
        pid = _handle_fig_creation_task((lambda : _create_mapped_figure(absDiffData,
                                                                        latitudeCommonData, longitudeCommonData,
                                                                        baseMapInstance, fullAxis,
                                                                        ("Absolute value of difference in\n" + variableDisplayName),
                                                                        invalidMask=(~ goodMask))),
                                        "\t\tsaving image of the absolute value of difference for " + variableDisplayName,
                                        outputPath + "/" + variableName + ".AbsDiff.png",
                                        outputPath + "/" + variableName + ".AbsDiff.small.png",
                                        makeSmall, doFork)
        if not (pid is 0) :
            childPids.append(pid)
            LOG.debug ("Started child process (pid: " + str(pid) + ") to create absolute value of difference image for " + variableDisplayName)
        
        # the subtraction of one data set from the other
        LOG.info("\t\tcreating image of the difference in " + variableDisplayName)
        pid = _handle_fig_creation_task((lambda : _create_mapped_figure(rawDiffData, latitudeCommonData, longitudeCommonData,
                                                                        baseMapInstance, fullAxis,
                                                                        ("Value of (Data File B - Data File A) for\n" + variableDisplayName),
                                                                        invalidMask=(~ goodMask))),
                                        "\t\tsaving image of the difference in " + variableDisplayName,
                                        outputPath + "/" + variableName + ".Diff.png",
                                        outputPath + "/" + variableName + ".Diff.small.png",
                                        makeSmall, doFork)
        if not (pid is 0) :
            childPids.append(pid)
            LOG.debug ("Started child process (pid: " + str(pid) + ") to create difference image for " + variableDisplayName)
        
        # mark the trouble points
        LOG.info("\t\tcreating image marking trouble data in " + variableDisplayName)
        # this is not an optimal solution, but we need to have at least somewhat valid data at any mismatched points so
        # that our plot won't be totally destroyed by missing or non-finite data from B
        bDataCopy = bData[:]
        tempMask = goodInAMask & (~goodInBMask) 
        bDataCopy[tempMask] = aData[tempMask]
        pid = _handle_fig_creation_task((lambda : _create_mapped_figure(bDataCopy, latitudeCommonData, longitudeCommonData,
                                                                        baseMapInstance, fullAxis,
                                                                        ("Areas of trouble data in\n" + variableDisplayName),
                                                                        (~(goodInAMask | goodInBMask)),
                                                                        mediumGrayColorMap, troubleMask,
                                                                        dataRanges=dataRanges,
                                                                        dataRangeNames=dataRangeNames)),
                                        "\t\tsaving image marking trouble data in " + variableDisplayName,
                                        outputPath + "/" + variableName + ".Trouble.png",
                                        outputPath + "/" + variableName + ".Trouble.small.png",
                                        makeSmall, doFork)
        if not (pid is 0) :
            childPids.append(pid)
            LOG.debug ("Started child process (pid: " + str(pid) + ") to create trouble image for " + variableDisplayName)
        
        # a histogram of the values of fileA - file B 
        LOG.info("\t\tcreating histogram of the amount of difference in " + variableDisplayName)
        numBinsToUse = 50
        valuesForHist = rawDiffData[goodMask]
        pid = _handle_fig_creation_task((lambda : _create_histogram(valuesForHist, numBinsToUse,
                                                                    ("Difference in\n" + variableDisplayName),
                                                                    ('Value of (Data File B - Data File A) at a Data Point'),
                                                                    ('Number of Data Points with a Given Difference'),
                                                                    True)),
                                        "\t\tsaving histogram of the amount of difference in " + variableDisplayName,
                                        outputPath + "/" + variableName + ".Hist.png",
                                        outputPath + "/" + variableName + ".Hist.small.png",
                                        makeSmall, doFork)
        if not (pid is 0) :
            childPids.append(pid)
            LOG.debug ("Started child process (pid: " + str(pid) + ") to create difference histogram image for " + variableDisplayName)
        
        # scatter plot of file a vs file b values
        LOG.info("\t\tcreating scatter plot of file a values vs file b values for " + variableDisplayName)
        
        pid = _handle_fig_creation_task((lambda : _create_scatter_plot(aData[goodMask], bData[goodMask],
                                                                       "Value in File A vs Value in File B",
                                                                       "File A Value", "File B Value",
                                                                       outsideEpsilonMask[goodMask],
                                                                       variableRunInfo['epsilon'])),
                                        "\t\tsaving scatter plot of file a values vs file b values in " + variableDisplayName,
                                        outputPath + "/" + variableName + ".Scatter.png",
                                        outputPath + "/" + variableName + ".Scatter.small.png",
                                        makeSmall, doFork)
        if not (pid is 0) :
            childPids.append(pid)
            LOG.debug ("Started child process (pid: " + str(pid) + ") to create scatter plot image for " + variableDisplayName)
    
    # now we need to wait for all of our child processes to terminate before returning
    if (isParent) : # just in case
        if len(childPids) > 0 :
            print ("waiting for completion of " + variableDisplayName + " images...")
        for pid in childPids:
            os.waitpid(pid, 0)
        print("... creation and saving of images for " + variableDisplayName + " completed")
    
    return

if __name__=='__main__':
    import doctest
    doctest.testmod()
