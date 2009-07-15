#!/usr/bin/env python
# encoding: utf-8
"""
Plotting routines for difference values using matplotlib

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging
from pylab import *
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.ticker import FormatStrFormatter

import keoni.map.graphics as maps

import glance.delta as delta
import glance.report as report

LOG = logging.getLogger(__name__)

# TODO this value is being used to work around a problem with the contourf
# and how it handles range boundaries. Find a better solution if at all possible.
offsetToRange = 0.00000000000000001

# the value that will denote "bad" longitudes and latitudes
badLonLat = 1.0E30

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
        numPts = len(data.ravel())
        
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

# create a figure including our data mapped onto a map at the lon/lat given
# the colorMap parameter can be used to control the colors the figure is drawn in
# if any masks are passed in in the tagData list they will be plotted as an overlays
# set on the existing image
def _create_mapped_figure(data, latitude, longitude, boundingAxes, title,
                          invalidMask=None, colorMap=None, tagData=None) :
    
    # this is very inefficient, TODO find a better way?
    latitudeCleaned = empty_like(latitude)
    latitudeCleaned[~invalidMask] = latitude[~invalidMask]
    latitudeCleaned[invalidMask] = badLonLat
    longitudeCleaned = empty_like(longitude)
    longitudeCleaned[~invalidMask] = longitude[~invalidMask]
    longitudeCleaned[invalidMask] = badLonLat
    
    # build the plot
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    # build extra info to go to the map plotting function
    kwargs = {}
    
    # figure the range for the color bars
    if not (data is None) :
        # todo, the use off the offset here is covering a problem with
        # contourf hiding data exactly at the end of the range and should
        # be removed if a better solution can be found
        minVal = _min_with_mask(data, invalidMask) - offsetToRange
        maxVal = _max_with_mask(data, invalidMask) + offsetToRange
        rangeForBar = np.linspace(minVal, maxVal, 50) 
        kwargs['levelsToUse'] = rangeForBar
    
    # if we've got a color map, pass it to the list of things we want to tell the plotting function
    if not (colorMap is None) :
        kwargs['cmap'] = colorMap
    
    # draw our data placed on a map
    bMap, x, y = maps.mapshow(longitudeCleaned, latitudeCleaned, data, boundingAxes, **kwargs)
    
    # and some informational stuff
    axes.set_title(title)
    # show a generic color bar
    if not (data is None) :
        colorbar(format='%.3f') 
    
    # if there are "tag" masks, plot them over the existing map
    if not (tagData is None) :
        
        # pick out the cooridinates of the points we want to plot
        newX = x[tagData]
        newY = y[tagData]
        
        # look at how many trouble points we have
        numTroublePoints = newX.shape[0]
        hasTrouble = False
        neededHighlighting = False
        
        if numTroublePoints > 0 :
            hasTrouble = True
            # figure out how many bad points there are
            totalNumPoints = x[~invalidMask].shape[0]
            percentBad = (float(numTroublePoints) / float(totalNumPoints)) * 100.0
            LOG.debug('\t\tnumber of trouble points: ' + str(numTroublePoints))
            LOG.debug('\t\tpercent of trouble points: ' + str(percentBad))
            
            # if there are very few points, make them easier to notice
            # by plotting some colored circles underneath them
            if (percentBad < 0.25) :
                neededHighlighting = True
                bMap.plot(newX, newY, 'o', color='#993399', markersize=5)
            elif (percentBad < 1.0) :
                neededHighlighting = True
                bMap.plot(newX, newY, 'o', color='#993399', markersize=3)
            
            # plot our point on top of the existing figure
            bMap.plot(newX, newY, '.', color='#00FF00', markersize=1)
            
        
        # display the number of trouble points on the report if we were passed a set of tag data
        # I'm not thrilled with this solution for getting it below the labels drawn by the basemap
        # but I don't think there's a better one at the moment given matplotlib's workings
        troublePtString = '\n\n\nShowing ' + str(numTroublePoints) + ' Trouble Points'
        # if our plot is more complex, add clarification
        if hasTrouble :
            troublePtString = troublePtString + ' in Green'
            if neededHighlighting :
                troublePtString = troublePtString + '\nwith Purple Circles for Visual Clarity'
        plt.xlabel(troublePtString)

    return figure

# get the min, ignoring the stuff in mask
def _min_with_mask(data, mask) :
    return data[~mask].ravel()[data[~mask].argmin()]
    
# get the max, ignoring the stuff in mask
def _max_with_mask(data, mask) :
    return data[~mask].ravel()[data[~mask].argmax()]
    
# figure out the bounding axes for the display given a set of
# longitude and latitude and possible a mask of invalid values
# that we should ignore in them
def _get_visible_axes(longitudeData, latitudeData, toIgnoreMask) :
    
    # calculate the bounding range for the display
    # this is in the form [longitude min, longitude max, latitude min, latitude max]
    visibleAxes = [_min_with_mask(longitudeData, toIgnoreMask),
                   _max_with_mask(longitudeData, toIgnoreMask),
                   _min_with_mask(latitudeData, toIgnoreMask),
                   _max_with_mask(latitudeData, toIgnoreMask)]
    
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
    # get bounding axes
    visibleAxes = _get_visible_axes(longitude, latitude, spaciallyInvalidMask)
    
    # make the figure
    spatialTroubleFig = _create_mapped_figure(None, latitude, longitude, visibleAxes, title,
                                              spaciallyInvalidMask, None, spacialTroubleMask)
    # save the figure
    spatialTroubleFig.savefig(outputPath + "/" + fileBaseName + "." + fileNameDiscriminator + ".png", dpi=200) 
    
    # we may also save a smaller versions of the figure
    if (makeSmall) :
        
        spatialTroubleFig.savefig(outputPath + "/" + fileBaseName + "." + fileNameDiscriminator + ".small.png", dpi=50)
    
    return
    
def plot_and_save_figure_comparison(aData, bData, variableName,
                                    fileAName, fileBName,
                                    latitudeData, longitudeData,
                                    spaciallyInvalidMaskA,
                                    spaciallyInvalidMaskB,
                                    spaciallyInvalidMaskBoth,
                                    outputPath, epsilon, missing,
                                    makeSmall=False, variableDisplayName=None) :
    """
    given two files, and information on what to compare, make comparison
    figures and save them to the given output graph.
    Four figures will be generated, one for each file, a comparison image
    that represents the amount of difference in that variable between the
    two files, and an image highlighting the trouble spots where the
    difference exceeds epsilon or there are missing or nan values in one
    or both of the files
    """
    # if we weren't given a variable display name,
    # just use the standard variable name instead
    if variableDisplayName is None :
        variableDisplayName = variableName
    
    print("Creating figures for: " + variableDisplayName)
    
    # compare the two data sets to get our difference data and trouble info
    rawDiffData, goodMask, troubleMask, (aNotFiniteMask, bNotFiniteMask), \
    (aMissingMask, bMissingMask) = delta.diff(aData, bData, epsilon, (missing,missing), spaciallyInvalidMaskBoth)
    diffData = np.abs(rawDiffData) # we want to show the distance between our two, rather than which one's bigger
    
    # mark where our invalid data is for each of the files (a and b) 
    invalidDataMaskA = aMissingMask | aNotFiniteMask
    invalidDataMaskB = bMissingMask | bNotFiniteMask
    # this mask potentially represents data we don't want to try to plot in our diff because it may be malformed
    everyThingWrongButEpsilon = spaciallyInvalidMaskBoth | invalidDataMaskA | invalidDataMaskB
    # use an exclusive or to get a mask for just the points deemed "bad" by epsilon comparison
    tooDifferentMask = everyThingWrongButEpsilon ^ troubleMask
    
    # calculate the bounding range for the display
    # this is in the form [longitude min, longitude max, latitude min, latitude max]
    visibleAxesA = _get_visible_axes(longitudeData, latitudeData, spaciallyInvalidMaskA)
    visibleAxesB = _get_visible_axes(longitudeData, latitudeData, spaciallyInvalidMaskB)
    visibleAxesBoth = _get_visible_axes(longitudeData, latitudeData, spaciallyInvalidMaskBoth)
    
    # make the three figures
    print("\tcreating image of file a")
    figureA = _create_mapped_figure(aData, latitudeData, longitudeData, visibleAxesA,
                                    (variableDisplayName + "\nin File A"),
                                    invalidMask=(spaciallyInvalidMaskA | invalidDataMaskA))
    print("\tcreating image of file b")
    figureB = _create_mapped_figure(bData, latitudeData, longitudeData, visibleAxesB,
                                    (variableDisplayName + "\nin File B"),
                                    invalidMask=(spaciallyInvalidMaskB | invalidDataMaskB))
    print("\tcreating image of the absolue value of difference")
    figureAbsDiff = _create_mapped_figure(diffData, latitudeData, longitudeData, visibleAxesBoth, 
                                          ("Absolute value of difference in\n" + variableDisplayName),
                                          invalidMask=(everyThingWrongButEpsilon))
    print("\tcreating image of the difference")
    figureDiff = _create_mapped_figure(rawDiffData, latitudeData, longitudeData, visibleAxesBoth, 
                                          ("Value of (Data File B - Data File A) for\n" + variableDisplayName),
                                          invalidMask=(everyThingWrongButEpsilon))
    # this figure is more complex because we want to mark the trouble points on it
    print("\tcreating image marking trouble data")
    figureBadDataInDiff = _create_mapped_figure(bData, latitudeData, longitudeData, visibleAxesBoth,
                                                ("Areas of trouble data in\n" + variableDisplayName),
                                                spaciallyInvalidMaskBoth | invalidDataMaskB,
                                                mediumGrayColorMap, troubleMask)
    # a histogram of the values of fileA - file B so that the distribution of error is visible (hopefully)
    print("\tcreating histogram of the amount of difference")
    numBinsToUse = 50
    diffHistogramFigure = _create_histogram(rawDiffData[~everyThingWrongButEpsilon].ravel(), numBinsToUse,
                                            ("Difference in\n" + variableDisplayName),
                                            ('Value of (Data File B - Data File A) at a Data Point'),
                                            ('Number of Data Points with a Given Difference'),
                                            True)
    
    # TEMP scatter plot test of file a and b comparison
    print("\tcreating scatter plot of file a values vs file b values")
    diffScatterPlot = _create_scatter_plot(aData[~everyThingWrongButEpsilon].ravel(), bData[~everyThingWrongButEpsilon].ravel(),
                                           "Value in File A vs Value in File B", "File A Value", "File B Value",
                                           tooDifferentMask[~everyThingWrongButEpsilon].ravel(), epsilon)
    
    # save the figures to disk
    print("Saving figures for: " + variableDisplayName)
    print("\tsaving image of file a")
    figureA.savefig(outputPath + "/" + variableName + ".A.png", dpi=200) 
    print("\tsaving image of file b")
    figureB.savefig(outputPath + "/" + variableName + ".B.png", dpi=200) 
    print("\tsaving image of the absolute value of difference")
    figureAbsDiff.savefig(outputPath + "/" + variableName + ".AbsDiff.png", dpi=200)
    print("\tsaving image of the difference")
    figureDiff.savefig(outputPath + "/" + variableName + ".Diff.png", dpi=200) 
    print("\tsaving image marking trouble data")
    figureBadDataInDiff.savefig(outputPath + "/" + variableName + ".Trouble.png", dpi=200) 
    print("\tsaving histogram of the amount of difference")
    diffHistogramFigure.savefig(outputPath + "/" + variableName + ".Hist.png", dpi=200) 
    print("\tsaving scatter plot of file a values vs file b values")
    diffScatterPlot.savefig(outputPath + "/" + variableName + ".Scatter.png", dpi=200) 
    
    # also save smaller versions of the figures if the parameter says the caller wants us to
    if (makeSmall) :
        print("\tsaving smaller versions of images")
        figureA.savefig(outputPath + "/" + variableName + ".A.small.png", dpi=50)
        figureB.savefig(outputPath + "/" + variableName + ".B.small.png", dpi=50)
        figureAbsDiff.savefig(outputPath + "/" + variableName + ".AbsDiff.small.png", dpi=50)
        figureDiff.savefig(outputPath + "/" + variableName + ".Diff.small.png", dpi=50)
        figureBadDataInDiff.savefig(outputPath + "/" + variableName + ".Trouble.small.png", dpi=50)
        diffHistogramFigure.savefig(outputPath + "/" + variableName + ".Hist.small.png", dpi=50)
        diffScatterPlot.savefig(outputPath + "/" + variableName + ".Scatter.small.png", dpi=50)
    
    return

if __name__=='__main__':
    import doctest
    doctest.testmod()
