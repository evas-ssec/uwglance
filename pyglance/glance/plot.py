#!/usr/bin/env python
# encoding: utf-8
"""
Plotting routines for difference values using matplotlib

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging
from pylab import *

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
import keoni.map.graphics as maps

import glance.delta as delta

LOG = logging.getLogger(__name__)

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
    
def plot_and_save_figure_comparison(aData, bData, variableName,
                                    fileAName, fileBName,
                                    latitudeData, longitudeData, outputPath,
                                    epsilon, missing) :
    """
    given two files, and information on what to compare, make comparison
    figures and save them to the given output graph.
    Four figures will be generated, one for each file, a comparison image
    that represents the amount of difference in that variable between the
    two files, and an image highlighting the trouble spots where the
    difference exceeds epsilon or there are missing or nan values in one
    or both of the files
    """
    print("Creating figures for: " + variableName)
    
    # build a mask of our spacially invalid data so we can ask our
    # comparison tool to ignore it
    invalidLatitude = (latitudeData < -90) | (latitudeData > 90)
    invalidLongitude = (longitudeData < -180)   | (longitudeData > 180)
    spaciallyInvalidMask = invalidLatitude | invalidLongitude
    
    # compare the two data sets to get our difference data and trouble info
    rawDiffData, goodMask, troubleMask, (aNotFiniteMask, bNotFiniteMask), \
    (aMissingMask, bMissingMask) = delta.diff(aData, bData, epsilon, (missing,missing), spaciallyInvalidMask)
    diffData = abs(rawDiffData) # we want to show the distance between our two, rather than which one's bigger
    # mark where our invalid data is for each of the files (a and b) 
    invalidDataMaskA = aMissingMask | aNotFiniteMask
    invalidDataMaskB = bMissingMask | bNotFiniteMask
    # this mask potentially represents data we don't want to try to plot in our diff because it may be malformed
    everyThingWrongButEpsilon = spaciallyInvalidMask | invalidDataMaskA | invalidDataMaskB
    # use an exclusive or to get a mask for just the points deemed "bad" by epsilon comparison
    tooDifferentMask = everyThingWrongButEpsilon ^ troubleMask

    # calculate the bounding range for the display
    # this is in the form [longitude min, longitude max, latitude min, latitude max]
    visibleAxes = [_min_with_mask(longitudeData, spaciallyInvalidMask),
                   _max_with_mask(longitudeData, spaciallyInvalidMask),
                   _min_with_mask(latitudeData, spaciallyInvalidMask),
                   _max_with_mask(latitudeData, spaciallyInvalidMask)]
    
    # make the three figures
    print("\tcreating image of file a")
    figureA = _create_mapped_figure(aData, latitudeData, longitudeData, visibleAxes,
                                    (variableName + "\nin File A"),
                                    invalidMask=(spaciallyInvalidMask | invalidDataMaskA))
    print("\tcreating image of file b")
    figureB = _create_mapped_figure(bData, latitudeData, longitudeData, visibleAxes,
                                    (variableName + "\nin File B"),
                                    invalidMask=(spaciallyInvalidMask | invalidDataMaskB))
    print("\tcreating image of the amount of difference")
    figureDiff = _create_mapped_figure(diffData, latitudeData, longitudeData, visibleAxes,
                                       ("Amount of difference in\n" + variableName),
                                       invalidMask=(everyThingWrongButEpsilon))
    # this figure is more complex because we want to mark the trouble points on it
    print("\tcreating image marking trouble data")
    figureBadDataInDiff = _create_mapped_figure(bData, latitudeData, longitudeData, visibleAxes,
                                                ("Areas of trouble data in\n" + variableName),
                                                spaciallyInvalidMask | invalidDataMaskB,
                                                mediumGrayColorMap, troubleMask)
    # a histogram of the values of fileA - file B so that the distribution of error is visible (hopefully)
    print("\tcreating histogram of the amount of difference")
    numBinsToUse = 50
    diffHistogramFigure = _create_histogram(rawDiffData[~everyThingWrongButEpsilon].ravel(), numBinsToUse,
                                            ("Difference in\n" + variableName),
                                            ('Value of (Data File B - Data File A) at a Data Point (in ' + str(numBinsToUse) + ' bins)'),
                                            ('Number of Data Points with a Given Difference'))
    
    # TEMP scatter plot test of file a and b comparison
    print("\tcreating scatter plot of file a values vs file b values")
    diffScatterPlot = _create_scatter_plot(aData[~everyThingWrongButEpsilon].ravel(), bData[~everyThingWrongButEpsilon].ravel(),
                                           "Value in File A vs Value in File B", "File A Value", "File B Value",
                                           tooDifferentMask[~everyThingWrongButEpsilon].ravel())
    
    # save the figures to disk
    print("Saving figures for: " + variableName)
    print("\tsaving image of file a")
    figureA.savefig(outputPath + "/" + variableName + ".A.png", dpi=200)
    print("\tsaving image of file b")
    figureB.savefig(outputPath + "/" + variableName + ".B.png", dpi=200)
    print("\tsaving image of the amount of difference")
    figureDiff.savefig(outputPath + "/" + variableName + ".Diff.png", dpi=200)
    print("\tsaving image marking trouble data")
    figureBadDataInDiff.savefig(outputPath + "/" + variableName + ".Trouble.png", dpi=200)
    print("\tsaving histogram of the amount of difference")
    diffHistogramFigure.savefig(outputPath + "/" + variableName + ".Hist.png", dpi=200)
    print("\tsaving scatter plot of file a values vs file b values")
    diffScatterPlot.savefig(outputPath + "/" + variableName + ".Scatter.png", dpi=200)
    
    return

# build a scatter plot of the x,y points
def _create_scatter_plot(dataX, dataY, title, xLabel, yLabel, badMask=None) :
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
    
    # the scatter plot of the data
    axes.plot(dataX, dataY, 'b,')
    if (badMask != None) :
        LOG.debug('\t\tnumber of trouble points: ' + str(badX.shape))
        axes.plot(badX, badY, 'r,')
    
    # and some informational stuff
    axes.set_title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    
    return figure

# build a histogram figure of the given data with the given title and number of bins
def _create_histogram(data, bins, title, xLabel, yLabel) :
    
    # make the figure
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    # the histogram of the data
    n, bins, patches = plt.hist(data, bins)
    
    # and some informational stuff
    axes.set_title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    
    return figure
    

# create a figure including our data mapped onto a map at the lon/lat given
# the colorMap parameter can be used to control the colors the figure is drawn in
# if tagData is passed in it will be plotted as an overlayed set on the existing
# image using the tagDataColorMap colors rather than those of the original image
def _create_mapped_figure(data, latitude, longitude, boundingAxes, title,
                          invalidMask=None, colorMap=None, tagDataToShowMask=None) :
    
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
    
    # figure the range for the color bars
    minVal = _min_with_mask(data, invalidMask)
    maxVal = _max_with_mask(data, invalidMask)
    rangeForBar = arange(minVal, maxVal, (maxVal - minVal) / 50)
    
    # draw our data placed on a map
    if (colorMap != None) :
        bMap, x, y = maps.mapshow(longitudeCleaned, latitudeCleaned, data, boundingAxes, levelsToUse=rangeForBar, cmap=colorMap)
    else :
        bMap, x, y = maps.mapshow(longitudeCleaned, latitudeCleaned, data, boundingAxes, levelsToUse=rangeForBar)
    
    # and some informational stuff
    axes.set_title(title)
    # show a generic color bar
    colorbar() 
    
    # if there's a"tag" mask, plot it over the existing map
    if (tagDataToShowMask != None) :
        
        # pick out the cooridinates of the points we want to plot
        newX = x[tagDataToShowMask]
        newY = y[tagDataToShowMask]
        
        # look at how many trouble points we have
        numTroublePoints = newX.shape[0]
        LOG.debug('\t\tnumber of trouble points: ' + str(numTroublePoints))
        # figure out how big to make the points when we plot them
        totalNumPoints = x[~invalidMask].shape[0]
        pointSize = 1
        percentBad = (float(numTroublePoints) / float(totalNumPoints)) * 100.0
        # if we have very few bad data points, make them bigger so they're easier to see
        LOG.debug('\t\tpercent of trouble points: ' + str(percentBad))
        if (percentBad > 1.0) :
            pointSize = 1
        elif (percentBad > 0.25) :
            pointSize = 3
        else :
            pointSize = 5
        
        # plot our point on top of the existing figure
        bMap.plot(newX, newY, '.', color='#00FF00', markersize=pointSize)

    return figure

# get the min, ignoring the stuff in mask
def _min_with_mask(data, mask) :
    return data[~mask].ravel()[data[~mask].argmin()]
    
# get the max, ignoring the stuff in mask
def _max_with_mask(data, mask) :
    return data[~mask].ravel()[data[~mask].argmax()]

if __name__=='__main__':
    import doctest
    doctest.testmod()
