#!/usr/bin/env python
# encoding: utf-8
"""
Plotting routines for different types of figures using matplotlib

Created by evas Dec 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

from pylab import *

import matplotlib.cm     as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from   matplotlib.ticker import FormatStrFormatter

import logging
import numpy as np
from   numpy import ma 

import glance.graphics as maps
import glance.delta    as delta
import glance.report   as report
import glance.stats    as statistics

LOG = logging.getLogger(__name__)

# TODO this value is being used to work around a problem with the contourf
# and how it handles range boundaries. Find a better solution if at all possible.
offsetToRange = 0.0000000000000000001

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
MEDIUM_GRAY_COLOR_MAP = colors.LinearSegmentedColormap('mediumGrayColorMap', mediumGrayColorMapData, 256)

# make an all green color map
greenColorMapData = {
    'red'   : ((0.0, 0.00, 0.00),
               (1.0, 0.00, 0.00)),
    'green' : ((0.0, 1.00, 1.00),
               (1.0, 1.00, 1.00)),
    'blue'  : ((0.0, 0.00, 0.00),
               (1.0, 0.00, 0.00))
}
greenColorMap = colors.LinearSegmentedColormap('greenColorMap', greenColorMapData, 256)

# todo, the use of the offset here is covering a problem with
# contourf hiding data exactly at the end of the range and should
# be removed if a better solution can be found
def _make_range(data_a, valid_a_mask, num_intervals, offset_to_range=0.0, data_b=None, valid_b_mask=None) :
    """
    get an array with numbers representing the bounds of a set of ranges
    that covers all the data present in data_a
    (these may be used for plotting the data)
    if an offset is passed, the outtermost range will be expanded by that much
    if the b data is passed, a total range that encompasses both sets of
    data will be used
    """
    minVal = delta.min_with_mask(data_a, valid_a_mask)
    maxVal = delta.max_with_mask(data_a, valid_a_mask)
    
    # if we have a second set of data, include it in the min/max calculations
    if (data_b is not None) :
        minVal = min(delta.min_with_mask(data_b, valid_b_mask), minVal)
        maxVal = max(delta.max_with_mask(data_b, valid_b_mask), maxVal)
    
    minVal = minVal - offset_to_range
    maxVal = maxVal + offset_to_range
    
    return np.linspace(minVal, maxVal, num_intervals)

def _plot_tag_data_simple(tagData) :
    """
    This method will plot tag data listed as true in the
    tagData mask on the current figure. It is assumed that
    the correlation between the mask and the pixel coordinates
    is exact (ie. no translation is needed).
    
    The return will be the number of points plotted or
    -1 if no valid tagData was given.
    """
    
    numMismatchPoints = -1
    
    # if there are "tag" masks, plot them over the existing map
    if not (tagData is None) :
        
        numMismatchPoints = sum(tagData)
        
        # if we have mismatch points, we need to show them
        if numMismatchPoints > 0:
            
            # figure out how many bad points there are
            totalNumPoints = tagData.size # the number of points
            percentBad = (float(numMismatchPoints) / float(totalNumPoints)) * 100.0
            LOG.debug('\t\tnumber  of mismatch points: ' + str(numMismatchPoints))
            LOG.debug('\t\tpercent of mismatch points: ' + str(percentBad))
            
            new_kwargs = {}
            new_kwargs['cmap'] = greenColorMap
            cleanTagData = ma.array(tagData, mask=~tagData)
            p = contourf(cleanTagData, **new_kwargs)
            # TODO, need to incorporate plot for small numbers of pts
        
        # display the number of mismatch points on the report if we were passed a set of tag data
        mismatchPtString = '\n\nShowing ' + str(numMismatchPoints) + ' Mismatch Points'
        # if our plot is more complex, add clarification
        if numMismatchPoints > 0 :
            mismatchPtString = mismatchPtString + ' in Green'
        plt.xlabel(mismatchPtString)
    
    return numMismatchPoints

def _plot_tag_data_mapped(bMap, tagData, x, y, addExplinationLabel=True) :
    """
    This method will plot the tagged data listed as true in the tagData mask
    on the current figure using the given basemap.
    
    A message will also be added below the map describing the number of
    points plotted, unless the addExplinationLabel variable is passed as False.
    
    The return will be the number of points plotted or
    -1 if no valid tagData was given.
    
    numMismatchPoints = _plot_tag_data_mapped(bMap, tagData, x, y)
    """
    
    numMismatchPoints = -1
    
    # if there are "tag" masks, plot them over the existing map
    if (tagData is not None) and (tagData.size > 0) :
        
        # look at how many mismatch points we have
        numMismatchPoints = sum(tagData)
        neededHighlighting = False
        
        if numMismatchPoints > 0 :
            
            # pick out the cooridinates of the points we want to plot
            newX = np.array(x[tagData])
            newY = np.array(y[tagData])
            
            # figure out how many bad points there are
            totalNumPoints = x.size # the number of points
            percentBad = (float(numMismatchPoints) / float(totalNumPoints)) * 100.0
            LOG.debug('\t\tnumber  of mismatch points: ' + str(numMismatchPoints))
            LOG.debug('\t\tpercent of mismatch points: ' + str(percentBad))
            
            # if there are very few points, make them easier to notice
            # by plotting some colored circles underneath them
            if (percentBad < 0.25) or (totalNumPoints < 20) :
                neededHighlighting = True
                p = bMap.plot(newX, newY, 'o', color='#993399', markersize=5)
            elif (percentBad < 1.0) or (totalNumPoints < 200) :
                neededHighlighting = True
                p = bMap.plot(newX, newY, 'o', color='#993399', markersize=3)
            
            # if there are way too many mismatch points, we can't use plot for this
            if (numMismatchPoints > 1000000) :
                new_kwargs = {}
                new_kwargs['cmap'] = greenColorMap
                p = maps.show_x_y_data(x, y, bMap, data=tagData, **new_kwargs)
            else :
                # plot our point on top of the existing figure
                p = bMap.plot(newX, newY, '.', color='#00FF00', markersize=1)
        
        if addExplinationLabel :
            # display the number of mismatch points on the report if we were passed a set of tag data
            # I'm not thrilled with this solution for getting it below the labels drawn by the basemap
            # but I don't think there's a better one at the moment given matplotlib's workings
            mismatchPtString = '\n\nShowing ' + str(numMismatchPoints) + ' Mismatch Points'
            # if our plot is more complex, add clarification
            if numMismatchPoints > 0 :
                mismatchPtString = mismatchPtString + ' in Green'
                if neededHighlighting :
                    mismatchPtString = mismatchPtString + '\nwith Purple Circles for Visual Clarity'
            plt.xlabel(mismatchPtString)
    
    return numMismatchPoints

# build a scatter plot of the x,y points
def create_scatter_plot(dataX, dataY, title, xLabel, yLabel, badMask=None, epsilon=None, units_x=None, units_y=None) :
    """
    build a scatter plot of the data
    if a bad mask is given the points selected by that mask will be plotted in a different color
    if an epsilon is given the lines for +/- epsilon will be drawn on the plot
    
    by default this plot uses blue for data points and red for data marked by the bad mask
    """
    
    return create_complex_scatter_plot ([(dataX, dataY, badMask,
                                          'b', 'r',
                                          'within\nepsilon', 'outside\nepsilon')],
                                        title, xLabel, yLabel, epsilon=epsilon, units_x=units_x, units_y=units_y)

def create_complex_scatter_plot(dataList, title, xLabel, yLabel, epsilon=None, units_x=None, units_y=None) :
    """
    build a scatter plot with multiple data sets in different colors
    the dataList parameter should be in the form:
    [(set1), (set2), ... , (setN)]
    
    where a set looks like:
    (x data, y data, mask of bad points or None, matlab color code for display, matlab color code for 'bad' points, good label, bad label)
    
    if a mask of bad points is given, it will be applyed to both the x and y data
    
    at least one data set must be given or no image will be created.
    """
    
    # make the figure
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    # if we have no data, stop now
    if (dataList is None) or (len(dataList) <= 0) :
        return figure;
    
    # look at the stuff in each of the data sets and plot that set
    for dataX, dataY, badMask, goodColor, badColor, goodLabel, badLabel in dataList :
        
        # if we have "bad" data to plot, pull it out
        badX = None
        badY = None
        if (badMask != None) :
            badX  = dataX[badMask]
            badY  = dataY[badMask]
            dataX = dataX[~badMask]
            dataY = dataY[~badMask]
        
        # the scatter plot of the good data
        axes.plot(dataX, dataY, ',', color=goodColor, label=goodLabel)
        
        # plot the bad data
        numMismatchPts = 0
        if (badX is not None) and (badY is not None) and (badMask is not None) :
            numMismatchPts = badX.size
            LOG.debug('\t\tplotting ' + str(numMismatchPts) + ' mismatch points in scatter plot.' )
            if numMismatchPts > 0 :
                axes.plot(badX, badY, ',', color=badColor, label=badLabel)
    
    # draw some extra informational lines
    _draw_x_equals_y_line(axes, epsilon=epsilon)
    
    # make a key to explain our plot
    # as long as things have been plotted with proper labels they should show up here
    axes.legend(loc=0, markerscale=3.0) # Note: at the moment markerscale doesn't seem to work
    
    # add the units to the x and y labels
    tempXLabel = xLabel
    tempYLabel = yLabel
    if str.lower(str(units_x)) != "none" :
        tempXLabel = tempXLabel + " in " + units_x
    if str.lower(str(units_y)) != "none" :
        tempYLabel = tempYLabel + " in " + units_y
    
    # and some informational stuff
    axes.set_title(title)
    plt.xlabel(tempXLabel)
    plt.ylabel(tempYLabel)
    
    # format our axes so they display gracefully
    yFormatter = FormatStrFormatter("%4.4g")
    axes.yaxis.set_major_formatter(yFormatter)
    xFormatter = FormatStrFormatter("%4.4g")
    axes.xaxis.set_major_formatter(xFormatter)
    
    return figure

# build a hexbin plot of the x,y points and show the density of the point distribution
def create_hexbin_plot(dataX, dataY, title, xLabel, yLabel, epsilon=None, units_x=None, units_y=None) :
    
    # make the figure
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    # for some reason, if you give the hexplot a data set that's all the same number it dies horribly
    if ( ((dataX is None) or (len(dataX) <= 0)) or ((dataY is None) or (len(dataY) <= 0)) or
         ((dataX.max() == dataX.min()) and (dataY.max() == dataY.min())) ):
        return figure
    
    # the hexbin plot of the good data 
    plt.hexbin(dataX, dataY, bins='log', cmap=cm.jet)
    plt.axis([dataX.min(), dataX.max(), dataY.min(), dataY.max()])
    
    # create a color bar
    cb = plt.colorbar()
    cb.set_label('log10 (count + 1)')
    
    # draw some extra informational lines
    _draw_x_equals_y_line(axes, color='w', epsilon=epsilon, epsilonColor='k')
    
    # add the units to the x and y labels
    tempXLabel = xLabel
    tempYLabel = yLabel
    if str.lower(str(units_x)) != "none" :
        tempXLabel = tempXLabel + " in " + units_x
    if str.lower(str(units_y)) != "none" :
        tempYLabel = tempYLabel + " in " + units_y
    
    # and some informational stuff
    axes.set_title(title)
    plt.xlabel(tempXLabel)
    plt.ylabel(tempYLabel)
    
    # format our axes so they display gracefully
    yFormatter = FormatStrFormatter("%4.4g")
    axes.yaxis.set_major_formatter(yFormatter)
    xFormatter = FormatStrFormatter("%4.4g")
    axes.xaxis.set_major_formatter(xFormatter)
    
    return figure

def _draw_x_equals_y_line(axes, color='k', style='--', epsilon=None, epsilonColor='#00FF00', epsilonStyle='--') :
    """
    Draw the x = y line using the axes and color/style given
    If epsilon is not None, also draw the +/- epsilon lines,
    if they fall in the graph
    """
    
    # get the bounds for our calculations and so we can reset the viewing window later
    xbounds = axes.get_xbound()
    ybounds = axes.get_ybound()
    
    # figure out the size of the ranges
    xrange = xbounds[1] - xbounds[0]
    yrange = ybounds[1] - ybounds[0]
    
    # draw the x=y line
    perfect = [max(xbounds[0], ybounds[0]), min(xbounds[1], ybounds[1])]
    axes.plot(perfect, perfect, style, color=color, label='A = B')
    
    # now draw the epsilon bound lines if they are visible and the lines won't be the same as A = B
    if (not (epsilon is None)) and (epsilon > 0.0) and (epsilon < xrange) and (epsilon < yrange):
        # plot the top line
        axes.plot([perfect[0], perfect[1] - epsilon], [perfect[0] + epsilon, perfect[1]], epsilonStyle, color=epsilonColor, label='+/-epsilon')
        # plot the bottom line
        axes.plot([perfect[0] + epsilon, perfect[1]], [perfect[0], perfect[1] - epsilon], epsilonStyle, color=epsilonColor)
    
    # reset the bounds
    axes.set_xbound(xbounds)
    axes.set_ybound(ybounds)

# build a histogram figure of the given data with the given title and number of bins
def create_histogram(data, bins, title, xLabel, yLabel, displayStats=False, units=None, rangeList=None) :
    
    # make the figure
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    if (data is None) or (len(data) <= 0) :
        return figure
    
    if rangeList is not None :
        assert len(rangeList) == 2
        assert rangeList[0] < rangeList[1]
    
    # the histogram of the data
    n, outBins, patches = plt.hist(data, bins, range=rangeList) # if rangeList is None the range won't be restricted
    
    # format our axes so they display gracefully
    yFormatter = FormatStrFormatter("%3.3g")
    axes.yaxis.set_major_formatter(yFormatter)
    xFormatter = FormatStrFormatter("%.4g")
    axes.xaxis.set_major_formatter(xFormatter)
    
    # add the units to the x and y labels
    tempXLabel = xLabel
    if str.lower(str(units)) != "none" :
        tempXLabel = tempXLabel + " in " + units
    
    # and some informational stuff
    axes.set_title(title)
    plt.xlabel(tempXLabel)
    plt.ylabel(yLabel)
    
    # if stats were passed in, put some of the information on the graph
    # the location is in the form x, y (I think)
    if displayStats :
        # info on the basic stats
        tempMask  = ones(data.shape, dtype=bool)
        tempStats = statistics.NumericalComparisonStatistics.basic_analysis(data, tempMask)
        medianVal = tempStats['median_delta']
        meanVal   = tempStats['mean_delta']
        stdVal    = tempStats['std_val']
        numPts    = data.size
        
        # info on the display of our statistics
        xbounds = axes.get_xbound()
        numBinsToUse = bins
        xrange = xbounds[1] - xbounds[0]
        binSize = xrange / float(numBinsToUse)
        
        # build the display string
        statText = ('%1.2e' % numPts) + ' data points'
        statText = statText + '\n' + 'mean: ' + report.make_formatted_display_string(meanVal)
        statText = statText + '\n' + 'median: ' + report.make_formatted_display_string(medianVal)
        statText = statText + '\n' + 'std: ' + report.make_formatted_display_string(stdVal)
        statText = statText + '\n\n' + 'intervals: ' + report.make_formatted_display_string(numBinsToUse)
        statText = statText + '\n' + 'interval size ' + report.make_formatted_display_string(binSize)
        
        # figure out where to place the text and put it on the figure
        centerOfDisplay = xbounds[0] + (float(xrange) / 2.0)
        xValToUse = 0.67
        # if most of the values will be on the right, move our text to the left...
        if (medianVal > centerOfDisplay) :
            xValToUse = 0.17
        figtext(xValToUse, 0.60, statText)
    
    # make sure we didn't mess up the range if it's being restricted
    # (this may be unnessicary, but it's a good extra precaution)
    if rangeList is not None:
        plt.xlim(rangeList)
    
    return figure

# create a figure including our data mapped onto a map at the lon/lat given
# the colorMap parameter can be used to control the colors the figure is drawn in
# if any masks are passed in the tagData list they will be plotted as an overlays
# set on the existing image
def create_mapped_figure(data, latitude, longitude, baseMapInstance, boundingAxes, title,
                          invalidMask=None, colorMap=None, tagData=None,
                          dataRanges=None, dataRangeNames=None, dataRangeColors=None, units=None, **kwargs) :
    
    # make a clean version of our lon/lat
    latitudeClean  = ma.array(latitude,  mask=~invalidMask)
    longitudeClean = ma.array(longitude, mask=~invalidMask)
    
    # build the plot
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    # build extra info to go to the map plotting function
    kwargs = { }
    
    # figure the range for the color bars
    # this is controllable with the "dataRanges" parameter for discrete data display
    if not (data is None) :
        if dataRanges is None :
            dataRanges = _make_range(data, ~invalidMask, 50, offset_to_range=offsetToRange)
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
    bMap, x, y = maps.show_lon_lat_data(longitudeClean, latitudeClean, baseMapInstance, data=data, **kwargs)
    
    # and some informational stuff
    axes.set_title(title)
    # show a generic color bar
    doLabelRanges = False
    if data is not None :
        cbar = colorbar(format='%.3g')
        # if there are specific requested labels, add them
        if not (dataRangeNames is None) :
            
            # if we don't have exactly the right number of range names to label the ranges
            # then label the tick marks
            if not (len(dataRangeNames) is (len(dataRanges) - 1)) :
                cbar.ax.set_yticklabels(dataRangeNames)
            else : # we will want to label the ranges themselves
                cbar.ax.set_yticklabels(dataRangeNames) # todo, this line is temporary
                doLabelRanges = True
        else :
            # add the units to the colorbar
            if str.lower(str(units)) != "none" :
                cbar.set_label(units)
    
    numMismatchPoints = _plot_tag_data_mapped(bMap, tagData, x, y)
    
    LOG.debug ('number of mismatch points: ' + str(numMismatchPoints))
    
    # if we still need to label the ranges, do it now that our fake axis won't mess the mismatch points up
    if doLabelRanges :
        """ TODO get this working properly
        fakeAx = plt.axes ([0.77, 0.05, 0.2, 0.9], frameon=False)
        fakeAx.xaxis.set_visible(False)
        fakeAx.yaxis.set_visible(False)
        
        testRect = Rectangle((0, 0), 1, 1, fc="r")
        legendKey = fakeAx.legend([testRect], ["r\n\n\n"], mode="expand", ncol=1, borderaxespad=0.)
        """
    
    return figure

# create a figure including a quiver plot of our vector data mapped onto a map at the lon/lat
# given, the colorMap parameter can be used to control the colors the figure is drawn.
# if any masks are passed in the tagData list they will be plotted as an overlays
# set on the existing image
# TODO, this method has not been throughly tested
# TODO, this method needs an input colormap so the mismatch plot can be the right color
def create_quiver_mapped_figure(data, latitude, longitude, baseMapInstance, boundingAxes, title,
                          invalidMask=None, tagData=None, uData=None, vData=None, units=None,  **kwargs) :
    
    # make a clean version of our lon/lat/data
    latitudeClean  =  latitude[~invalidMask]
    longitudeClean = longitude[~invalidMask]
    colorData      = None
    if (data is not None) :
        colorData  =      data[~invalidMask]
    uDataClean     = None
    vDataClean     = None
    if (uData is not None) and (vData is not None) :
        uDataClean =     uData[~invalidMask]
        vDataClean =     vData[~invalidMask]
    tagDataClean   = None
    if tagData is not None :
        tagDataClean = tagData[~invalidMask]
    
    # build the plot
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    # draw our data placed on a map
    maps.draw_basic_features(baseMapInstance, boundingAxes)
    bMap, x, y = maps.show_quiver_plot (longitudeClean, latitudeClean, baseMapInstance, (uDataClean, vDataClean), colordata=colorData)
    
    # show the title
    axes.set_title(title)
    
    # make a color bar if we have color data
    """ todo, currently this crashes glance because of an error in matplot lib
    if colorData is not None :
        cbar = plt.colorbar(format='%.3g')
        # add the units to the colorbar
        if str.lower(str(units)) != "none" :
            cbar.set_label(units)
    """
    
    numMismatchPoints = _plot_tag_data_mapped(bMap, tagDataClean, x, y)
    
    return figure

def create_simple_figure(data, figureTitle, invalidMask=None, tagData=None, colorMap=None, colorbarLimits=None, units=None) :
    """
    create a simple figure showing the data given masked by the invalid mask
    any tagData passed in will be interpreted as mismatch points on the image and plotted as a
    filled contour overlay in green on the image
    if a colorMap is given it will be used to plot the data,
    if not the default colorMap for imshow will be used
    """
    
    cleanData = ma.array(data, mask=invalidMask)
    
    # build the plot
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    # build extra info to go to the map plotting function
    kwargs = { } 
    
    # if we've got a color map, pass it to the list of things we want to tell the plotting function
    if not (colorMap is None) :
        kwargs['cmap'] = colorMap
    
    if (data is not None) and (np.sum(invalidMask) < invalidMask.size) :
        # draw our data
        im = imshow(cleanData, **kwargs)
        
        # if our colorbar has limits set those
        if colorbarLimits is not None :
            LOG.debug("setting colorbar limits: " + str(colorbarLimits))
            clim(vmin=colorbarLimits[0], vmax=colorbarLimits[-1])
        # make a color bar
        cbar = colorbar(format='%.3g')
        # add the units to the colorbar
        if str.lower(str(units)) != "none" :
            cbar.set_label(units)
    
    # and some informational stuff
    axes.set_title(figureTitle)
    
    numMismatchPoints = _plot_tag_data_simple(tagData)
    
    return figure

def create_line_plot_figure(dataList, figureTitle) :
    """
    create a basic line plot of the data vs. it's index, ignoring any invalid data
    if tagData is given, under-label those points with green circles
    
    Each entry in the dataList should be a tupple containing:
            (data, invalidMask, colorString, labelName, tagData, units)
    
    The color string describes a color for plotting in matplotlib.
    The label names will be used for the legend, which will be shown if there is
    more than one set of data plotted or if there is tag data plotted. Invalid
    masks, colors, and label names may be given as None, in which case no data
    will be masked and a default label of "data#" (where # is an arbitrary
    unique counter) will be used.
    tagData may also be passed as None if tagging is not desired in the output.
    units describes the units used to measure the data (such as mm or degrees)
    and will be used to label the plot. units may be passed as None.
    """
    
    # build the plot
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    # plot each of the data sets
    dataSetLabelNumber = 1
    minTagPts  = -1
    maxTagPts  = -1
    plottedTagData = False
    for dataSet, invalidMask, colorString, labelName, tagData, units in dataList :
        
        # if we don't have these, set them to defaults
        if invalidMask is None :
            invalidMask = zeros(dataSet.size, dtype=bool)
        if labelName is None :
            labelName = 'data' + str(dataSetLabelNumber)
            dataSetLabelNumber = dataSetLabelNumber + 1
        if colorString is None:
            colorString = ''
        
        if (dataSet is not None) and (sum(invalidMask) < invalidMask.size) :
            
            # if we don't have a real min yet, set it based on the size
            if minTagPts < 0 :
                minTagPts = dataSet.size + 1
            
            indexData = ma.array(range(dataSet.size), mask=invalidMask.ravel())
            cleanData = ma.array(dataSet.ravel(),     mask=invalidMask.ravel())
            
            # plot the tag data and gather information about it
            if tagData is not None :
                
                plottedTagData = True
                numMismatchPoints = sum(tagData)
                LOG.debug('\t\tnumber of mismatch points: ' + str(numMismatchPoints))
                if numMismatchPoints < minTagPts:
                    minTagPts = numMismatchPoints
                if numMismatchPoints > maxTagPts :
                    maxTagPts = numMismatchPoints
                
                # if we have mismatch points, we need to show them
                if numMismatchPoints > 0:
                    
                    cleanTagData = ma.array(dataSet.ravel(), mask=~tagData.ravel() | invalidMask.ravel())
                    axes.plot(indexData, cleanTagData, 'yo', label='mismatch point')
            
            if str.lower(str(units)) !="none" :
                labelName = labelName + " in " + units
            
            axes.plot(indexData, cleanData, '-' + colorString, label=labelName)
    
    # display the number of mismatch points on the report if we were passed
    # a set of tag data and we were able to compare it to some actual data
    if (plottedTagData and (minTagPts >= 0) and (maxTagPts >=0)) :
        
        mismatchPtString = '\nMarking '
        
        if (minTagPts == maxTagPts) :
            mismatchPtString = mismatchPtString + str(minTagPts) + ' Mismatch Points with Yellow Circles'
        else :
            mismatchPtString = (mismatchPtString + 'between ' + str(minTagPts) + ' and ' + str(maxTagPts) + ' Mismatch Points'
                               + '\non the various data sets (using Yellow Circles)')
        
        plt.xlabel(mismatchPtString)
    
    if (len(dataList) > 1) or plottedTagData :
        # make a key to explain our plot
        # as long as things have been plotted with proper labels they should show up here
        axes.legend(loc=0, markerscale=3.0) # Note: at the moment markerscale doesn't seem to work 
        pass
    
    # and some informational stuff
    axes.set_title(figureTitle)
    
    return figure

if __name__=='__main__':
    import doctest
    doctest.testmod()
