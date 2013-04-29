#!/usr/bin/env python
# encoding: utf-8
"""
This module manages creating figures for the Glance GUI.

Created by evas Nov 2011.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
"""

# these first two lines must stay before the pylab import
import matplotlib
# Note: it's assumed that you've already set up this use previously
#matplotlib.use('Qt4Agg') # use the Qt Anti-Grain Geometry rendering engine

from pylab import *

import matplotlib.cm     as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import logging
import numpy as np

import glance.data      as dataobjects
import glance.figures   as figures
import glance.gui_model as model
from   glance.gui_constants import *

LOG = logging.getLogger(__name__)

#temp_dict = {'blue': [(0.0, 0.75, 0.75), (0.11, 0.99955436720142599, 0.99955436720142599), (0.34000000000000002, 0.99810246679316883, 0.99810246679316883), (0.34999999999999998, 0.98545224541429477, 0.98545224541429477), (0.375, 0.94117647058823528, 0.94117647058823528), (0.64000000000000001, 0.51739405439595187, 0.51739405439595187), (0.65000000000000002, 0.5, 0.5), (0.66000000000000003, 0.5, 0.5), (0.89000000000000001, 0.5, 0.5), (0.91000000000000003, 0.5, 0.5), (1.0, 0.5, 0.5)], 'green': [(0.0, 0.5, 0.5), (0.11, 0.5, 0.5), (0.125, 0.50098039215686274, 0.50098039215686274), (0.34000000000000002, 0.93235294117647061, 0.93235294117647061), (0.34999999999999998, 0.94803921568627447, 0.94803921568627447), (0.375, 1.0, 1.0), (0.64000000000000001, 1.0, 1.0), (0.65000000000000002, 0.97966594045025435, 0.97966594045025435), (0.66000000000000003, 0.96514161220043593, 0.96514161220043593), (0.89000000000000001, 0.53667392883079168, 0.53667392883079168), (0.91000000000000003, 0.50036310820624552, 0.50036310820624552), (1.0, 0.5, 0.5)], 'red': [(0.0, 0.5, 0.5), (0.11, 0.5, 0.5), (0.125, 0.5, 0.5), (0.34000000000000002, 0.5, 0.5), (0.34999999999999998, 0.5, 0.5), (0.375, 0.54269449715370022, 0.54269449715370022), (0.64000000000000001, 0.96647691334598351, 0.96647691334598351), (0.65000000000000002, 0.98545224541429466, 0.98545224541429466), (0.66000000000000003, 0.99810246679316883, 0.99810246679316883), (0.89000000000000001, 0.99955436720142621, 0.99955436720142621), (0.91000000000000003, 0.9549910873440286, 0.9549910873440286), (1.0, 0.75, 0.75)]}
temp_dict = {'blue': [(0.0, 0.58333333333333326, 0.58333333333333326), (0.11, 0.91607248960190135, 0.91607248960190135), (0.125, 0.91666666666666663, 0.91666666666666663), (0.34000000000000002, 0.91413662239089188, 0.91413662239089188), (0.34999999999999998, 0.89726966055239299, 0.89726966055239299), (0.375, 0.83823529411764708, 0.83823529411764708), (0.64000000000000001, 0.27319207252793593, 0.27319207252793593), (0.65000000000000002, 0.25, 0.25), (0.66000000000000003, 0.25, 0.25), (0.89000000000000001, 0.25, 0.25), (0.91000000000000003, 0.25, 0.25), (1.0, 0.25, 0.25)], 'green': [(0.0, 0.25, 0.25), (0.11, 0.25, 0.25), (0.125, 0.25130718954248366, 0.25130718954248366), (0.34000000000000002, 0.82647058823529418, 0.82647058823529418), (0.34999999999999998, 0.84738562091503267, 0.84738562091503267), (0.375, 0.91666666666666663, 0.91666666666666663), (0.64000000000000001, 0.91666666666666663, 0.91666666666666663), (0.65000000000000002, 0.88955458726700576, 0.88955458726700576), (0.66000000000000003, 0.87018881626724787, 0.87018881626724787), (0.89000000000000001, 0.29889857177438889, 0.29889857177438889), (0.91000000000000003, 0.25048414427499405, 0.25048414427499405), (1.0, 0.25, 0.25)], 'red': [(0.0, 0.25, 0.25), (0.11, 0.25, 0.25), (0.125, 0.25, 0.25), (0.34000000000000002, 0.25, 0.25), (0.34999999999999998, 0.25, 0.25), (0.375, 0.30692599620493355, 0.30692599620493355), (0.64000000000000001, 0.87196921779464465, 0.87196921779464465), (0.65000000000000002, 0.89726966055239288, 0.89726966055239288), (0.66000000000000003, 0.91413662239089177, 0.91413662239089177), (0.89000000000000001, 0.91607248960190157, 0.91607248960190157), (0.91000000000000003, 0.85665478312537158, 0.85665478312537158), (1.0, 0.58333333333333326, 0.58333333333333326)]}
DESAT_MAP = matplotlib.colors.LinearSegmentedColormap('colormap', temp_dict, 1024)

# colormaps that are available in the GUI
# TODO, if this changes the list of colormap names in the constants module needs to be kept up
AVAILABLE_COLORMAPS = {CM_RAINBOW:       cm.jet,
                       CM_RAINBOW_REV:   cm.jet_r,
                       CM_RAINBOW_DESAT: DESAT_MAP,
                       CM_GRAY:          cm.bone,
                       CM_GRAY_REV:      cm.bone_r,
                       CM_SPECTRAL:      cm.spectral}

class GlanceGUIFigures (object) :
    """
    This class handles creating figures for the glance gui.
    
    (in future it may manage them more actively)
    
    it includes:
    
    self.dataModel      - the GlanceGUIModel object that contains the main data
                          model for the GUI
    self.errorHandlers  - objects that want to be notified when there's a serious error
    """
    
    def __init__ (self, dataModelToSave) :
        """
        create a figure manager, hanging on to the data model, for use in creating figures
        """
        
        self.dataModel     = dataModelToSave
        self.errorHandlers = [ ]
    
    def registerErrorHandler (self, objectToRegister) :
        """
        add the given object to our list of error handlers
        """
        
        if objectToRegister not in self.errorHandlers :
            self.errorHandlers.append(objectToRegister)
    
    def _getVariableInformation (self, filePrefix, variableName=None) :
        """
        Pull the name, data, and units for the variable currently selected in the given file prefix
        """
        varNameToUse = variableName
        if varNameToUse is None :
            varNameToUse = self.dataModel.getVariableName(filePrefix) # get the currently selected variable
        
        dataObject           = self.dataModel.getVariableData(filePrefix, varNameToUse,  doCorrections=True)
        unitsText            = self.dataModel.getUnitsText   (filePrefix, varNameToUse)
        
        if dataObject is not None :
            dataObject.self_analysis()
        
        return varNameToUse, dataObject, unitsText
    
    def _getVariableInfoSmart (self, filePrefix, imageType) :
        """
        if appropriate for the image type, get information on the variable, otherwise return None's
        """
        
        varName, dataObject, unitsText = None, None, None
        
        # only load the data if it will be needed for the plot
        if ( self.dataModel.getShouldShowOriginalPlotsInSameRange() or
             ((imageType == ORIGINAL_A)  and (filePrefix == "A") or
              (imageType == ORIGINAL_B)  and (filePrefix == "B") or
              (imageType == HISTOGRAM_A) and (filePrefix == "A") or
              (imageType == HISTOGRAM_B) and (filePrefix == "B") or
              (imageType in COMPARISON_IMAGES))) :
            varName, dataObject, unitsText = self._getVariableInformation(filePrefix)
        
        return varName, dataObject, unitsText
    
    def _buildDiffInfoObjectSmart (self, imageType, dataObjectA, dataObjectB, varNameA, varNameB,
                                   epsilon_value=None, epsilon_percent=None) :
        """
        if appropriate for the image type, build the difference object, otherwise return None
        
        this method may rase an IncompatableDataObjects exception if the two data objects it's given can't be compared
        """
        
        diffObject = None
        
        # only build the difference if we need to compare the data
        if imageType in COMPARISON_IMAGES :
            
            # check to see if our data is minimally compatable; this call may raise an IncompatableDataObjects exception
            dataobjects.DiffInfoObject.verifyDataCompatability (dataObjectA, dataObjectB, varNameA, varNameB)
            
            # compare our data
            diffObject = dataobjects.DiffInfoObject(dataObjectA, dataObjectB,
                                                    epsilonValue=epsilon_value, epsilonPercent=epsilon_percent)
        
        return diffObject
    
    def _load_and_analyse_lonlat (self, listOfFilePrefixes=["A", "B"], lonNames=None, latNames=None, stopIfComparisonFails=False) :
        """
        load information on the longidue and latitude,
        if there are multiple file prefixes given:
            find the shared range
            analyse how different the navigation is between the files
            (if there is a lon/lat epsilon defined and the difference is more than that, either stop with an error or log a warning)
        
        lonNames and latNames should be dictionaries giving the names of the longitude and latitude variables indexed by the file prefixes
        
        This method may raise an IncompatableDataObjects exception if multiple file prefixes are passed in the listOfFilePrefixes
        and the longitude and latidues for those files can not be compared.
        """
        
        lonlatData = { }
        lonRange   = None
        latRange   = None
        
        # load and compare stuff for each file prefix
        for filePrefix in listOfFilePrefixes :
            
            # get information on the lon/lat from the current file
            currentLonObj, currentLatObj, currentLonRange, currentLatRange = self._load_lonlat(filePrefix, lonNames[filePrefix], latNames[filePrefix])
            
            # TODO, this will currently crash if there's a problem, we don't really want that
            assert currentLonObj.data.shape == currentLatObj.data.shape
            
            # expand our lon/lat ranges if we need to
            if lonRange is None :
                lonRange = currentLonRange
            else :
                lonRange[0] = min(currentLonRange[0], lonRange[0])
                lonRange[1] = max(currentLonRange[1], lonRange[1])
            if latRange is None:
                latRange = currentLatRange
            else :
                latRange[0] = min(currentLatRange[0], latRange[0])
                latRange[1] = max(currentLatRange[1], latRange[1])
            
            # compare this file to whatever other data we have
            for filePrefixToCompare in lonlatData.keys() :
                lonToCompare, latToCompare = lonlatData[filePrefixToCompare]
                # TODO, this is going to crash if there's a problem, we don't really want that
                assert lonToCompare.data.shape == currentLatObj.data.shape
                assert lonToCompare.data.shape == currentLonObj.data.shape
            
            # add this data to the list of lonlat data
            lonlatData[filePrefix] = [currentLonObj, currentLatObj]
        
        # return longitude and latitude information and the shared ranges
        return lonlatData, lonRange, latRange
    
    def _load_lonlat (self, filePrefix, lonName, latName) :
        """
        load the longitude and latitude information for the file and determine the ranges
        present in both
        """
        
        _, lonObject, _ = self._getVariableInformation(filePrefix, lonName)
        _, latObject, _ = self._getVariableInformation(filePrefix, latName)
        
        lonRange = [lonObject.get_min(), lonObject.get_max()]
        latRange = [latObject.get_min(), latObject.get_max()]
        
        return lonObject, latObject, lonRange, latRange
    
    def spawnPlot (self) :
        """
        create a matplotlib plot using the current model information
        
        this method may raise an IncompatableDataObjects exception if the a and b data
        are completely incomparable
        this method may also raise a ValueError if the data could not be plotted
        for reasons not encompassed by an IncompatableDataObjects exception
        """
        
        # retrieve some plotting settings
        imageType     = self.dataModel.getImageType()
        dataForm      = self.dataModel.getDataForm()
        colorMapToUse = AVAILABLE_COLORMAPS[self.dataModel.getColormapName()]
        
        LOG.info ("Preparing variable data for plotting...")
        
        # load the variable data
        aVarName, aDataObject, aUnitsText = self._getVariableInfoSmart("A", imageType)
        bVarName, bDataObject, bUnitsText = self._getVariableInfoSmart("B", imageType)
        # compare the variables 
        diffData = self._buildDiffInfoObjectSmart(imageType,
                                                  aDataObject, bDataObject,
                                                  aVarName,    bVarName,
                                                  epsilon_value=self.dataModel.getEpsilon(),
                                                  epsilon_percent=self.dataModel.getEpsilonPercent())
        
        # if we need to build a shared range, do that now
        rangeInfo = None
        if (self.dataModel.getShouldShowOriginalPlotsInSameRange() and (aDataObject is not None) and (bDataObject is not None)) :
            rangeInfo = [min(aDataObject.get_min(), bDataObject.get_min()), max(aDataObject.get_max(), bDataObject.get_max())]
        
        # if the user asked for a mapped plotting format and type of plot that is mapped
        if ((dataForm == MAPPED_2D) and (imageType != HISTOGRAM) and
                                        (imageType != HISTOGRAM_A) and
                                        (imageType != HISTOGRAM_B) and 
                                        (imageType != model.SCATTER) and
                                        (imageType != model.HEX_PLOT)) :
            lonNames = {
                        "A": self.dataModel.getLongitudeName("A"),
                        "B": self.dataModel.getLongitudeName("B")
                        }
            latNames = {
                        "A": self.dataModel.getLatitudeName("A"),
                        "B": self.dataModel.getLatitudeName("B")
                        }
            lonlatData, lonRange, latRange = self._load_and_analyse_lonlat(listOfFilePrefixes=["A", "B"],
                                                                           lonNames=lonNames, latNames=latNames)
            
            # double check that lon/lat are compatable with the data
            if aDataObject is not None :
                assert(lonlatData["A"][0].shape == aDataObject.shape)
            if bDataObject is not None :
                assert(lonlatData["B"][0].shape == bDataObject.shape)
            # make composite valid mask
            allValidMask = ( lonlatData["A"][0].masks.valid_mask & lonlatData["A"][1].masks.valid_mask &
                             lonlatData["B"][0].masks.valid_mask & lonlatData["B"][1].masks.valid_mask )
            
            # build basemap, FUTURE, don't hard code so much of this stuff
            basemapObject = Basemap(llcrnrlon=lonRange[0], llcrnrlat=latRange[0], urcrnrlon=lonRange[1], urcrnrlat=latRange[1],
                                    resolution='i', area_thresh=10000, projection="merc")
            # TODO get all these variables outside the if statement
            
        
        LOG.info("Spawning plot window: " + imageType)
        
        plt.ion() # make sure interactive plotting is on
        
        # create whichever type of plot was asked for
        
        if (imageType == ORIGINAL_A) or (imageType == ORIGINAL_B) :
            
            # sort out some values based on which of the data sets we're showing
            data_object_to_use = aDataObject if (imageType == ORIGINAL_A) else bDataObject
            var_name_to_use    = aVarName    if (imageType == ORIGINAL_A) else bVarName
            file_char_to_use   = "A"         if (imageType == ORIGINAL_A) else "B"
            units_text_to_use  = aUnitsText  if (imageType == ORIGINAL_A) else bUnitsText
            oneD_color_to_use  = 'b'         if (imageType == ORIGINAL_A) else 'c'
            
            # if the data doesn't exist, we can't make this plot
            if data_object_to_use is None :
                
                raise ValueError(NO_DATA_MESSAGE)
            
            if dataForm == SIMPLE_2D :
                tempFigure = figures.create_simple_figure(data_object_to_use.data, var_name_to_use + "\nin File " + file_char_to_use,
                                                          invalidMask=~data_object_to_use.masks.valid_mask, colorMap=colorMapToUse,
                                                          colorbarLimits=rangeInfo, units=units_text_to_use)
            elif dataForm == MAPPED_2D :
                #_, tempLatObj, _ = self._getVariableInformation(file_char_to_use, variableName=self.dataModel.getLatitudeName (file_char_to_use))
                #_, tempLonObj, _ = self._getVariableInformation(file_char_to_use, variableName=self.dataModel.getLongitudeName(file_char_to_use))
                # TODO ***
                #tempFigure = figures.create_mapped_figure(data_object_to_use.data, tempLatObj.data, tempLonObj.data, baseMapInstance, boundingAxes, title,
                #          invalidMask=None, colorMap=None, tagData=None,
                #          dataRanges=None, dataRangeNames=None, dataRangeColors=None, units=None, **kwargs)
                pass
                
            elif dataForm == ONLY_1D :
                temp = [(data_object_to_use.data, ~data_object_to_use.masks.valid_mask, oneD_color_to_use, None, None, None)]
                tempFigure = figures.create_line_plot_figure(temp, var_name_to_use + "\n in File " + file_char_to_use)
            else :
                raise ValueError(UNKNOWN_DATA_FORM)
            
        elif (imageType == HISTOGRAM_A) or (imageType == HISTOGRAM_B) :
            
            # Note: histograms don't care about data format requested, they are histogram formatted
            
            # select the things that are file A or B specific
            file_desc_to_use   = "A"         if (imageType == HISTOGRAM_A) else "B"
            var_name_to_use    = aVarName    if (imageType == HISTOGRAM_A) else bVarName
            data_object_to_use = aDataObject if (imageType == HISTOGRAM_A) else bDataObject
            units_text_to_use  = aUnitsText  if (imageType == HISTOGRAM_A) else bUnitsText
            
            # if the data doesn't exist, we can't make this plot
            if data_object_to_use is None :
                
                raise ValueError(NO_DATA_MESSAGE)
            
            # build the histogram
            clean_data = data_object_to_use.data[data_object_to_use.masks.valid_mask]
            # TODO, should the range option be added here?
            tempFigure = figures.create_histogram(clean_data, DEFAULT_NUM_BINS, var_name_to_use + "\nin File " + file_desc_to_use,
                                                  "Value of data at a given point", "Number of points with a given value", units=units_text_to_use)
            
        elif imageType in COMPARISON_IMAGES :
            
            # if we're making the absolute or raw difference image, do that
            if (imageType == ABS_DIFF) or (imageType == RAW_DIFF) :
                
                # now choose between the raw and abs diff
                dataToUse   = diffData.diff_data_object.data
                titlePrefix = "Value of (Data File B - Data File A)\nfor "
                if imageType == ABS_DIFF :
                    dataToUse   = np.abs(dataToUse)
                    titlePrefix = "Absolute value of difference\nin "
                
                if dataForm == SIMPLE_2D :
                    tempFigure = figures.create_simple_figure(dataToUse, titlePrefix + aVarName,
                                                              invalidMask=~diffData.diff_data_object.masks.valid_mask,
                                                              colorMap=colorMapToUse, units=aUnitsText)
                elif dataForm == MAPPED_2D :
                    pass # TODO
                elif dataForm == ONLY_1D :
                    tempTitle = titlePrefix + aVarName
                    if aVarName != bVarName :
                        tempTitle = tempTitle + " / " + bVarName
                    temp = [(dataToUse, ~diffData.diff_data_object.masks.valid_mask, 'm', None, None, None)]
                    tempFigure = figures.create_line_plot_figure(temp, tempTitle)
                else :
                    raise ValueError(UNKNOWN_DATA_FORM)
                
            elif imageType == MISMATCH :
                
                mismatchMask = diffData.diff_data_object.masks.mismatch_mask
                if dataForm == SIMPLE_2D :
                    tempFigure = figures.create_simple_figure(aDataObject.data, "Areas of mismatch data\nin " + aVarName,
                                                              invalidMask=~aDataObject.masks.valid_mask, tagData=mismatchMask,
                                                              colorMap=figures.MEDIUM_GRAY_COLOR_MAP, units=aUnitsText)
                                                                # TODO, change colormap?
                elif dataForm == MAPPED_2D :
                    pass # TODO
                elif dataForm == ONLY_1D :
                    
                    temp = [(aDataObject.data, ~aDataObject.masks.valid_mask, 'k', None, mismatchMask, None)]
                    tempFigure = figures.create_line_plot_figure(temp, "Areas of mismatch data\nin " + aVarName)
                    
                else :
                    raise ValueError(UNKNOWN_DATA_FORM)
                
            elif imageType == HISTOGRAM :
                
                # Note: histograms don't care about data format requested, they are histogram formatted
                
                rawDiffDataClean = diffData.diff_data_object.data[diffData.diff_data_object.masks.valid_mask]
                titleText = ("Difference in\n" + aVarName) if (aVarName == bVarName) else str( "Value of\n" + bVarName + " - " + aVarName )
                # TODO, should the range option be added here?
                tempFigure = figures.create_histogram(rawDiffDataClean, DEFAULT_NUM_BINS, titleText,
                                                      "Value of (B - A) at each data point", "Number of points with a given difference", units=aUnitsText)
                
            elif (imageType == SCATTER) or (imageType == HEX_PLOT) :
                
                # Note: scatter and hex plots don't care about data format requested, they're scatter or hex plots
                
                tempCleanMask     = aDataObject.masks.valid_mask & bDataObject.masks.valid_mask
                aDataClean        = aDataObject.data[tempCleanMask]
                bDataClean        = bDataObject.data[tempCleanMask]
                
                if imageType == SCATTER :
                    
                    cleanMismatchMask = diffData.diff_data_object.masks.mismatch_mask[tempCleanMask]
                    figures.create_scatter_plot(aDataClean, bDataClean, "Value in File A vs Value in File B", 
                                        "File A Value in " + aVarName,
                                        "File B Value in " + bVarName,
                                        badMask=cleanMismatchMask, epsilon=self.dataModel.getEpsilon(),
                                        units_x=aUnitsText, units_y=bUnitsText)
                
                else:
                    
                    tempFigure = figures.create_hexbin_plot(aDataClean, bDataClean,
                                                    "Value in File A vs Value in File B",
                                                    "File A Value in " + aVarName,
                                                    "File B Value in " + bVarName,
                                                    epsilon=self.dataModel.getEpsilon(),
                                                    units_x=aUnitsText, units_y=bUnitsText)
            
        
        plt.draw()
