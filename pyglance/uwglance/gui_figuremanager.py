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

LOG = logging.getLogger(__name__)

# the number of bins to use for histograms
DEFAULT_NUM_BINS = 50

NO_DATA_MESSAGE = "Requested data was not available or did not exist."

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
    
    def _getVariableInformation (self, filePrefix) :
        """
        Pull the name, data, and units for the variable currently selected in the given file prefix
        """
        
        selectedVariableName = self.dataModel.getVariableName(filePrefix)
        dataObject           = self.dataModel.getVariableData(filePrefix, selectedVariableName)
        unitsText            = self.dataModel.getUnitsText(filePrefix, selectedVariableName)
        
        if dataObject is not None :
            dataObject.self_analysis()
        
        return selectedVariableName, dataObject, unitsText
    
    def _getVariableInfoSmart (self, filePrefix, imageType) :
        """
        if appropriate for the image type, get information on the variable, otherwise return None's
        """
        
        varName, dataObject, unitsText = None, None, None
        
        # only load the data if it will be needed for the plot
        if ((imageType == model.ORIGINAL_A) and (filePrefix == "A") or
            (imageType == model.ORIGINAL_B) and (filePrefix == "B") or
            (imageType in model.COMPARISON_IMAGES)) :
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
        if imageType in model.COMPARISON_IMAGES :
            
            # check to see if our data is minimally compatable; this call may raise an IncompatableDataObjects exception
            dataobjects.DiffInfoObject.verifyDataCompatability (dataObjectA, dataObjectB, varNameA, varNameB)
            
            # compare our data
            diffObject = dataobjects.DiffInfoObject(dataObjectA, dataObjectB,
                                                    epsilonValue=epsilon_value, epsilonPercent=epsilon_percent)
        
        return diffObject
    
    def spawnPlot (self) :
        """
        create a matplotlib plot using the current model information
        
        this method may raise an IncompatableDataObjects exception if the a and b data
        are completely incomparable
        this method may also raise a ValueError if the data could not be plotted
        for reasons not encompassed by an IncompatableDataObjects exception
        """
        
        imageType = self.dataModel.getImageType()
        dataForm  = self.dataModel.getDataForm()
        
        LOG.info ("Preparing variable data for plotting...")
        
        aVarName, aDataObject, aUnitsText = self._getVariableInfoSmart("A", imageType)
        bVarName, bDataObject, bUnitsText = self._getVariableInfoSmart("B", imageType)
        
        diffData = self._buildDiffInfoObjectSmart(imageType,
                                                  aDataObject, bDataObject,
                                                  aVarName,    bVarName,
                                                  epsilon_value=self.dataModel.getEpsilon(),
                                                  epsilon_percent=self.dataModel.getEpsilonPercent())
        
        LOG.info("Spawning plot window: " + imageType)
        
        plt.ion() # make sure interactive plotting is on
        
        # create the plot
        
        if   imageType == model.ORIGINAL_A :
            
            # if the data doesn't exist, we can't make this plot
            if aDataObject is None :
                
                raise ValueError(NO_DATA_MESSAGE)
            
            tempFigure = figures.create_simple_figure(aDataObject.data, aVarName + "\nin File A",
                                                      invalidMask=~aDataObject.masks.valid_mask, colorMap=cm.jet, units=aUnitsText)
            # TODO, this is a hack to show AWIPS data, make an option for this at some point on the second tab
            #tempFigure = figures.create_simple_figure(aDataObject.data.astype(np.uint8), aVarName + "\nin File A",
            #                                          invalidMask=~aDataObject.masks.valid_mask, colorMap=cm.bone, units=aUnitsText)
            
            
        elif imageType == model.ORIGINAL_B :
            
            # if the data doesn't exist, we can't make this plot
            if bDataObject is None :
                
                raise ValueError(NO_DATA_MESSAGE)
            
            tempFigure = figures.create_simple_figure(bDataObject.data, bVarName + "\nin File B",
                                                      invalidMask=~bDataObject.masks.valid_mask, colorMap=cm.jet, units=bUnitsText)
            
        elif imageType in model.COMPARISON_IMAGES :
            
            # if we're making the absolute or raw difference image, do that
            if (imageType == model.ABS_DIFF) or (imageType == model.RAW_DIFF) :
                
                # now choose between the raw and abs diff
                dataToUse   = diffData.diff_data_object.data
                titlePrefix = "Value of (Data File B - Data File A)\nfor "
                if imageType == model.ABS_DIFF :
                    dataToUse   = np.abs(dataToUse)
                    titlePrefix = "Absolute value of difference\nin "
                
                tempFigure = figures.create_simple_figure(dataToUse, titlePrefix + aVarName,
                                                          invalidMask=~diffData.diff_data_object.masks.valid_mask,
                                                          colorMap=cm.jet, units=aUnitsText)
            elif imageType == model.MISMATCH :
                
                mismatchMask = diffData.diff_data_object.masks.mismatch_mask
                tempFigure = figures.create_simple_figure(aDataObject.data, "Areas of mismatch data\nin " + aVarName,
                                                          invalidMask=~aDataObject.masks.valid_mask, tagData=mismatchMask,
                                                          colorMap=figures.MEDIUM_GRAY_COLOR_MAP, units=aUnitsText)
                
            elif imageType == model.HISTOGRAM :
                
                rawDiffDataClean = diffData.diff_data_object.data[diffData.diff_data_object.masks.valid_mask]
                tempFigure = figures.create_histogram(rawDiffDataClean, DEFAULT_NUM_BINS, "Difference in\n" + aVarName,
                                                      "Value of (B - A) at each data point", "Number of points with a given difference", units=aUnitsText)
                
            elif (imageType == model.SCATTER) or (imageType == model.HEX_PLOT) :
                
                tempCleanMask     = aDataObject.masks.valid_mask & bDataObject.masks.valid_mask
                aDataClean        = aDataObject.data[tempCleanMask]
                bDataClean        = bDataObject.data[tempCleanMask]
                
                if imageType == model.SCATTER :
                    
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
