#!/usr/bin/env python
# encoding: utf-8
"""
This module manages creating figures for the Glance GUI.

Created by evas Nov 2011.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
"""

# these first two lines must stay before the pylab import
import matplotlib
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
    
    def spawnPlot (self) :
        """
        create a matplotlib plot using the current model information
        """
        
        imageType = self.dataModel.getImageType()
        
        LOG.info ("Preparing variable data for plotting...")
        
        # get Variable names
        aVarName    = self.dataModel.getVariableName("A")
        bVarName    = self.dataModel.getVariableName("B")
        
        # get Data objects
        aDataObject = self.dataModel.getVariableData("A", aVarName)
        bDataObject = self.dataModel.getVariableData("B", bVarName)
        
        # TODO, this ignores the fact that the "original" plots don't need two sets of data
        message = dataobjects.DiffInfoObject.verifyDataCompatability (aDataObject, bDataObject, aVarName, bVarName)
        
        # if the data isn't valid, stop now
        if message is not None :
            for errorHandler in self.errorHandlers :
                errorHandler.handleWarning(message)
            # we can't make any images from this data, so just return
            return
        
        # compare our data
        diffData = dataobjects.DiffInfoObject(aDataObject, bDataObject, epsilonValue=self.dataModel.getEpsilon())
        
        # get units text for display
        aUnitsText  = self.dataModel.getUnitsText("A", aVarName)
        bUnitsText  = self.dataModel.getUnitsText("B", bVarName)
        
        LOG.info("Spawning plot window: " + imageType)
        
        plt.ion() # make sure interactive plotting is on
        
        # create the plot
        
        if   imageType == model.ORIGINAL_A :
            
            tempFigure = figures.create_simple_figure(aDataObject.data, aVarName + "\nin File A",
                                                      invalidMask=aDataObject.masks.missing_mask, colorMap=cm.jet, units=aUnitsText)
            
        elif imageType == model.ORIGINAL_B :
            
            tempFigure = figures.create_simple_figure(bDataObject.data, bVarName + "\nin File B",
                                                      invalidMask=bDataObject.masks.missing_mask, colorMap=cm.jet, units=bUnitsText)
            
        elif imageType == model.ABS_DIFF :
            
            tempFigure = figures.create_simple_figure(np.abs(diffData.diff_data_object.data), "Absolute value of difference\nin " + aVarName,
                                                      invalidMask=~diffData.diff_data_object.masks.valid_mask, colorMap=cm.jet, units=aUnitsText)
            
        elif imageType == model.RAW_DIFF :
            
            tempFigure = figures.create_simple_figure(diffData.diff_data_object.data, "Value of (Data File B - Data File A)\nfor " + aVarName,
                                                      invalidMask=~diffData.diff_data_object.masks.valid_mask, colorMap=cm.jet, units=aUnitsText)
            
        elif imageType == model.HISTOGRAM :
            
            rawDiffDataClean = diffData.diff_data_object.data[diffData.diff_data_object.masks.valid_mask]
            tempFigure = figures.create_histogram(rawDiffDataClean, DEFAULT_NUM_BINS, "Difference in\n" + aVarName,
                                                  "Value of (B - A) at each data point", "Number of points with a given difference", units=aUnitsText)
            
        elif imageType == model.MISMATCH :
            
            mismatchMask = diffData.diff_data_object.masks.mismatch_mask
            tempFigure = figures.create_simple_figure(aDataObject.data, "Areas of mismatch data\nin " + aVarName,
                                                      invalidMask=aDataObject.masks.missing_mask, tagData=mismatchMask,
                                                      colorMap=figures.MEDIUM_GRAY_COLOR_MAP, units=aUnitsText)
            
        elif imageType == model.SCATTER :
            
            tempCleanMask     = aDataObject.masks.missing_mask | bDataObject.masks.missing_mask
            aDataClean        = aDataObject.data[~tempCleanMask]
            bDataClean        = bDataObject.data[~tempCleanMask]
            cleanMismatchMask = diffData.diff_data_object.masks.mismatch_mask[~tempCleanMask]
            figures.create_scatter_plot(aDataClean, bDataClean, "Value in File A vs Value in File B", 
                                        "File A Value in " + aVarName,
                                        "File B Value in " + bVarName,
                                        badMask=cleanMismatchMask, epsilon=self.dataModel.getEpsilon(),
                                        units_x=aUnitsText, units_y=bUnitsText)
            
        elif imageType == model.HEX_PLOT :
            
            tempCleanMask     = aDataObject.masks.missing_mask | bDataObject.masks.missing_mask
            aDataClean        = aDataObject.data[~tempCleanMask]
            bDataClean        = bDataObject.data[~tempCleanMask]
            tempFigure = figures.create_hexbin_plot(aDataClean, bDataClean,
                                                    "Value in File A vs Value in File B",
                                                    "File A Value in " + aVarName,
                                                    "File B Value in " + bVarName,
                                                    epsilon=self.dataModel.getEpsilon(),
                                                    units_x=aUnitsText, units_y=bUnitsText)
        
        plt.draw()
