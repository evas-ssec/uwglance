#!/usr/bin/env python
# encoding: utf-8
"""
The model portion of the glance GUI code. 

Created by evas Oct 2011.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
"""

# these first two lines must stay before the pylab import
import matplotlib
matplotlib.use('Qt4Agg') # use the Qt Anti-Grain Geometry rendering engine

from pylab import *

import matplotlib.cm     as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import sys, os.path, logging
import numpy as np

import glance.data    as dataobjects
import glance.figures as figures
import glance.io      as io

LOG = logging.getLogger(__name__)

"""
The model handles the data in the Glance GUI. It not only stores data and handles
updates, it's also responsible for caching data and handling some logic related to
what data is used with overrides. 

The model allows outside objects to register to listen for data updates or for
errors. It's expected that error handlers will manage either logging or displaying
errors appropriately. 
"""

# constants for the possible image types
ORIGINAL_A = "Original A Data"
ORIGINAL_B = "Original B Data"
ABS_DIFF   = "Abs. Difference"
RAW_DIFF   = "Raw Difference"
HISTOGRAM  = "Histogram"
MISMATCH   = "Mismatch Areas"
SCATTER    = "Scatter Plot"
HEX_PLOT   = "Hex Plot"

# a list of all the image types, for convenience
IMAGE_TYPES = [ORIGINAL_A,
               ORIGINAL_B,
               ABS_DIFF,
               RAW_DIFF,
               HISTOGRAM,
               MISMATCH,
               SCATTER,
               HEX_PLOT
              ]

# the number of bins to use for histograms
DEFAULT_NUM_BINS = 50

class _FileModelData (object) :
    """
    This object is meant to be used internally by the GUI model. The model is going to mess with the
    data directly. This is not the best practice ever, but I like it more than using a dictionary straight.
    This data object includes the following:
    
    self.file             - the FileInfo object representing this file, can be used to load more information later
    self.variable         - the name of the selected variable
    self.doOverride       - false if the default fill value should be used, true otherwise
    self.fillValue        - the fill value (may be different than the default), should be used when override is true
    self.defaultFillValue - the default fill value, should be used when override is false
    self.ALL_VARIABLES    - a list of all the variable names in the file
    self.var_data         - the data contained in the currently selected variable (or None if no variable is selected)
    self.var_attrs        - a dictionary of the variable attributes, keyed with the attribute names and containing their values
    
    TODO, eventually replace this with a better data object and/or some combination of data and settings objects?
    """
    
    def __init__(self, file_object=None, variable_selection=None, do_override=False, fill_value=None, default_fill_value=None,
                 variables_list=None, variable_data=None, variable_attributes={ }) :
        """
        create a set of model data, using the data passed in
        
        Note: the file_object is intended to be a FileInfo object from glance.data
        """
        
        self.file             = file_object
        self.variable         = variable_selection
        self.doOverride       = do_override
        self.fillValue        = fill_value
        self.defaultFillValue = default_fill_value
        self.ALL_VARIABLES    = variables_list
        
        self.var_data         = variable_data
        self.var_attrs        = variable_attributes

class GlanceGUIModel (object) :
    """
    This is the main model that handles the information behind the glance GUI.
    It includes:
    
    self.fileData      - a dictionary with _FileModelData objects indexed by the two file prefixes
    self.epsilon       - the epsilon value for comparison
    self.imageType     - the image type that should be created when files are compared
    self.dataListeners - objects that want to be notified when data changes
    self.errorHandlers - objects that want to be notified when there's a serious error
    """
    
    def __init__ (self) :
        """
        set up the basic model with our initial default values and empty listener lists
        """
        
        # set up the file related data structures
        self.fileData      = { }
        self.fileData["A"] = _FileModelData( )
        self.fileData["B"] = _FileModelData( )
        
        # general settings
        self.epsilon       = 0.0
        self.imageType     = None
        
        # select the first image type as a default
        self.imageType     = IMAGE_TYPES[0]
        
        # this represents all the people who want to hear about data updates
        # these people can register and will get data related messages
        self.dataListeners = [ ]
        # this represents all the people who want to hear when we have errors
        # that the user should know about
        # these people can register and will get error related messages
        self.errorHandlers = [ ]
    
    def loadNewFile (self, filePrefix, newFilePath) :
        """
        load up a new file based on the prefix and path given
        """
        
        # check to see if we were given a valid file path
        if (newFilePath is None) or (len(newFilePath) is 0) :
            LOG.debug("No file selected. Aborting load.")
            return # if there's no path, we can't load anything
        
        # attempt to open the file
        try :
            newFile      = dataobjects.FileInfo(str(newFilePath))
        except KeyError :
            newFile = None
            messageTemp = "Unable to open file. Glance was not able to determine the file type for\n " + newFilePath + "\nor cannot process that type of file."
            for errorHandler in self.errorHandlers :
                errorHandler.handleWarning(messageTemp)
        
        # if we couldn't get a valid file, stop now
        if newFile is None :
            return
        
        # get the list of variables, and pick one
        variableList = sorted(newFile.file_object()) # gets a list of all the variables in the file
        tempVariable = variableList[0]
        fillValue    = newFile.file_object.missing_value(tempVariable)
        
        # save all of the data related to this file for later use
        self.fileData[filePrefix].file               = newFile
        self.fileData[filePrefix].variable           = tempVariable
        self.fileData[filePrefix].doOverride         = False
        self.fileData[filePrefix].defaultFillValue   = fillValue
        self.fileData[filePrefix].fillValue          = fillValue
        self.fileData[filePrefix].ALL_VARIABLES      = variableList
        
        # get the size of the currently selected variable TODO, is it possible to do this without loading the variable?
        tempShape = self._load_variable_data(filePrefix, str(self.fileData[filePrefix].variable))
        
        # get the variable's attributes TODO, does this work on all types of files? (FIXME no, make a general method in io!)
        self.fileData[filePrefix].var_attrs          = newFile.file_object.get_variable_object(str(tempVariable)).attributes()
        
        # Now tell our data listeners that the file data changed
        for dataListener in self.dataListeners :
            LOG.debug("Sending update for file " + filePrefix + " with loaded data.")
            dataListener.fileDataUpdate(filePrefix, newFile.path, tempVariable, False, fillValue, str(tempShape),
                                        variable_list=variableList, attribute_list=self.fileData[filePrefix].var_attrs)
    
    def _load_variable_data (self, file_prefix, variable_name) :
        """
        Load up new variable data, saving it to our fileData structure
        return the shape of the data for convenience
        
        TODO, can this be handled as a background task in the future?
        """
        
        self.fileData[file_prefix].var_data = self.fileData[file_prefix].file.file_object[variable_name]
        return self.fileData[file_prefix].var_data.shape
    
    def sendGeneralSettingsData (self) :
        """
        send off the general settings data that's not related to the individual files
        """
        
        # let each of our listeners know about the general data
        for dataListener in self.dataListeners :
            dataListener.updateEpsilon(self.epsilon)
            dataListener.updateImageTypes(self.imageType, list=IMAGE_TYPES)
    
    def updateFileDataSelection (self, file_prefix, newVariableText=None, newOverrideValue=None, newFillValue=np.nan) :
        """
        someone has updated one or more of the file related data selections
        
        Note: if an input value is left at it's default (None or nan) then it's assumed that it was not externally changed
        """
        
        didUpdate = False
        
        # update the variable selection if needed
        if (newVariableText is not None) and (newVariableText != self.fileData[file_prefix].variable) :
            if newVariableText in self.fileData[file_prefix].ALL_VARIABLES :
                LOG.debug("Setting file " + file_prefix + " variable selection to: " + newVariableText)
                self.fileData[file_prefix].variable = newVariableText
                didUpdate = True
                
                # load the data for this new variable
                self._load_variable_data(file_prefix, str(newVariableText))
                
                # get the variable's attributes TODO, does this work on all types of files? (FIXME nope, make a general method in io!)
                self.fileData[file_prefix].var_attrs = self.fileData[file_prefix].file.file_object.get_variable_object(str(newVariableText)).attributes()
                
                # the new fill value should be loaded and the override should be cleared
                self.fileData[file_prefix].doOverride        = False
                self.fileData[file_prefix].defaultFillValue  = self.fileData[file_prefix].file.file_object.missing_value(str(newVariableText))
                self.fileData[file_prefix].fillValue         = self.fileData[file_prefix].defaultFillValue
        
        # update the override selection if needed
        if newOverrideValue is not None :
            LOG.debug("Setting file " + file_prefix + " override selection to: " + str(newOverrideValue))
            self.fileData[file_prefix].doOverride = newOverrideValue
            didUpdate = True
        
        # update the fill value if needed
        if newFillValue is not np.nan :
            LOG.debug("Setting file " + file_prefix + " fill value to: " + str(newFillValue))
            self.fileData[file_prefix].fillValue = newFillValue
            didUpdate = True
        
        # let our data listeners know about any changes
        if didUpdate :
            for listener in self.dataListeners :
                listener.fileDataUpdate(file_prefix, self.fileData[file_prefix].file.path,  self.fileData[file_prefix].variable,
                                                     self.fileData[file_prefix].doOverride, self._select_fill_value(file_prefix),
                                                     str(self.fileData[file_prefix].var_data.shape), attribute_list=self.fileData[file_prefix].var_attrs)
    
    def _select_fill_value (self, file_prefix) :
        """
        which fill value should currently be used?
        """
        return self.fileData[file_prefix].fillValue if self.fileData[file_prefix].doOverride else self.fileData[file_prefix].defaultFillValue
    
    def updateSettingsDataSelection (self, newEpsilonValue=np.nan, newImageType=None) :
        """
        someone has changed one or more of the general settings related data values
        
        Note: if an input value is left at it's default (None or nan) then it's assumed that it was not externally changed
        """
        
        didUpdate = False
        
        # update the epsilon if needed
        if newEpsilonValue is not np.nan :
            LOG.debug("Setting epsilon to: " + str(newEpsilonValue))
            self.epsilon = newEpsilonValue
            didUpdate = True
        
        # update the image type if needed
        if (newImageType is not None) and (newImageType != self.imageType) :
            if newImageType in IMAGE_TYPES :
                LOG.debug("Setting image type to: " + newImageType)
                self.imageType = newImageType
                didUpdate = True
        
        # let our data listeners know about any changes
        if didUpdate :
            for listener in self.dataListeners :
                listener.updateEpsilon(self.epsilon)
                listener.updateImageTypes(self.imageType)
    
    def spawnPlotWithCurrentInfo (self) :
        """
        create a matplotlib plot using the current model information
        
        TODO, move this into some sort of figure manager model object/module?
        """
        
        LOG.info ("Preparing variable data for plotting...")
        
        aData = self.fileData["A"].var_data
        bData = self.fileData["B"].var_data
        
        message = None
        
        # TODO, right now this does not take into account the fact that the two "original" plots only need one of the two variables
        # minimally validate the data
        if (aData is None) or (bData is None) :
            message = ("Data for requested files was not available. " +
                       "Please load or reload files and try again.")
        # check to see if the two variables have the same shape of data
        elif aData.shape != bData.shape :
            message = (self.fileData["A"].variable + ' / ' + self.fileData["B"].variable + ' ' + 
                       'could not be compared because the data for these variables does not match in shape ' +
                       'between the two files (file A data shape: ' + str(aData.shape) + '; file B data shape: '
                       + str(bData.shape) + ').')
        # if the data isn't valid, stop now
        if message is not None :
            for errorHandler in self.errorHandlers :
                errorHandler.handleWarning(message)
            # we can't make any images from this data, so just return
            return
        
        tempAFillMask = np.zeros(aData.shape, dtype=np.bool)
        tempAFillMask[aData == float(self._select_fill_value("A"))] = True
        tempBFillMask = np.zeros(bData.shape, dtype=np.bool)
        tempBFillMask[bData == float(self._select_fill_value("B"))] = True
        tempCleanMask = tempAFillMask | tempBFillMask
        
        aDataClean = aData[~tempCleanMask]
        bDataClean = bData[~tempCleanMask]
        cleanMismatchMask = np.zeros(aDataClean.shape, dtype=np.bool)
        cleanMismatchMask[abs(aDataClean - bDataClean) > float(self.epsilon)] = True
        
        rawDiffDataClean = bDataClean - aDataClean
        
        # pull the units information
        aUnits = self.fileData["A"].file.file_object.get_attribute(str(self.fileData["A"].variable), io.UNITS_CONSTANT)
        bUnits = self.fileData["B"].file.file_object.get_attribute(str(self.fileData["B"].variable), io.UNITS_CONSTANT)
        
        LOG.info("Spawning plot window: " + self.imageType)
        
        plt.ion() # make sure interactive plotting is on
        
        # create the plot
        
        if   self.imageType == ORIGINAL_A :
            
            tempFigure = figures.create_simple_figure(aData, self.fileData["A"].variable + "\nin File A",
                                                      invalidMask=tempAFillMask, colorMap=cm.jet, units=aUnits)
            
        elif self.imageType == ORIGINAL_B :
            
            tempFigure = figures.create_simple_figure(bData, self.fileData["B"].variable + "\nin File B",
                                                      invalidMask=tempBFillMask, colorMap=cm.jet, units=bUnits)
            
        elif self.imageType == ABS_DIFF :
            
            tempFigure = figures.create_simple_figure(np.abs(bData - aData), "Absolute value of difference\nin " + self.fileData["A"].variable,
                                                      invalidMask=tempCleanMask, colorMap=cm.jet, units=aUnits)
            
        elif self.imageType == RAW_DIFF :
            
            tempFigure = figures.create_simple_figure(bData - aData, "Value of (Data File B - Data File A)\nfor " + self.fileData["A"].variable,
                                                      invalidMask=tempCleanMask, colorMap=cm.jet, units=aUnits)
            
        elif self.imageType == HISTOGRAM :
            
            tempFigure = figures.create_histogram(rawDiffDataClean, DEFAULT_NUM_BINS, "Difference in\n" + self.fileData["A"].variable,
                                                  "Value of (B - A) at each data point", "Number of points with a given difference", units=aUnits)
            
        elif self.imageType == MISMATCH :
            
            mismatchMask = np.zeros(aData.shape, dtype=bool)
            mismatchMask[abs(aData - bData) > float(self.epsilon)] = True
            mismatchMask[tempCleanMask] = False
            mismatchMask[tempBFillMask ^ tempAFillMask] = True
            
            tempFigure = figures.create_simple_figure(aData, "Areas of mismatch data\nin " + self.fileData["A"].variable,
                                                      invalidMask=tempCleanMask, tagData=mismatchMask, colorMap=figures.MEDIUM_GRAY_COLOR_MAP, units=aUnits)
            
        elif self.imageType == SCATTER :
            
            figures.create_scatter_plot(aDataClean, bDataClean, "Value in File A vs Value in File B", 
                                        "File A Value in " + self.fileData["A"].variable,
                                        "File B Value in " + self.fileData["B"].variable,
                                        badMask=cleanMismatchMask, epsilon=float(self.epsilon),
                                        units_x=aUnits, units_y=bUnits)
            
        elif self.imageType == HEX_PLOT :
            
            tempFigure = figures.create_hexbin_plot(aDataClean, bDataClean,
                                                    "Value in File A vs Value in File B",
                                                    "File A Value in " + self.fileData["A"].variable,
                                                    "File B Value in " + self.fileData["B"].variable,
                                                    epsilon=float(self.epsilon),
                                                    units_x=aUnits, units_y=bUnits)
        
        plt.draw()
    
    def registerDataListener (self, objectToRegister) :
        """
        add the given object to our list of data listeners
        """
        
        if objectToRegister not in self.dataListeners :
            self.dataListeners.append(objectToRegister)
    
    def registerErrorHandler (self, objectToRegister) :
        """
        add the given object to our list of error handlers
        """
        
        if objectToRegister not in self.errorHandlers :
            self.errorHandlers.append(objectToRegister)

if __name__=='__main__':
    import doctest
    doctest.testmod()
