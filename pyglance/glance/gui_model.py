#!/usr/bin/env python
# encoding: utf-8
"""
The model portion of the glance GUI code. 

Created by evas Oct 2011.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
"""

import logging
import numpy as np
from   os import path

import glance.data as dataobjects
import glance.io   as io
from   glance.gui_constants import *

LOG = logging.getLogger(__name__)

"""
The model handles the data in the Glance GUI. It not only stores data and handles
updates, it's also responsible for caching data and handling some logic related to
what data is used with overrides. 

The model allows outside objects to register to listen for data updates or for
errors. It's expected that error handlers will manage either logging or displaying
errors appropriately. 
"""

class UnableToReadFile(Exception):
    """
    An exception to be used when glance could not read a file
    """
    
    def __init__(self, filePath):
        """
        create this exception, giving a message based on the file path
        """
        
        self.message = str("Unable to open file. Glance was not able to determine the file type for "
                           + filePath + " or cannot process that type of file.")
    
    def __str__(self):
        return self.message

class _FileModelData (object) :
    """
    This object is meant to be used internally by the GUI model. The model is going to mess with the
    data directly. This is not the best practice ever, but I like it more than using a dictionary straight.
    This data object includes the following:
    
    self.file             - the FileInfo object representing this file, can be used to load more information later
    self.variable         - the name of the selected variable
    self.var_data_cache   - a cache of all the variable data that has been loaded for this file
                            (stored in dataobjects.DataObject objects), keyed by variable name
    self.var_attrs_cache  - a cache of variable attributes (keyed by variable name), each set of attributes is a dictionary,
                            keyed with the attribute names and containing their values
    self.ALL_VARIABLES    - a list of all the variable names in the file
    """
    
    def __init__(self, file_object=None, variable_selection=None, do_override=False, fill_value=None, default_fill_value=None,
                 latitude_name=None, longitude_name=None,
                 variables_list=None, variable_data=None, variable_attributes=None) :
        """
        create a set of model data, using the data passed in
        
        Note: the file_object is intended to be a FileInfo object from glance.data
        """
        
        self.file             = file_object
        self.variable         = str(variable_selection)
        self.latitude         = latitude_name
        self.longitude        = longitude_name
        self.ALL_VARIABLES    = variables_list
        
        self.var_data_cache   = { }
        self.var_attrs_cache  = { }
        if variable_selection is not None :
            self.var_data_cache[variable_selection]  = dataobjects.DataObject(variable_data, fillValue= fill_value,
                                                                              overrideFillValue=do_override,
                                                                              defaultFillValue=default_fill_value)
            self.var_attrs_cache[variable_selection] = variable_attributes

class GlanceGUIModel (object) :
    """
    This is the main model that handles the information behind the glance GUI.
    It includes:
    
    self.fileData       - a dictionary with _FileModelData objects indexed by the two file prefixes
    self.epsilon        - the epsilon value for comparison
    self.epsilonPercent - the epsilon percent value for comparison
    self.llEpsilon      - the epsilon used for judging the longitude and latitude data
    self.imageType      - the image type that should be created when files are compared
    self.colormap       - the colormap to use for plotting (if one is needed)
    self.dataForm       - the form the data should be considered to be
    self.useSharedRange - True if images for the original data should be displayed in a range
                          that includeds the data from both files, False if not
    self.plotGeoTiffAsRGB - True if multi-channel geotiff images should be treated as RGB or RGBA images
    self.hideMismatchNav - True if data corresponding to navigation that's more different
                           than self.llEpsilon should be hidden when plotting, False if not
    
    self.fileSettings - a dictionary of settings specific to a file, in the form:
                                {
                                    "doRange":  True or False for whether or not the range of the data should be restricted,
                                    "minRange": The minimum acceptable value if the range is restricted,
                                    "maxRange": The maximum acceptable value if the range is restricted,
                                    "isAWIPS":  True or False for whether or not the data is in AWIPS format
                                }
    
    self.dataListeners  - objects that want to be notified when data changes
    self.errorHandlers  - objects that want to be notified when there's a serious error
    """
    
    DO_RANGE  = "doRange"
    MIN_RANGE = "minRange"
    MAX_RANGE = "maxRange"
    IS_AWIPS  = "isAWIPS"
    
    def __init__ (self) :
        """
        set up the basic model with our initial default values and empty listener lists
        """
        
        # set up the file related data structures
        self.fileData      = { }
        self.fileData[A_CONST] = _FileModelData( )
        self.fileData[B_CONST] = _FileModelData( )
        
        # general settings
        self.epsilon          = 0.0
        self.epsilonPercent   = None
        self.llEpsilon        = None
        self.imageType        = IMAGE_TYPES[0]
        self.colormap         = COLORMAP_NAMES[0]
        self.dataForm         = SIMPLE_2D
        self.useSharedRange   = False
        self.plotGeoTiffAsRGB = False
        self.hideMismatchNav  = False
        
        # This is obviously only going to work for these two prefixes, would need
        # to add a fully formed sub-class to make this more general
        self.fileSettings = { }
        self.fileSettings[A_CONST] = \
                               {
                                GlanceGUIModel.DO_RANGE:  False,
                                GlanceGUIModel.MIN_RANGE: None,
                                GlanceGUIModel.MAX_RANGE: None,
                                GlanceGUIModel.IS_AWIPS:  False
                               }
        self.fileSettings[B_CONST] = \
                               {
                                GlanceGUIModel.DO_RANGE:  False,
                                GlanceGUIModel.MIN_RANGE: None,
                                GlanceGUIModel.MAX_RANGE: None,
                                GlanceGUIModel.IS_AWIPS:  False
                               }
        
        # this represents all the people who want to hear about data updates
        # these people can register and will get data related messages
        self.dataListeners  = [ ]
        # this represents all the people who want to hear when we have errors
        # that the user should know about
        # these people can register and will get error related messages
        self.errorHandlers  = [ ]
    
    def loadNewFile (self, filePrefix, newFilePath) :
        """
        load up a new file based on the prefix and path given
        
        may raise an UnableToReadFile exception if the given file path is non-null
        but the resulting file cannot be parsed by glance
        """
        
        # check to see if we were given a valid file path
        if (newFilePath is None) or (len(newFilePath) is 0) :
            LOG.debug("No file selected. Aborting load.")
            return # if there's no path, we can't load anything
        
        # attempt to open the file
        try :
            newFile = dataobjects.FileInfo(str(newFilePath))
        except KeyError :
            raise UnableToReadFile(newFilePath)
        
        # reset our caches
        self._resetCaches(filePrefix)
        
        # get the list of variables, and pick one
        variableList = sorted(newFile.file_object()) # gets a list of all the variables in the file
        tempVariable = str(variableList[0])
        tempIndex = 0
        # don't automatically select the nwp variables if possible
        while (tempVariable.find("nwp_") >= 0) :
            tempIndex    = tempIndex + 1
            tempVariable = variableList[tempIndex]
        
        
        LOG.debug ("selected variable: " + str(tempVariable))
        
        # save all of the data related to this file for later use
        self.fileData[filePrefix].file          = newFile
        self.fileData[filePrefix].variable      = tempVariable
        self.fileData[filePrefix].ALL_VARIABLES = variableList
        
        # set the longitude and latitude names, using the defaults if they exist
        self.fileData[filePrefix].latitude      = tempVariable
        self.fileData[filePrefix].longitude     = tempVariable
        if DEFAULT_LATITUDE  in variableList :
            self.fileData[filePrefix].latitude  = DEFAULT_LATITUDE
        if DEFAULT_LONGITUDE in variableList :
            self.fileData[filePrefix].longitude = DEFAULT_LONGITUDE
        # make sure the longitude and latitude are loaded into our local cache
        # TODO, it would be good to background this task at some point
        _ = self._load_variable_data(filePrefix, self.fileData[filePrefix].latitude)
        _ = self._load_variable_data(filePrefix, self.fileData[filePrefix].longitude) 
        
        # load info on the current variable
        tempDataObj = self._load_variable_data(filePrefix, tempVariable)
        
        # get the variable's attributes
        tempAttrs = self._load_variable_attributes (filePrefix, tempVariable)
        
        # Now tell our data listeners that the file data changed
        for dataListener in self.dataListeners :
            LOG.debug("Sending update for file " + filePrefix + " with loaded data.")
            dataListener.fileDataUpdate(filePrefix, newFile.path, tempVariable, tempDataObj.override_fill_value,
                                        self._select_fill_value(filePrefix), str(tempDataObj.data.shape),
                                        variable_list=variableList, attribute_list=tempAttrs)
            dataListener.updateSelectedLatLon(filePrefix,
                                              self.fileData[filePrefix].latitude,
                                              self.fileData[filePrefix].longitude,
                                              lonlatList=variableList)
    
    def _load_variable_attributes (self, file_prefix, variable_name) :
        """
        Load the attributes for for a given file name, saving them to the
        file data structure as appropriate
        
        return the loaded attributes for convenience
        """
        
        variable_name = str(variable_name)
        
        tempAttrs = None
        # if we can load the attributes from the cache, do that otherwise get them from the file
        if variable_name in  self.fileData[file_prefix].var_attrs_cache.keys() :
            tempAttrs = self.fileData[file_prefix].var_attrs_cache[variable_name]
        else :
            tempAttrs = self.fileData[file_prefix].file.file_object.get_variable_attributes(variable_name)
            # cache these for later use
            self.fileData[file_prefix].var_attrs_cache[variable_name] = tempAttrs
        
        return tempAttrs
    
    def _load_variable_data (self, file_prefix, variable_name) :
        """
        Load up new variable data, saving it to our fileData structure
        return the shape of the data for convenience
        
        TODO, can this be handled as a background task in the future?
        """
        
        variable_name = str(variable_name)
        
        tempData = None
        # if we have a cached version of this variable, use that, otherwise, load it from the file
        if variable_name in self.fileData[file_prefix].var_data_cache.keys() :
            LOG.debug ("Loading " + str(file_prefix) + " file cached variable: " + str(variable_name))
            tempData = self.fileData[file_prefix].var_data_cache[variable_name]
        else :
            LOG.debug ("Loading " + str(file_prefix) + " file variable from file: " + str(variable_name))
            tempRawData  = self.fileData[file_prefix].file.file_object[variable_name]
            tempFillVal  = self.fileData[file_prefix].file.file_object.missing_value(variable_name)
            tempOverride = False
            # also save this new data in our cache TODO, this won't save the other data we need?
            tempData = dataobjects.DataObject(tempRawData, fillValue=tempFillVal,
                                              overrideFillValue=tempOverride,
                                              defaultFillValue=tempFillVal)
            self.fileData[file_prefix].var_data_cache[variable_name] = tempData
        
        return tempData
    
    def _resetCaches (self, file_prefix) :
        """
        Clear the two internal caches
        """
        self.fileData[file_prefix].var_data_cache  = { }
        self.fileData[file_prefix].var_attrs_cache = { }
    
    def sendGeneralSettingsData (self) :
        """
        send off the general settings data that's not related to the individual files
        """
        
        # let each of our listeners know about the general data
        for dataListener in self.dataListeners :
            dataListener.updateEpsilon(self.epsilon)
            dataListener.updateEpsilonPercent(self.epsilonPercent)
            dataListener.updateLLEpsilon(self.llEpsilon)
            dataListener.updateImageTypes(self.imageType, list=IMAGE_TYPES)
            dataListener.updateColormaps(self.colormap, list=COLORMAP_NAMES)
            dataListener.updateDataForms(self.dataForm, list=DATA_FORMS)
            dataListener.updateUseSharedRange(self.useSharedRange)
            dataListener.updatePlotGeoTiffAsRGB(self.plotGeoTiffAsRGB)
            dataListener.updateHideMismatchNav(self.hideMismatchNav)
        
        self.sendFileSettings(A_CONST)
        self.sendFileSettings(B_CONST)
    
    def sendFileSettings (self, file_prefix) :
        """
        send out settings data that's related to the individual files but not data selections
        """
        
        # let our data listeners know about these values
        for listener in self.dataListeners :
            listener.updateDoRestrictRange (file_prefix, self.fileSettings[file_prefix][GlanceGUIModel.DO_RANGE ])
            listener.updateRestrictRangeMin(file_prefix, self.fileSettings[file_prefix][GlanceGUIModel.MIN_RANGE])
            listener.updateRestrictRangeMax(file_prefix, self.fileSettings[file_prefix][GlanceGUIModel.MAX_RANGE])
            listener.updateIsAWIPS         (file_prefix, self.fileSettings[file_prefix][GlanceGUIModel.IS_AWIPS ])
            pass
    
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
                self.fileData[file_prefix].variable = str(newVariableText)
                didUpdate = True
                
                # load the data for this new variable
                self._load_variable_data(file_prefix, str(newVariableText))
                
                # get the variable's attributes
                self._load_variable_attributes (file_prefix, str(newVariableText))
        
        # for convenience hang on to this
        tempVariableName = self.fileData[file_prefix].variable
        
        # update the override selection if needed
        if (newOverrideValue is not None) and (tempVariableName in self.fileData[file_prefix].var_data_cache.keys()) :
            LOG.debug("Setting file " + file_prefix + " override selection to: " + str(newOverrideValue))
            self.fileData[file_prefix].var_data_cache[self.fileData[file_prefix].variable].override_fill_value = newOverrideValue
            didUpdate = True
        
        # update the fill value if needed
        if (newFillValue is not np.nan) and (tempVariableName in self.fileData[file_prefix].var_data_cache.keys()) :
            LOG.debug("Setting file " + file_prefix + " fill value to: " + str(newFillValue))
            self.fileData[file_prefix].var_data_cache[self.fileData[file_prefix].variable].fill_value = newFillValue
            didUpdate = True
        
        # let our data listeners know about any changes
        if didUpdate :
            tempDataObject = self.fileData[file_prefix].var_data_cache[tempVariableName]
            tempAttrsList  = self.fileData[file_prefix].var_attrs_cache[tempVariableName]
            for listener in self.dataListeners :
                listener.fileDataUpdate(file_prefix, self.fileData[file_prefix].file.path, tempVariableName,
                                                     tempDataObject.override_fill_value,   self._select_fill_value(file_prefix),
                                                     str(tempDataObject.data.shape),       attribute_list=tempAttrsList)
    
    def _select_fill_value (self, file_prefix) :
        """
        which fill value should currently be used?
        """
        return self.fileData[file_prefix].var_data_cache[self.fileData[file_prefix].variable].select_fill_value()
    
    def updateSettingsDataSelection (self, newEpsilonValue=np.nan, newImageType=None, newDataForm=None,
                                     newEpsilonPercent=np.nan, newllEpsilon=np.nan,
                                     useSharedRangeForOriginals=None, newColormap=None,
                                     doHideDataFromMismatchedNav=None,
                                     doPlotGeoTiffAsRGB=None) :
        """
        someone has changed one or more of the general settings related data values
        
        Note: if an input value is left at it's default (None or nan)
        then it's assumed that it was not externally changed
        """
        
        didUpdate = False
        
        # update the epsilon if needed
        if (newEpsilonValue is not np.nan) and (newEpsilonValue != self.epsilon) :
            LOG.debug("Setting epsilon to: " + str(newEpsilonValue))
            self.epsilon = newEpsilonValue
            didUpdate = True
        
        # update the epsilon %
        if (newEpsilonPercent is not np.nan) and (newEpsilonPercent != self.epsilonPercent) :
            LOG.debug("Setting epsilon percent to: " + str(newEpsilonPercent))
            self.epsilonPercent = newEpsilonPercent
            didUpdate = True
        
        # update the lon/lat epsilon if needed
        if (newllEpsilon is not np.nan) and (newllEpsilon != self.llEpsilon) :
            LOG.debug("Setting lon/lat epsilon to: " + str(newllEpsilon))
            self.llEpsilon = newllEpsilon
            didUpdate = True
        
        # update the image type if needed
        if (newImageType is not None) and (newImageType != self.imageType) :
            if newImageType in IMAGE_TYPES :
                LOG.debug("Setting image type to: " + newImageType)
                self.imageType = str(newImageType)
                didUpdate = True
        
        # update the colormap if needed
        if (newColormap is not None) and (newColormap != self.colormap) :
            if newColormap in COLORMAP_NAMES :
                LOG.debug("Setting colormap to: " + newColormap)
                self.colormap = str(newColormap)
                didUpdate = True
        
        # update the data form if needed
        if (newDataForm is not None) and (newDataForm != self.dataForm) :
            if newDataForm in DATA_FORMS :
                LOG.debug("Setting data form to: " + newDataForm)
                self.dataForm = str(newDataForm)
                didUpdate = True
        
        # update the shared range settings if needed
        if (useSharedRangeForOriginals is not None) and (useSharedRangeForOriginals != self.useSharedRange) :
            if useSharedRangeForOriginals is True or useSharedRangeForOriginals is False :
                LOG.debug("Setting use shared range for originals to: " + str(useSharedRangeForOriginals))
                self.useSharedRange = useSharedRangeForOriginals
                didUpdate = True
        
        # update the geotiff plotting settings if needed
        if (doPlotGeoTiffAsRGB is not None) and (doPlotGeoTiffAsRGB != self.plotGeoTiffAsRGB) :
            if doPlotGeoTiffAsRGB is True or doPlotGeoTiffAsRGB is False :
                LOG.debug("Setting plot geoTiff data as RGB images to: " + str(doPlotGeoTiffAsRGB))
                self.plotGeoTiffAsRGB = doPlotGeoTiffAsRGB
                didUpdate = True
        
        # update whether or not we'll hide data based on the lon/lat comparison
        if (doHideDataFromMismatchedNav is not None) and (doHideDataFromMismatchedNav != self.hideMismatchNav) :
            if doHideDataFromMismatchedNav is True or doHideDataFromMismatchedNav is False :
                LOG.debug("Setting hide data based on mismatched navigation to: " + str(doHideDataFromMismatchedNav))
                self.hideMismatchNav = doHideDataFromMismatchedNav
                didUpdate = True
        
        # let our data listeners know about any changes
        if didUpdate :
            for listener in self.dataListeners :
                listener.updateEpsilon(self.epsilon)
                listener.updateEpsilonPercent(self.epsilonPercent)
                listener.updateLLEpsilon(self.llEpsilon)
                listener.updateImageTypes(self.imageType)
                listener.updateColormaps(self.colormap)
                listener.updateDataForms(self.dataForm)
                listener.updateUseSharedRange(self.useSharedRange)
                listener.updatePlotGeoTiffAsRGB(self.plotGeoTiffAsRGB)
                listener.updateHideMismatchNav(self.hideMismatchNav)
    
    def updateFileSettings (self, file_prefix, doRestrictRange=None,
                            newRangeMin=np.nan, newRangeMax=np.nan,
                            doCorrectForAWIPS=None) :
        """
        someone has changed one or more of the file specific settings
        
        Note: if an input value is left at it's default (None or nan)
        then it's assumed that it was not externally changed
        """
        
        didUpdate = False
        
        if file_prefix not in self.fileSettings.keys() :
            LOG.warn("Unknown file prefix (" + str(file_prefix) + ") in updateFileSettings.")
            return
        
        # update the range restriction setting if needed
        if (doRestrictRange is not None) and (self.fileSettings[file_prefix][GlanceGUIModel.DO_RANGE] != doRestrictRange) :
            LOG.debug("Setting use range restriction for file " + str(file_prefix) + " to: " + str(doRestrictRange))
            self.fileSettings[file_prefix][GlanceGUIModel.DO_RANGE] = doRestrictRange
            didUpdate = True
        
        if (newRangeMin is not np.nan) and (self.fileSettings[file_prefix][GlanceGUIModel.MIN_RANGE] != newRangeMin) :
            LOG.debug("Setting minimum value for range restriction in file " + str(file_prefix) + " to: " + str(newRangeMin))
            self.fileSettings[file_prefix][GlanceGUIModel.MIN_RANGE] = newRangeMin
            didUpdate = True
        
        if (newRangeMax is not np.nan) and (self.fileSettings[file_prefix][GlanceGUIModel.MAX_RANGE] != newRangeMax) :
            LOG.debug("Setting maximum value for range restriction in file " + str(file_prefix) + " to: " + str(newRangeMax))
            self.fileSettings[file_prefix][GlanceGUIModel.MAX_RANGE] = newRangeMax
            didUpdate = True
        
        if (doCorrectForAWIPS is not None) and (self.fileSettings[file_prefix][GlanceGUIModel.IS_AWIPS] != doCorrectForAWIPS) :
            if (doCorrectForAWIPS is True) or (doCorrectForAWIPS is False) :
                LOG.debug("Setting do AWIPS data correction for file " + str(file_prefix) + " to: " + str(doCorrectForAWIPS))
                self.fileSettings[file_prefix][GlanceGUIModel.IS_AWIPS] = doCorrectForAWIPS
                didUpdate = True
        
        # let our data listeners know about any changes
        if didUpdate :
            self.sendFileSettings(file_prefix)
    
    def updateLonLatSelections (self, file_prefix, new_latitude_name=None, new_longitude_name=None) :
        """
        someone has changed one or more of the file specific longitude/latitude related values
        
        Note: if an input value is left at it's default (None) then it's assumed that it was not externally changed
        """
        
        didUpdate = False
        
        # update the latitude name
        if (new_latitude_name is not None) and (new_latitude_name != self.fileData[file_prefix].latitude) :
            LOG.debug ("Setting latitude name to: " + new_latitude_name)
            self.fileData[file_prefix].latitude = str(new_latitude_name)
            # make sure that this variable is in the cache for use later
            _ = self._load_variable_data (file_prefix, str(new_latitude_name))
            didUpdate = True
        
        # update the longitude name
        if new_longitude_name is not None :
            LOG.debug ("Setting longitude name to: " + new_longitude_name)
            self.fileData[file_prefix].longitude = str(new_longitude_name)
            # make sure that this variable is in the cache for use later
            _ = self._load_variable_data (file_prefix, str(new_longitude_name))
            didUpdate = True
        
        # let our listeners know if we did any updating
        if didUpdate :
            for listener in self.dataListeners :
                listener.updateSelectedLatLon(file_prefix, self.fileData[file_prefix].latitude, self.fileData[file_prefix].longitude)
    
    def getVariableName (self, filePrefix) :
        """
        get the name of the variable loaded for the given file prefix
        or None if no file is loaded
        """
        toReturn = None
        
        if filePrefix in self.fileData.keys() :
            toReturn = self.fileData[filePrefix].variable
        
        return toReturn
    
    def getLongitudeName (self, filePrefix) :
        """
        get the name of the longitude variable selected for the given file prefix
        or None if no file is loaded
        """
        
        toReturn = None
        
        if filePrefix in self.fileData.keys() :
            toReturn = self.fileData[filePrefix].longitude
        
        return toReturn
    
    def getLatitudeName (self, filePrefix) :
        """
        get the name of the latitude variable selected for the given file prefix
        or None if no file is loaded
        """
        
        toReturn = None
        
        if filePrefix in self.fileData.keys() :
            toReturn = self.fileData[filePrefix].latitude
        
        return toReturn
    
    def getVariableData (self, filePrefix, variableName, doCorrections=True) :
        """
        get the data object for the variable of variableName associated with the file prefix
        or None if that variable is not loaded
        
        If doCorrections is True, data filtering for AWIPS and range corrections will be done
        by this function based on the currently selected settings for that file.
        
        Note: this is not a copy, but the original object, so any manipulations
        done to it will be reflected in the model
        """
        toReturn = None
        
        if (filePrefix in self.fileData) and (variableName in self.fileData[filePrefix].var_data_cache) :
            toReturn = self.fileData[filePrefix].var_data_cache[variableName].copy()
            
            # if we should do automatic corrections, do those
            if doCorrections :
                if self.fileSettings[filePrefix][GlanceGUIModel.IS_AWIPS] :
                    fill_mask = toReturn.data == toReturn.fill_value if toReturn.fill_value is not None else np.zeros(toReturn.data.shape, dtype=np.bool)
                    toReturn.data = toReturn.data.astype(np.uint8) # TODO, will setting this break anything?
                    toReturn.data = toReturn.data.astype(np.int32) # make the range larger so we can do comparisons without overflow
                    toReturn.data[fill_mask] = toReturn.fill_value
                
                if self.fileSettings[filePrefix][GlanceGUIModel.DO_RANGE] :
                    if self.fileSettings[filePrefix][GlanceGUIModel.MIN_RANGE] is not None :
                        toReturn.data[toReturn.data < self.fileSettings[filePrefix][GlanceGUIModel.MIN_RANGE]] = toReturn.fill_value
                    if self.fileSettings[filePrefix][GlanceGUIModel.MAX_RANGE] is not None :
                        toReturn.data[toReturn.data > self.fileSettings[filePrefix][GlanceGUIModel.MAX_RANGE]] = toReturn.fill_value
        
        return toReturn
    
    def getUnitsText (self, filePrefix, variableName) :
        """
        get the text describing the units of the variable if the variable exists and that
        attribute exists, otherwise return None
        """
        
        toReturn = None
        
        if (filePrefix in self.fileData) and (self.fileData[filePrefix].file is not None) :
            toReturn = self.fileData[filePrefix].file.file_object.get_attribute(variableName, io.UNITS_CONSTANT)
        
        return toReturn
    
    def getEpsilon (self) :
        """
        get the current value of epsilon
        """
        return self.epsilon
    
    def getEpsilonPercent (self) :
        """
        get the current value of epsilon percent
        """
        return self.epsilonPercent
    
    def getLLEpsilon (self) :
        """
        get the current value of the lon/lat epsilon
        """
        return self.llEpsilon
    
    def getImageType (self) :
        """
        get the text string describing the image type currently selected
        
        the return will correspond to one of the constants from this module:
        
            ORIGINAL_A,
            ORIGINAL_B,
            ABS_DIFF,
            RAW_DIFF,
            HISTOGRAM,
            MISMATCH,
            SCATTER,
            HEX_PLOT
        """
        return self.imageType
    
    def getColormapName (self) :
        """
        get the name of the colormap to use
        
        the return will be on of the constants from gui_constants:
            COLORMAP_NAMES = [CM_RAINBOW, CM_RAINBOW_REV, CM_RAINBOW_DESAT, CM_GRAY, CM_GRAY_REV, CM_SPECTRAL]
        """
        
        return str(self.colormap)
    
    def getIsAWIPS (self, filePrefix) :
        """
        get whether or not the data is in AWIPS format
        """
        
        return self.fileSettings[filePrefix][GlanceGUIModel.IS_AWIPS]
    
    def getDataForm (self) :
        """
        get the text string describing the data form currently selected
        
        the return will correspond to one of the constants from glance.gui_constants:
        
            SIMPLE_2D
            MAPPED_2D
            ONLY_1D
        """
        
        return self.dataForm
    
    def getShouldShowOriginalPlotsInSameRange (self) :
        """
        get whether or not the original plots should be shown in a shared range
        """
        
        return self.useSharedRange
    
    def getDoPlotAsRGB (self, filePrefix) :
        """
        get whether or not multi-channel geotiffs should be treated as RGB or RGBA
        """
        
        isGeoTiff = False
        if ( self.fileData[filePrefix].file is not None ) :
            extension_temp = path.splitext(self.fileData[filePrefix].file.path)[-1]
            isGeoTiff      = (extension_temp == '.tiff') or (extension_temp == '.tif') or (extension_temp == '.tifa')
        
        LOG.debug ("Checking for file " + str(filePrefix) + " geoTiff status... ")
        toReturn = self.plotGeoTiffAsRGB and isGeoTiff
        if toReturn :
            LOG.debug ("file " + str(filePrefix) + " is a geoTiff and should be plotted as an RGB image.")
        
        return toReturn
    
    def getDoHideDataBasedOnMismatchedNavigation (self) :
        """
        get whether or not mapped plots should hide data in places where the
        navigation is more different than llEpsilon
        """
        
        return self.hideMismatchNav
    
    # FUTURE, use this to do more last minute loading?
    def makeSureVariablesAreAvailable (self, filePrefix, listOfVariableNames) :
        """
        given a list of variable names, make sure that they're all loaded
        and available for use.
        
        If the model is able to load all the requested variables (or they're
        already loaded) it will return True, otherwise it will load whichever
        ones it can and return False.
        """
        
        couldLoadAll = True
        
        for varName in listOfVariableNames :
            varDataObj = self._load_variable_data(filePrefix, varName) if varName in self.fileData[filePrefix].ALL_VARIABLES else None
            
            couldLoadAll = couldLoadAll and (varDataObj is not None)
        
        return couldLoadAll
    
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
