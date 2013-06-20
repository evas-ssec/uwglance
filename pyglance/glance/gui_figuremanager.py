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

from mpl_toolkits.basemap import Basemap

import logging
import numpy as np

import glance.data      as dataobjects
import glance.figures   as figures
import glance.gui_model as model
from   glance.gui_constants import *
from   glance.plotcreatefns import select_projection

LOG = logging.getLogger(__name__)

#temp_dict = {'blue': [(0.0, 0.75, 0.75), (0.11, 0.99955436720142599, 0.99955436720142599), (0.34000000000000002, 0.99810246679316883, 0.99810246679316883), (0.34999999999999998, 0.98545224541429477, 0.98545224541429477), (0.375, 0.94117647058823528, 0.94117647058823528), (0.64000000000000001, 0.51739405439595187, 0.51739405439595187), (0.65000000000000002, 0.5, 0.5), (0.66000000000000003, 0.5, 0.5), (0.89000000000000001, 0.5, 0.5), (0.91000000000000003, 0.5, 0.5), (1.0, 0.5, 0.5)], 'green': [(0.0, 0.5, 0.5), (0.11, 0.5, 0.5), (0.125, 0.50098039215686274, 0.50098039215686274), (0.34000000000000002, 0.93235294117647061, 0.93235294117647061), (0.34999999999999998, 0.94803921568627447, 0.94803921568627447), (0.375, 1.0, 1.0), (0.64000000000000001, 1.0, 1.0), (0.65000000000000002, 0.97966594045025435, 0.97966594045025435), (0.66000000000000003, 0.96514161220043593, 0.96514161220043593), (0.89000000000000001, 0.53667392883079168, 0.53667392883079168), (0.91000000000000003, 0.50036310820624552, 0.50036310820624552), (1.0, 0.5, 0.5)], 'red': [(0.0, 0.5, 0.5), (0.11, 0.5, 0.5), (0.125, 0.5, 0.5), (0.34000000000000002, 0.5, 0.5), (0.34999999999999998, 0.5, 0.5), (0.375, 0.54269449715370022, 0.54269449715370022), (0.64000000000000001, 0.96647691334598351, 0.96647691334598351), (0.65000000000000002, 0.98545224541429466, 0.98545224541429466), (0.66000000000000003, 0.99810246679316883, 0.99810246679316883), (0.89000000000000001, 0.99955436720142621, 0.99955436720142621), (0.91000000000000003, 0.9549910873440286, 0.9549910873440286), (1.0, 0.75, 0.75)]}
temp_dict = {'blue': [(0.0, 0.58333333333333326, 0.58333333333333326), (0.11, 0.91607248960190135, 0.91607248960190135), (0.125, 0.91666666666666663, 0.91666666666666663), (0.34000000000000002, 0.91413662239089188, 0.91413662239089188), (0.34999999999999998, 0.89726966055239299, 0.89726966055239299), (0.375, 0.83823529411764708, 0.83823529411764708), (0.64000000000000001, 0.27319207252793593, 0.27319207252793593), (0.65000000000000002, 0.25, 0.25), (0.66000000000000003, 0.25, 0.25), (0.89000000000000001, 0.25, 0.25), (0.91000000000000003, 0.25, 0.25), (1.0, 0.25, 0.25)], 'green': [(0.0, 0.25, 0.25), (0.11, 0.25, 0.25), (0.125, 0.25130718954248366, 0.25130718954248366), (0.34000000000000002, 0.82647058823529418, 0.82647058823529418), (0.34999999999999998, 0.84738562091503267, 0.84738562091503267), (0.375, 0.91666666666666663, 0.91666666666666663), (0.64000000000000001, 0.91666666666666663, 0.91666666666666663), (0.65000000000000002, 0.88955458726700576, 0.88955458726700576), (0.66000000000000003, 0.87018881626724787, 0.87018881626724787), (0.89000000000000001, 0.29889857177438889, 0.29889857177438889), (0.91000000000000003, 0.25048414427499405, 0.25048414427499405), (1.0, 0.25, 0.25)], 'red': [(0.0, 0.25, 0.25), (0.11, 0.25, 0.25), (0.125, 0.25, 0.25), (0.34000000000000002, 0.25, 0.25), (0.34999999999999998, 0.25, 0.25), (0.375, 0.30692599620493355, 0.30692599620493355), (0.64000000000000001, 0.87196921779464465, 0.87196921779464465), (0.65000000000000002, 0.89726966055239288, 0.89726966055239288), (0.66000000000000003, 0.91413662239089177, 0.91413662239089177), (0.89000000000000001, 0.91607248960190157, 0.91607248960190157), (0.91000000000000003, 0.85665478312537158, 0.85665478312537158), (1.0, 0.58333333333333326, 0.58333333333333326)]}
DESAT_MAP = matplotlib.colors.LinearSegmentedColormap('colormap', temp_dict, 1024)

# colormaps that are available in the GUI
# if this changes the list of colormap names in the constants module needs to be kept up
AVAILABLE_COLORMAPS = {CM_RAINBOW:       cm.jet,
                       CM_RAINBOW_REV:   cm.jet_r,
                       CM_RAINBOW_DESAT: DESAT_MAP,
                       CM_GRAY:          cm.bone,
                       CM_GRAY_REV:      cm.bone_r,
                       CM_SPECTRAL:      cm.spectral}

# whether or not the plot can be drawn on a map
CAN_BE_MAPPED = {
                    ORIGINAL_A      : True,
                    ORIGINAL_B      : True,
                    ABS_DIFF        : True,
                    RAW_DIFF        : True,
                    HISTOGRAM_A     : False,
                    HISTOGRAM_B     : False,
                    HISTOGRAM       : False,
                    MISMATCH        : True,
                    SCATTER         : False,
                    HEX_PLOT        : False,
                }

# which data sets the plot needs
NEEDED_DATA_PER_PLOT = \
                {
                    ORIGINAL_A      : set([A_CONST]),
                    ORIGINAL_B      : set([         B_CONST]),
                    ABS_DIFF        : set([A_CONST, B_CONST]),
                    RAW_DIFF        : set([A_CONST, B_CONST]),
                    HISTOGRAM_A     : set([A_CONST]),
                    HISTOGRAM_B     : set([         B_CONST]),
                    HISTOGRAM       : set([A_CONST, B_CONST]),
                    MISMATCH        : set([A_CONST, B_CONST]),
                    SCATTER         : set([A_CONST, B_CONST]),
                    HEX_PLOT        : set([A_CONST, B_CONST]),
                }

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
    
    def _getVariableInformation (self, filePrefix, variableName=None, doCorrections=True) :
        """
        Pull the name, data, and units for the variable currently selected in the given file prefix
        """
        varNameToUse = variableName
        if varNameToUse is None :
            varNameToUse = self.dataModel.getVariableName(filePrefix) # get the currently selected variable
        
        dataObject           = self.dataModel.getVariableData(filePrefix, varNameToUse,  doCorrections=doCorrections)
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
             ( filePrefix in NEEDED_DATA_PER_PLOT[imageType] ) ) :
            
            shouldUseRGBVersion = self.dataModel.getDoPlotAsRGB(filePrefix) and ( (imageType == ORIGINAL_A) or (imageType == ORIGINAL_B) )
            varName, dataObject, unitsText = self._getVariableInformation(filePrefix) if not shouldUseRGBVersion else self._makeRGBdata(filePrefix)
        
        return varName, dataObject, unitsText
    
    def _makeRGBdata (self, filePrefix) :
        """
        build an RGB or RGBA version of the data
        """
        
        # get the red, green, and blue data
        canGetData = self.dataModel.makeSureVariablesAreAvailable(filePrefix, [RED_VAR_NAME, GREEN_VAR_NAME, BLUE_VAR_NAME])
        if not canGetData : # if the basic rgb data doesn't exist, stop now
            "", None, ""
        _, rDataObj, _ = self._getVariableInformation(filePrefix, variableName=RED_VAR_NAME,   doCorrections=False)
        _, gDataObj, _ = self._getVariableInformation(filePrefix, variableName=GREEN_VAR_NAME, doCorrections=False)
        _, bDataObj, _ = self._getVariableInformation(filePrefix, variableName=BLUE_VAR_NAME,  doCorrections=False)
        
        # if possible get alpha data
        _ = self.dataModel.makeSureVariablesAreAvailable(filePrefix, [ALPHA_VAR_NAME]) # we need to make sure the model loads the data, but it's optional
        _, aDataObj, _ = self._getVariableInformation(filePrefix, variableName=ALPHA_VAR_NAME, doCorrections=False)
        
        # build the finished rgb set
        rawData = [rDataObj.data, gDataObj.data, bDataObj.data] if aDataObj is None else [rDataObj.data, gDataObj.data, bDataObj.data, aDataObj.data]
        rawData = np.rot90(np.fliplr(np.transpose(np.array(rawData))))
        # now that the data is in the right shape/orientation make the data object
        newDataObj = dataobjects.DataObject(rawData, fillValue=rDataObj.fill_value) # TODO, need to fix the fill values if they differ
        newDataObj.self_analysis()
        
        # return varName, dataObject, unitsText
        return "rgb data", newDataObj, ""
    
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
    
    def _load_and_analyse_lonlat (self, listOfFilePrefixes=[A_CONST, B_CONST], lonNames=None, latNames=None, stopIfComparisonFails=False) :
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
            currentLonObj.self_analysis()
            currentLatObj.self_analysis()
            
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
            
            # we can't use longitude and latitude that don't match in size
            if currentLonObj.data.shape != currentLatObj.data.shape :
                raise ValueError ("Longitude and Latitude for file " + filePrefix + " are different shapes." +
                                  "\nCannot match differently shaped navigation data.")
            
            # compare this file to whatever other data we have
            for filePrefixToCompare in lonlatData.keys() :
                lonToCompare, latToCompare = lonlatData[filePrefixToCompare]
                # make sure the files are the same shape
                if (currentLonObj.data.shape != lonToCompare.data.shape) :
                    raise ValueError ("Navigation data for file " + filePrefix +
                                      " is a different shape than that for file " + filePrefixToCompare + "." +
                                      "\nCannot match differently shaped navigation data.")
            
            # add this data to the list of lonlat data
            lonlatData[filePrefix] = [currentLonObj, currentLatObj]
        
        # return longitude and latitude information and the shared ranges
        return lonlatData, lonRange, latRange
    
    def _load_lonlat (self, filePrefix, lonName, latName) :
        """
        load the longitude and latitude information for the file and determine the ranges
        present in both
        """
        
        _, lonObject, _ = self._getVariableInformation(filePrefix, lonName, doCorrections=False)
        _, latObject, _ = self._getVariableInformation(filePrefix, latName, doCorrections=False)
        
        lonRange = [lonObject.get_min(), lonObject.get_max()]
        latRange = [latObject.get_min(), latObject.get_max()]
        
        return lonObject, latObject, lonRange, latRange
    
    def _find_common_lonlat (self, lonlatData, doUnion=False) :
        """
        given lonlatData like that created by _load_and_analyse_lonlat
        find a common set of longitude and latitude
        
        If doUnion is True, create a set that contains valid
        longitudes and latitudes in as many places as possible.
        Navigation data will be chosen preferentially based on
        the sorting order of the keys in lonlatData.
        If doUnion is False, the intersection of the data will
        be produced instead (using the first data set by key
        order and masking by data placement in later sets).
        """
        
        commonLon = None
        commonLat = None
        validMask = None
        
        # look through each of the possible data sets
        for file_prefix in sorted(lonlatData.keys()) :
            tempLonObj, tempLatObj = lonlatData[file_prefix]
            if commonLon is None :
                commonLon = tempLonObj.copy()
                commonLat = tempLatObj.copy()
                commonLon.self_analysis()
                commonLat.self_analysis()
                validMask = commonLon.masks.valid_mask & commonLat.masks.valid_mask
            else :
                tempLonObj.self_analysis()
                tempLatObj.self_analysis()
                if doUnion :
                    newValid = (tempLatObj.masks.valid_mask & tempLonObj.masks.valid_mask) & ~ validMask
                    commonLon.data[newValid] = tempLonObj.data[newValid]
                    commonLat.data[newValid] = tempLatObj.data[newValid]
                    validMask |= newValid
                else:
                    newInvalid = ~(tempLatObj.masks.valid_mask & tempLonObj.masks.valid_mask) & validMask
                    commonLon.data[newInvalid] = commonLon.fill_value
                    commonLat.data[newInvalid] = commonLat.fill_value
                    validMask &= ~newInvalid
        
        # since we changed the data, rebuild the internal analysis
        commonLat.self_analysis(re_do_analysis=True)
        commonLon.self_analysis(re_do_analysis=True)
        
        LOG.debug("common lon/lat validMask.shape: " + str(validMask.shape))
        LOG.debug("common lon/lat sum(validMask):  " + str(sum(validMask)))
        
        return commonLon, commonLat, validMask
    
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
        aVarName, aDataObject, aUnitsText = self._getVariableInfoSmart(A_CONST, imageType)
        bVarName, bDataObject, bUnitsText = self._getVariableInfoSmart(B_CONST, imageType)
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
        lonlatData     = None
        basemapObject  = None
        lonlatWarnings = ""
        if ((dataForm == MAPPED_2D) and CAN_BE_MAPPED[imageType]) :
            
            # get the longitude and latitude information for the files, as needed
            dataNeeded = list(NEEDED_DATA_PER_PLOT[imageType]) # this is naturally a set, use a list here
            lonNames = { }
            latNames = { }
            for file_const in dataNeeded :
                lonNames[file_const] = self.dataModel.getLongitudeName(file_const)
                latNames[file_const] = self.dataModel.getLatitudeName (file_const)
            lonlatData, lonRange, latRange = self._load_and_analyse_lonlat(listOfFilePrefixes=dataNeeded,
                                                                           lonNames=lonNames, latNames=latNames)
            
            # double check that lon/lat are compatable with the data
            if (aDataObject is not None) and (A_CONST in dataNeeded) :
                if lonlatData[A_CONST][0].data.shape != aDataObject.data.shape :
                    raise ValueError("Unable to use selected navigation variables for file " + A_CONST +
                                     "\nbecause they differ in size from the selected data variable for that file.")
            if (bDataObject is not None) and (B_CONST in dataNeeded) :
                if lonlatData[B_CONST][0].data.shape != bDataObject.data.shape :
                    raise ValueError("Unable to use selected navigation variables for file " + B_CONST +
                                     "\nbecause they differ in size from the selected data variable for that file.")
            # FUTURE if there were ever more data sets, they'd need to be checked individually or make this more general?
            
            # build basemap and axes,
            # FUTURE, don't hard code so much of this stuff, let the projection and possibly others be selected
            # FUTURE, some of this is in graphics.py, but needs to be refactored so I can call it in a different way
            # FUTURE (may go with the axis finding changes from Graeme)
            boundingAxes  = [lonRange[0], lonRange[1], latRange[0], latRange[1]]
            projToUse     = select_projection(boundingAxes)
            LOG.debug("Selecting projection: " + projToUse)
            midLat        = (latRange[0] + latRange[1]) / 2.0 # this will fail horribly where we cross discontinious lines
            midLon        = (lonRange[0] + lonRange[1]) / 2.0 # this will fail horribly where we cross discontinious lines
            if projToUse is 'ortho' :
                basemapObject = Basemap(lat_0=midLat, lon_0=midLon, resolution='i', area_thresh=10000., projection=projToUse)
            else :
                basemapObject = Basemap(llcrnrlon=lonRange[0], urcrnrlon=lonRange[1],
                                        llcrnrlat=latRange[0], urcrnrlat=latRange[1],
                                        lat_1=midLat, lon_0=midLon,
                                        resolution='i', area_thresh=10000., projection=projToUse)
            
            # do a rough comparison of the longitude and latitude
            if (aDataObject is not None) and (bDataObject is not None) :
                llEpsilon = self.dataModel.getLLEpsilon()
                lonDiffInfo = dataobjects.DiffInfoObject(lonlatData[A_CONST][0],
                                                         lonlatData[B_CONST][0],
                                                         epsilonValue=llEpsilon)
                latDiffInfo = dataobjects.DiffInfoObject(lonlatData[A_CONST][1],
                                                         lonlatData[B_CONST][1],
                                                         epsilonValue=llEpsilon)
                validA = lonlatData[A_CONST][0].masks.valid_mask & lonlatData[A_CONST][1].masks.valid_mask
                validB = lonlatData[B_CONST][0].masks.valid_mask & lonlatData[B_CONST][1].masks.valid_mask
                
                if sum(validA ^ validB) > 0 :
                    lonlatWarnings += "Valid areas in the two files do not match.\n"
                    lonlatWarnings += ("File " + A_CONST + " contains " + str(sum(validA & ~ validB)) +
                                       " points which are not valid in file " + B_CONST + ".\n")
                    lonlatWarnings += ("File " + B_CONST + " contains " + str(sum(validB & ~ validA)) +
                                       " points which are not valid in file " + A_CONST + ".\n")
                
                if sum(lonDiffInfo.diff_data_object.masks.outside_epsilon_mask) > 0 :
                    lonlatWarnings += (str(sum(lonDiffInfo.diff_data_object.masks.outside_epsilon_mask)) +
                                       " longitude points differed by more than the epsilon of " +
                                       str(llEpsilon) + " between the two files.\n")
                if sum(latDiffInfo.diff_data_object.masks.outside_epsilon_mask) > 0 :
                    lonlatWarnings += (str(sum(latDiffInfo.diff_data_object.masks.outside_epsilon_mask)) +
                                       " latitude points differed by more than the epsilon of " +
                                       str(llEpsilon) + " between the two files.\n")
        
        LOG.info("Spawning plot window: " + imageType)
        
        plt.ion() # make sure interactive plotting is on
        
        # create whichever type of plot was asked for
        
        if (imageType == ORIGINAL_A) or (imageType == ORIGINAL_B) :
            
            # sort out some values based on which of the data sets we're showing
            data_object_to_use = aDataObject if (imageType == ORIGINAL_A) else bDataObject
            var_name_to_use    = aVarName    if (imageType == ORIGINAL_A) else bVarName
            file_char_to_use   = A_CONST     if (imageType == ORIGINAL_A) else B_CONST
            units_text_to_use  = aUnitsText  if (imageType == ORIGINAL_A) else bUnitsText
            oneD_color_to_use  = 'b'         if (imageType == ORIGINAL_A) else 'c'
            
            plotAsRGB          = self.dataModel.getDoPlotAsRGB(A_CONST if imageType == ORIGINAL_A else B_CONST)
            
            # if the data doesn't exist, we can't make this plot
            if data_object_to_use is None :
                
                raise ValueError(NO_DATA_MESSAGE)
            
            if dataForm == SIMPLE_2D :
                if plotAsRGB :
                    figures.create_raw_image_plot(data_object_to_use.data, "RGB image in File " + file_char_to_use)
                else :
                    tempFigure = figures.create_simple_figure(data_object_to_use.data, var_name_to_use + "\nin File " + file_char_to_use,
                                                              invalidMask=~data_object_to_use.masks.valid_mask, colorMap=colorMapToUse,
                                                              colorbarLimits=rangeInfo, units=units_text_to_use)
                
            elif dataForm == MAPPED_2D :
                tempLonObj = lonlatData[file_char_to_use][0]
                tempLatObj = lonlatData[file_char_to_use][1]
                tempValid  = data_object_to_use.masks.valid_mask
                tempValid  &= tempLonObj.masks.valid_mask
                tempValid  &= tempLatObj.masks.valid_mask
                tempFigure = figures.create_mapped_figure(data_object_to_use.data,
                                                          tempLatObj.data, tempLonObj.data,
                                                          basemapObject, boundingAxes, 
                                                          var_name_to_use + "\nin File " + file_char_to_use,
                                                          invalidMask=~tempValid, colorMap=colorMapToUse,
                                                          units=units_text_to_use)
                
            elif dataForm == ONLY_1D :
                temp = [(data_object_to_use.data, ~data_object_to_use.masks.valid_mask, oneD_color_to_use, None, None, None)]
                tempFigure = figures.create_line_plot_figure(temp, var_name_to_use + "\n in File " + file_char_to_use)
            else :
                raise ValueError(UNKNOWN_DATA_FORM)
            
        elif (imageType == HISTOGRAM_A) or (imageType == HISTOGRAM_B) :
            
            # Note: histograms don't care about data format requested, they are histogram formatted
            
            # select the things that are file A or B specific
            file_desc_to_use   = A_CONST     if (imageType == HISTOGRAM_A) else B_CONST
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
                    
                    tempLonObj, tempLatObj, tempValid = self._find_common_lonlat(lonlatData)
                    tempValid &= diffData.diff_data_object.masks.valid_mask
                    tempFigure = figures.create_mapped_figure(dataToUse,
                                                              tempLatObj.data, tempLonObj.data,
                                                              basemapObject, boundingAxes, 
                                                              titlePrefix + aVarName,
                                                              invalidMask=~tempValid, colorMap=colorMapToUse,
                                                              units=aUnitsText)
                    
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
                elif dataForm == MAPPED_2D :
                    
                    tempLonObj, tempLatObj, tempValid = self._find_common_lonlat(lonlatData, doUnion=True)
                    tempValid &= (aDataObject.masks.valid_mask | bDataObject.masks.valid_mask)
                    tempData = aDataObject.copy()
                    tempMask = bDataObject.masks.valid_mask & ~aDataObject.masks.valid_mask
                    tempData.data[tempMask] = bDataObject.data[tempMask]
                    tempFigure = figures.create_mapped_figure(tempData.data,
                                                              tempLatObj.data, tempLonObj.data,
                                                              basemapObject, boundingAxes, 
                                                              "Areas of mismatch data\nin " + aVarName,
                                                              invalidMask=~tempValid,
                                                              tagData=mismatchMask,
                                                              colorMap=figures.MEDIUM_GRAY_COLOR_MAP,
                                                              units=aUnitsText)
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
        
        if lonlatWarnings != "" :
            raise ValueError(lonlatWarnings)
