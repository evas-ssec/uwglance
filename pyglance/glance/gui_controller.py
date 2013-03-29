#!/usr/bin/env python
# encoding: utf-8
"""
The controller portion of the glance GUI code. 

Created by evas Oct 2011.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
"""

import sys, os.path, logging

from PyQt4 import QtGui

import glance.gui_view          as gui_view
import glance.gui_model         as gui_model
import glance.gui_statsprovider as gui_stats
import glance.gui_figuremanager as gui_figs

from glance.data import IncompatableDataObjects

LOG = logging.getLogger(__name__)

"""
The controller is resposible for facilitating communication between the GUI and the model.
It listens to both the userUpdates from the GUI and the errors from the model and passes that
information on (or handles it itself) as needed.
"""

class GlanceGUIController (object) :
    """
    This class is responsible for arranging communication in the Glance GUI.
    It includes:
    
    self.view  - the view object (see glance.gui_view)
    self.model - the model object (see glance.gui_model)
    self.stats - the stats provider object (see glance.gui_statsprovider)
    self.figs  - the figure manager object (see glance.gui_figuremanager)
    self.qtApp - an application object, used to start QT
    """
    
    def __init__ (self, version_string) :
        """
        build the various objects required for the GUI and be ready to launch them
        """
        
        # create the various other entities
        self.qtApp = QtGui.QApplication(sys.argv)
        self.view  = gui_view.GlanceGUIView(version_string)
        self.model = gui_model.GlanceGUIModel()
        self.stats = gui_stats.GlanceGUIStats(self.model)
        self.figs  = gui_figs.GlanceGUIFigures(self.model)
        
        # set things up to talk to each other
        self.model.registerErrorHandler(self)
        self.model.registerDataListener(self.view)
        self.stats.registerErrorHandler(self)
        self.stats.registerStatsListener(self.view)
        self.figs.registerErrorHandler(self)
        self.view.registerUserUpdateListener(self)
        
        # force the initial info load from the model
        # (this is stuff like default epsilon and image types)
        self.model.sendGeneralSettingsData()
    
    def launch_gui (self) :
        """
        start up the GUI
        """
        self.view.show()
        return(self.qtApp.exec_()) # this is needed to keep the gui from exiting prematurely
    
    ################# methods to handle user input reporting #################
    
    # a method that comes from the view
    def newFileSelected(self, file_prefix, new_file_path) :
        """
        a new file has been selected by the user
        """
        
        try :
            self.model.loadNewFile(file_prefix, new_file_path)
        except (gui_model.UnableToReadFile, ValueError) as utrf :
            self.handleWarning(str(utrf))
    
    def userSelectedVariable (self, file_prefix, newSelection) :
        """
        the user selected a new variable
        """
        
        try :
            self.model.updateFileDataSelection(file_prefix, newVariableText=newSelection)
        except ValueError as ve :
            self.handleWarning(str(ve))
    
    def userChangedOverload (self, file_prefix, new_override_value) :
        """
        the user checked or unchecked the override box
        """
        
        self.model.updateFileDataSelection(file_prefix, newOverrideValue=new_override_value)
    
    def userChangedFillValue (self, file_prefix, new_fill_value) :
        """
        the user has entered a new fill value
        """
        
        self.model.updateFileDataSelection(file_prefix, newFillValue=new_fill_value)
    
    def userChangedEpsilon (self, new_epsilon) :
        """
        the user has entered a new epsilon value
        """
        
        self.model.updateSettingsDataSelection(newEpsilonValue=new_epsilon)
    
    def userChangedEpsilonPercent (self, new_epsilon_percent) :
        """
        the user has entered a new epsilon percent value
        """
        
        self.model.updateSettingsDataSelection(newEpsilonPercent=new_epsilon_percent)
    
    def userChangedLLEpsilon (self, new_lonlat_epsilon) :
        """
        the user has entered a new lon/lat epsilon
        """
        
        self.model.updateSettingsDataSelection(newllEpsilon=new_lonlat_epsilon)
    
    def userSelectedLongitude (self, file_prefix, newSelection) :
        """
        the user selected a new longitude variable
        """
        
        try:
            self.model.updateLonLatSelections(file_prefix, new_longitude_name=newSelection)
        except ValueError as ve :
            self.handleWarning(str(ve))
    
    def userSelectedLatitude (self, file_prefix, newSelection) :
        """
        the user selected a new latitude variable
        """
        
        try :
            self.model.updateLonLatSelections(file_prefix, new_latitude_name=newSelection)
        except ValueError as ve :
            self.handleWarning(str(ve))
    
    def userSelectedImageType (self, new_image_type) :
        """
        the user has selected a new image type
        """
        
        self.model.updateSettingsDataSelection(newImageType=new_image_type)
    
    def userSelectedColormap (self, new_colormap) :
        """
        the user has selected a new colormap
        """
        
        self.model.updateSettingsDataSelection(newColormap=new_colormap)
    
    def userSelectedDataForm (self, new_data_form) :
        """
        the user has selected a new data form
        """
        
        self.model.updateSettingsDataSelection(newDataForm=new_data_form)
    
    def userToggledSharedRange(self, should_use_shared_range) :
        """
        the user has toggled whether or not the original data should use a shared range
        """
        
        self.model.updateSettingsDataSelection(useSharedRangeForOriginals=should_use_shared_range)
    
    def userToggledRestrictRange(self, file_prefix, should_restrict_range) :
        """
        the user has toggled whether or not to restrict the data to a fixed range
        """
        
        self.model.updateFileSettings(file_prefix, doRestrictRange=should_restrict_range)
    
    def userChangedRangeMin(self, file_prefix, new_range_min) :
        """
        the user changed the minimum of the acceptable data range
        """
        
        self.model.updateFileSettings(file_prefix, newRangeMin=new_range_min)
    
    def userChangedRangeMax(self, file_prefix, new_range_max) :
        """
        the user changed the maximum of the acceptable data range
        """
        
        self.model.updateFileSettings(file_prefix, newRangeMax=new_range_max)
    
    def userToggledIsAWIPS(self, file_prefix, data_is_AWIPS) :
        """
        the user has toggled whether or not the file should be treated as AWIPS formatted data
        """
        
        self.model.updateFileSettings(file_prefix, doCorrectForAWIPS=data_is_AWIPS)
    
    def userRequestsStats (self) :
        """
        the user has asked for stats information
        """
        
        try :
            self.stats.sendStatsInfo()
        except IncompatableDataObjects as ido :
            self.handleWarning(str(ido))
    
    def userRequestsPlot (self) :
        """
        the user has asked for a plot
        """
        
        try :
            self.figs.spawnPlot()
        except (IncompatableDataObjects, ValueError) as idove :
            self.handleWarning(str(idove))
            #raise
    
    ################# end of methods to handle user input reporting #################
    
    # handle warnings reported from the model
    def handleWarning (self, warningMessage) :
        """
        this warning needs to go to the user
        """
        
        LOG.warn(warningMessage)
        self.view.showWarning(warningMessage)
    

if __name__=='__main__':
    import doctest
    doctest.testmod()
