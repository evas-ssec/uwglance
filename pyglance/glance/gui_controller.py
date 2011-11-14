#!/usr/bin/env python
# encoding: utf-8
"""
The controller portion of the glance GUI code. 

Created by evas Oct 2011.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
"""

import sys, os.path, logging

PYQT4_HAX = '/sw/lib/qt4-mac/lib/python2.6/site-packages'
if os.path.isdir(PYQT4_HAX):
    sys.path.append(PYQT4_HAX)

from PyQt4 import QtGui

import glance.gui_view  as gui_view
import glance.gui_model as gui_model

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
        
        # set things up to talk to each other
        self.model.registerErrorHandler(self)
        self.model.registerDataListener(self.view)
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
        
        self.model.loadNewFile(file_prefix, new_file_path)
    
    def userSelectedVariable (self, file_prefix, newSelection) :
        """
        the user selected a new variable
        """
        
        self.model.updateFileDataSelection(file_prefix, newVariableText=newSelection)
    
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
    
    def userSelectedImageType (self, new_image_type) :
        """
        the user has selected a new image type
        """
        
        self.model.updateSettingsDataSelection(newImageType=new_image_type)
    
    def userRequestsPlot (self) :
        """
        the user has asked for a plot
        """
        
        self.model.spawnPlotWithCurrentInfo()
    
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
