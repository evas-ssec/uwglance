#!/usr/bin/env python
# encoding: utf-8
"""
The view portion of the glance GUI code. 

Created by evas Oct 2011.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
"""

import sys, os.path, logging

PYQT4_HAX = '/sw/lib/qt4-mac/lib/python2.6/site-packages'
if os.path.isdir(PYQT4_HAX):
    sys.path.append(PYQT4_HAX)

from PyQt4 import QtGui, QtCore

LOG = logging.getLogger(__name__)

"""
The GUI view is pretty simple.
All it needs to do is send messages when data is changed and
accept UI updates (of both values and lock state) from the model and controller.

The GUI's most complex responsibility is keeping the user from entering
garbage data into the UI.
"""

class GlanceGUIView (QtGui.QWidget) :
    """
    The main view object that will create the gui that's visible to the user.
    
    It includes:
    
    self.userUpdateListeners    - objects that want to be notified when the user
                                  changes data in the gui
    self.lastFilePath           - the last path that was loaded on this program run,
                                  starts as './'
    
    self.epsilonWidget          - the widget that handles the input of epsilon
    self.imageSelectionDropDown - the drop down that lets the user select the type
                                  of image to be created
    self.displayButton          - the button that lets the user create a plot
    self.widgetInfo             - a dictionary of the widgets that need to be
                                  classifed by file; these are indexed in the form
                                  self.widgetInfo[file_prefix][widget_name] and
                                  include the following widget names:
            'path' - the file path text display for that file
            'load' - the load button for that file
            'variable' - the variable drop down selector for that file
            'dims' - the label to display dimensions for that file
            'attrs' - the table to display attributes for that file
            'override' - the override check box for that file
            'fillValue' - the fill value for that file
    """
    
    def __init__ (self, versionString, parent=None) :
        """
        build the various Qt controls that make up this application
        """
        
        QtGui.QWidget.__init__(self,parent)
        
        # set our title with the version string
        self.setWindowTitle(versionString)
        
        # setup the rest of our window
        self._setup_qt_controls()
        
        # this will represent those who want to be notified
        # when the user changes things
        self.userUpdateListeners = [ ]
        
        # we will use this to remember were the user wanted to load files from last
        # TODO, can we remember this between program runs? something like preferences?
        self.lastFilePath = './'
        
        # hang on to stats windows so they don't vanish
        self.statsWindows = { }
        self.statsCounter = 1
    
    def _add_file_related_controls (self, file_prefix, grid_layout, currentRow) :
        """
        add a set of file related controls to the widget, using the
        given name prefix and grid layout and starting on "currentRow"
        
        return the index of the next empty "currentRow" after we're finished adding controls
        """
        
        # where we'll hang onto the controls we need to access later
        self.widgetInfo[file_prefix] = { }
        
        # set up the file loading control
        grid_layout.addWidget(QtGui.QLabel("File " + file_prefix + ":"), currentRow, 0)
        filePath   = QtGui.QLineEdit()
        filePath.setReadOnly(True) # this is mostly for displaying the file path, the load button selects it
        self.widgetInfo[file_prefix]['path'] = filePath
        grid_layout.addWidget(filePath, currentRow, 1, 1, 3)
        loadButton = QtGui.QPushButton("Load")
        # set some tooltip text
        loadButton.setToolTip("Load a file: glance can handle NetCDF, HDF4, HDF5, and AERI files")
        # connect the button to an action
        if   file_prefix is "A" :
            loadButton.clicked.connect(self.clickedALoad)
        elif file_prefix is "B" :
            loadButton.clicked.connect(self.clickedBLoad)
        self.widgetInfo[file_prefix]['load'] = loadButton
        grid_layout.addWidget(loadButton, currentRow, 4)
        
        currentRow = currentRow + 1
        
        # set up the drop down for the variable select
        grid_layout.addWidget(QtGui.QLabel("variable name:"), currentRow, 1)
        variableSelection = QtGui.QComboBox()
        variableSelection.setDisabled(True)
        if   file_prefix is "A" :
            variableSelection.activated.connect(self.selectedAVariable)
        elif file_prefix is "B" :
            variableSelection.activated.connect(self.selectedBVariable)
        self.widgetInfo[file_prefix]['variable'] = variableSelection
        grid_layout.addWidget(variableSelection, currentRow, 2, 1, 3)
        
        currentRow = currentRow + 1
        
        # set up a label to display the variable dimension information
        tempShapeLabel = QtGui.QLabel("data shape:")
        tempShapeLabel.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        grid_layout.addWidget(tempShapeLabel, currentRow, 1)
        dimensionsLabel = QtGui.QLabel(" ")
        self.widgetInfo[file_prefix]['dims'] = dimensionsLabel
        grid_layout.addWidget(dimensionsLabel, currentRow, 2, 1, 3)
        
        currentRow = currentRow + 1
        
        # set up a table to display variable attribute information
        tempAttributesTable = QtGui.QTableWidget()
        tempAttributesTable.setColumnCount(2)
        # set up the table headers
        tempAttributesTable.setHorizontalHeaderLabels(["variable attribute", "value"])
        tempAttributesTable.horizontalHeader().setResizeMode(0, QtGui.QHeaderView.Stretch)
        tempAttributesTable.horizontalHeader().setResizeMode(1, QtGui.QHeaderView.Stretch)
        tempAttributesTable.verticalHeader().hide()
        # save the widget and put it in our layout
        self.widgetInfo[file_prefix]['attrs'] = tempAttributesTable
        grid_layout.addWidget(tempAttributesTable, currentRow, 1, 1, 4)
        
        currentRow = currentRow + 1
        
        # set up a check box to override the fill value loaded from the file
        overrideFillButton = QtGui.QCheckBox("override fill value")
        overrideFillButton.setDisabled(True)
        if   file_prefix is "A" :
            overrideFillButton.stateChanged.connect(self.toggledAOverride)
        elif file_prefix is "B" :
            overrideFillButton.stateChanged.connect(self.toggledBOverride)
        self.widgetInfo[file_prefix]['override'] = overrideFillButton
        grid_layout.addWidget(overrideFillButton, currentRow, 1)
        
        # now set up the input of the fill value that will be used
        grid_layout.addWidget(QtGui.QLabel("fill value:"), currentRow+1, 1)
        fillValue = QtGui.QLineEdit()
        tempValidator = QtGui.QDoubleValidator(fillValue)
        tempValidator.setNotation(QtGui.QDoubleValidator.StandardNotation)
        fillValue.setValidator(tempValidator)
        fillValue.setDisabled(True)
        if   file_prefix is "A" :
            fillValue.editingFinished.connect(self.fillValueEditedA)
        elif file_prefix is "B" :
            fillValue.editingFinished.connect(self.fillValueEditedB)
        self.widgetInfo[file_prefix]['fillValue'] = fillValue
        grid_layout.addWidget(fillValue, currentRow+1, 2, 1, 3)
        
        currentRow = currentRow + 2
        
        return currentRow
    
    def _setup_qt_controls (self) :
        """
        built the basic input boxes / labels / buttons / etc and lay them out on the window
        """
        
        # create the layout and set up some of the overall record keeping
        layoutToUse = QtGui.QGridLayout()
        currentRow = 0
        self.widgetInfo = { }
        
        # set up the file info for the A file
        currentRow = self._add_file_related_controls("A", layoutToUse, currentRow)
        # set up the file info for the B file
        currentRow = self._add_file_related_controls("B", layoutToUse, currentRow)
        
        # set up the epsilon input box
        layoutToUse.addWidget(QtGui.QLabel("epsilon:"), currentRow, 0)
        self.epsilonWidget = QtGui.QLineEdit()
        tempValidator = QtGui.QDoubleValidator(self.epsilonWidget)
        tempValidator.setBottom(0.0) # only accept positive epsilons
        tempValidator.setNotation(QtGui.QDoubleValidator.StandardNotation)
        self.epsilonWidget.setValidator(tempValidator)
        self.epsilonWidget.editingFinished.connect(self.reportEpsilonChanged)
        layoutToUse.addWidget(self.epsilonWidget, currentRow, 1, 1, 2)
        
        currentRow = currentRow + 1
        
        # set up the drop down to allow image type selection
        layoutToUse.addWidget(QtGui.QLabel("Image Type:"), currentRow, 0)
        self.imageSelectionDropDown = QtGui.QComboBox()
        self.imageSelectionDropDown.activated.connect(self.reportImageTypeSelected)
        layoutToUse.addWidget(self.imageSelectionDropDown, currentRow, 1, 1, 2)
        
        currentRow = currentRow + 1
        
        # TODO should I add a drop down to select how to visualize the data? (ie 2D, line, on a map, etc?)
        
        # set up a button that shows stats
        self.statsButton = QtGui.QPushButton("Display Statistics")
        self.statsButton.clicked.connect(self.reportDisplayStatsClicked)
        layoutToUse.addWidget(self.statsButton, currentRow, 1, 1, 2)
        
        # set up the button at the bottom that creates plots
        self.displayButton = QtGui.QPushButton("Display Plot")
        self.displayButton.clicked.connect(self.reportDisplayPlotClicked)
        layoutToUse.addWidget(self.displayButton, currentRow, 3, 1, 2)
        
        # set up the overall window geometry
        self.setLayout(layoutToUse)
        self.setGeometry(600, 600, 625, 700)
    
    ################# start methods related to user input #################
    
    def clickedALoad (self) :
        
        self.selectFileToLoad ("A")
    
    def clickedBLoad (self) :
        
        self.selectFileToLoad ("B")
    
    def selectFileToLoad (self, file_prefix) :
        """
        when the load button is pressed, let the user pick a file to load
        """
        
        # get the file from the user
        tempFilePath = QtGui.QFileDialog.getOpenFileName(self, 'Open ' + str(file_prefix) + ' file', self.lastFilePath)
        
        LOG.debug ("New " + str(file_prefix) + " file selected in GUI: " + str(tempFilePath))
        
        if (tempFilePath is not None) and (len(tempFilePath) > 0) :
            self.lastFilePath = os.path.dirname(str(tempFilePath))
        
        # let our listeners know that the user picked a file
        for listener in self.userUpdateListeners :
            listener.newFileSelected(file_prefix, tempFilePath)
    
    def selectedAVariable (self) :
        
        self.reportVariableSelected("A")
    
    def selectedBVariable (self) :
        
        self.reportVariableSelected("B")
    
    def reportVariableSelected (self, file_prefix) :
        """
        when a variable is selected for one of the files, report it to any user update listeners
        """
        
        selectionText = self.widgetInfo[file_prefix]['variable'].currentText()
        
        # let our listeners know the user selected a variable
        for listener in self.userUpdateListeners :
            listener.userSelectedVariable(file_prefix, selectionText)
    
    def toggledAOverride (self) :
        
        self.reportOverrideChange("A")
    
    def toggledBOverride (self) :
        
        self.reportOverrideChange("B")
    
    def reportOverrideChange (self, file_prefix) :
        """
        when the user checks or unchecks one of the override checkboxes, report it to user update listeners
        """
        
        # this must be recorded before we tamper with the focus, because that will
        # trigger other events that may erase this information temporarily
        shouldDoOverride = self.widgetInfo[file_prefix]['override'].isChecked()
        
        # first we need to clean up focus in case it's in one of the line-edit boxes
        self.setFocus()
        
        # let our listeners know the user changed an overload setting
        for listener in self.userUpdateListeners :
            listener.userChangedOverload(file_prefix, shouldDoOverride)
    
    def fillValueEditedA (self) :
        
        self.fillValueChanged("A")
    
    def fillValueEditedB (self) :
        
        self.fillValueChanged("B")
    
    def fillValueChanged (self, file_prefix) :
        """
        when the user edits a fill value, report it to user update listeners
        """
        
        newFillValue = self.widgetInfo[file_prefix]['fillValue'].text()
        # it's still possible for this to not be a number, so fix that
        newFillValue = self._extra_number_validation(newFillValue)
        self.widgetInfo[file_prefix]['fillValue'].setText(str(newFillValue))
        
        # let our user update listeners know the fill value changed
        for listener in self.userUpdateListeners :
            listener.userChangedFillValue(file_prefix, newFillValue)
    
    def reportEpsilonChanged (self) :
        """
        when the epsilon changes, report it to user update listeners
        """
        
        newEpsilon = self.epsilonWidget.text()
        # it's still possible for epsilon to not be a number, so fix that
        newEpsilon = self._extra_number_validation(newEpsilon)
        self.epsilonWidget.setText(str(newEpsilon))
        
        # let our user update listeners know the epsilon changed
        for listener in self.userUpdateListeners :
            listener.userChangedEpsilon(newEpsilon)
    
    def reportImageTypeSelected (self) :
        """
        the user selected a new image type, so let our user update listeners know that
        """
        
        newImageType = self.imageSelectionDropDown.currentText()
        
        # report the new image type to our user update listeners
        for listener in self.userUpdateListeners :
            listener.userSelectedImageType(newImageType)
    
    def reportDisplayStatsClicked (self) :
        """
        the user clicked the display stats button
        """
        
        # make sure the focus isn't in a line-edit box
        self.statsButton.setFocus()
        
        # now report to our listeners that the user wants stats
        for listener in self.userUpdateListeners :
            listener.userRequestsStats()
    
    def reportDisplayPlotClicked (self) :
        """
        the user clicked the display plot button
        """
        
        # first we need to clean up focus in case it's in one of the line-edit boxes
        self.displayButton.setFocus()
        
        # now report to our listeners that the user wants a plot
        for listener in self.userUpdateListeners :
            listener.userRequestsPlot()
    
    def _extra_number_validation (self, string_that_should_be_a_number) :
        """
        try to validate the string that should be a number
        """
        
        toReturn = None
        
        try :
            toReturn = int(string_that_should_be_a_number)
        except ValueError :
            try :
                toReturn = float(string_that_should_be_a_number)
            except ValueError :
                pass # in this case we can't convert it, so just toss it
        
        return toReturn
    
    #################     end methods related to user input   #################
    
    ################# start data model update related methods #################
    
    def displayStatsData (self, aVariableName, bVariableName, statsAnalysis) :
        """
        given the names of the two variables and the statistical analysis,
        display this to the user
        """
        
        tempID            = self.statsCounter
        self.statsCounter = self.statsCounter + 1
        
        # I don't like this solution, but it would allow me to show multiple sets of stats at a time
        self.statsWindows[tempID] = StatisticsDisplayWindow(tempID,
                                                            aVariableName, variable_name_b=bVariableName,
                                                            statsTextToDisplay=str(statsAnalysis.dictionary_form()), stored_in=self.statsWindows)
                                                            #TODO, this is a terrible way to display this info, but shows that it is there
    
    def fileDataUpdate (self, file_prefix, file_path, selected_variable, use_fill_override, new_fill_value, variable_dimensions,
                        variable_list=None, attribute_list=None) :
        """
        The file data for one of the two files has changed. Update the GUI to reflect the change in data.
        """
        
        # set the path
        self.widgetInfo[file_prefix]['path'].setText(file_path)
        
        # if we got a new variable list, set up the list of variables
        if variable_list is not None :
            self.widgetInfo[file_prefix]['variable'].clear()
            self.widgetInfo[file_prefix]['variable'].addItems(variable_list)
        # set the selected variable
        tempPosition = self.widgetInfo[file_prefix]['variable'].findText(selected_variable)
        self.widgetInfo[file_prefix]['variable'].setCurrentIndex(tempPosition)
        
        # set the override
        self.widgetInfo[file_prefix]['override'].setChecked(use_fill_override)
        
        # set the fill value that's going to be used
        self.widgetInfo[file_prefix]['fillValue'].setText(str(new_fill_value)) # TODO, this should not be a string
        
        # show the variable dimensions
        self.widgetInfo[file_prefix]['dims'].setText(str(variable_dimensions))
        
        # if needed, show the attribute list
        if attribute_list is not None :
            temp_table = self.widgetInfo[file_prefix]['attrs']
            temp_table.clearContents()
            temp_table.setRowCount(len(attribute_list.keys()))
            rowCounter = 0
            for attributeKey in sorted(attribute_list.keys()) :
                temp_table.setCellWidget(rowCounter, 0, QtGui.QLabel(str(attributeKey)))
                temp_table.setCellWidget(rowCounter, 1, QtGui.QLabel(str(attribute_list[attributeKey])))
                rowCounter = rowCounter + 1
        
        # if there is a file selected, enable some of the other controls
        if file_path != "" :
            self.widgetInfo[file_prefix]['variable'].setDisabled(False)
            self.widgetInfo[file_prefix]['override'].setDisabled(False)
            self.widgetInfo[file_prefix]['fillValue'].setDisabled(not use_fill_override)
    
    def updateEpsilon (self, epsilon) :
        """
        update the comparison epsilon displayed to the user
        """
        
        self.epsilonWidget.setText(str(epsilon))
        
    
    def updateImageTypes (self, imageType, list=None) :
        """
        update the image type that's selected,
        if the list is given, clear and reset the list of possible image types
        """
        
        # replace the list if needed
        if list is not None :
            self.imageSelectionDropDown.clear()
            self.imageSelectionDropDown.addItems(list)
        
        # change the currently selected image type
        tempPosition = self.imageSelectionDropDown.findText(imageType)
        self.imageSelectionDropDown.setCurrentIndex(tempPosition)
    
    ################# end data model update related methods #################
    
    def showWarning (self, warningMessage):
        """
        show the user a warning dialog to the user
        """
        
        tempMessageBox = QtGui.QMessageBox()
        tempMessageBox.setText(warningMessage)
        tempMessageBox.exec_()
    
    def registerUserUpdateListener (self, objectToRegister) :
        """
        add the given object to our list of listeners
        """
        
        if objectToRegister not in self.userUpdateListeners :
            self.userUpdateListeners.append(objectToRegister)

class StatisticsDisplayWindow (QtGui.QWidget) :
    """
    This class represents a window that displays a statistics comparison between two variables.
    This window is intended to be displayed in a non-modal manner.
    """
    
    def __init__ (self, id_number, variable_name_a, variable_name_b=None,
                  statsTextToDisplay="", stored_in=None, parent=None) :
        """
        set up a window to display stats
        """
        
        QtGui.QWidget.__init__(self, parent)
        
        self.id     = id_number
        self.stored = stored_in
        
        # build and set the window title
        tempTitle = "Statistics Comparing " + str(variable_name_a)
        if (variable_name_b is not None) and (variable_name_b != variable_name_a) :
            tempTitle = tempTitle + " / " + str(variable_name_b)
        tempTitle = tempTitle + " data"
        self.setWindowTitle(tempTitle)
        
        # create the layout and set up some of the overall record keeping
        layoutToUse = QtGui.QGridLayout()
        
        # set up the button at the bottom that creates plots
        self.statsText = QtGui.QTextEdit(statsTextToDisplay)
        layoutToUse.addWidget(self.statsText, 1, 1)
        
        # set up the overall window geometry
        self.setLayout(layoutToUse)
        self.setGeometry(400, 400, 500, 500)
        
        self.show()
    
    def closeEvent (self, event) :
        """
        we need to clean some stuff up when the window wants to close
        """
        
        if self.stored is not None :
            del self.stored[self.id]
        
        event.accept()

if __name__=='__main__':
    import doctest
    doctest.testmod()
