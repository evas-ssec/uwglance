#!/usr/bin/env python
# encoding: utf-8
"""
The view portion of the glance GUI code. 

Created by evas Oct 2011.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
"""

import sys, os.path, logging

from PyQt4 import QtGui, QtCore

from functools import partial

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
        
        # a place to hang onto our file specific widgets
        self.widgetInfo = { }
        
        # setup the rest of our window
        # TODO, this section is temporary, set up the tabs in a sub function?
        self.tabWidget = QtGui.QTabWidget()
        
        # add the main tab that allows for basic control
        tempWidget = QtGui.QWidget()
        tempWidget.setLayout(self._build_data_tab())
        self.tabWidget.addTab(tempWidget, "basic")
        
        # TODO uncomment this to work on the settings tab
        # add a tab that allows more detailed, optional settings
        #tempWidget = QtGui.QWidget()
        #tempWidget.setLayout(self._build_settings_tab())
        #self.tabWidget.addTab(tempWidget, "settings")
        
        tempLayout = QtGui.QGridLayout()
        tempLayout.addWidget(self.tabWidget)
        self.setLayout(tempLayout)
        self.setGeometry(600, 600, 625, 700)
        
        # this will represent those who want to be notified
        # when the user changes things
        self.userUpdateListeners = [ ]
        
        # we will use this to remember were the user wanted to load files from last
        # TODO, can we remember this between program runs? something like preferences?
        self.lastFilePath = './'
        
        # hang on to stats windows so they don't vanish
        self.statsWindows = { }
        self.statsCounter = 1
    
    def _build_data_tab (self) :
        """
        built the qt controls for the basic data tab and lay them out in a grid layout
        """
        
        # create the layout and set up some of the overall record keeping
        layoutToUse = QtGui.QGridLayout()
        currentRow = 0
        
        # set up the file info for the A file
        currentRow = self._add_file_related_controls("A", layoutToUse, currentRow)
        # set up the file info for the B file
        currentRow = self._add_file_related_controls("B", layoutToUse, currentRow)
        
        # set up the epsilon input box
        layoutToUse.addWidget(QtGui.QLabel("epsilon:"), currentRow, 0)
        self.epsilonWidget = QtGui.QLineEdit()
        self.epsilonWidget.setToolTip("Maximum acceptible difference between the variable data in the two files.")
        tempValidator = QtGui.QDoubleValidator(self.epsilonWidget)
        tempValidator.setBottom(0.0) # only accept positive epsilons
        tempValidator.setNotation(QtGui.QDoubleValidator.StandardNotation)
        self.epsilonWidget.setValidator(tempValidator)
        self.epsilonWidget.editingFinished.connect(self.reportEpsilonChanged)
        layoutToUse.addWidget(self.epsilonWidget, currentRow, 1, 1, 2)
        
        currentRow = currentRow + 1
        
        # set up the epsilon percent input box
        layoutToUse.addWidget(QtGui.QLabel("epsilon percent:"), currentRow, 0)
        self.epsilonPerWidget = QtGui.QLineEdit()
        self.epsilonPerWidget.setToolTip("Maximum acceptible difference between the variable data in terms of % of each data point in the A file.")
        tempValidator = QtGui.QDoubleValidator(self.epsilonPerWidget)
        tempValidator.setBottom(0.0) # only accept positive epsilon percents
        tempValidator.setNotation(QtGui.QDoubleValidator.StandardNotation)
        self.epsilonPerWidget.setValidator(tempValidator)
        self.epsilonPerWidget.editingFinished.connect(self.reportEpsilonPercentChanged)
        layoutToUse.addWidget(self.epsilonPerWidget, currentRow, 1, 1, 2)
        
        currentRow = currentRow + 1
        
        # set up the drop down to allow image type selection
        layoutToUse.addWidget(QtGui.QLabel("Image Type:"), currentRow, 0)
        self.imageSelectionDropDown = QtGui.QComboBox()
        self.imageSelectionDropDown.activated.connect(self.reportImageTypeSelected)
        layoutToUse.addWidget(self.imageSelectionDropDown, currentRow, 1, 1, 2)
        
        currentRow = currentRow + 1
        
        # set up a button that shows stats
        self.statsButton = QtGui.QPushButton("Display Statistics")
        self.statsButton.clicked.connect(self.reportDisplayStatsClicked)
        layoutToUse.addWidget(self.statsButton, currentRow, 1, 1, 2)
        
        # set up the button at the bottom that creates plots
        self.displayButton = QtGui.QPushButton("Display Plot")
        self.displayButton.clicked.connect(self.reportDisplayPlotClicked)
        layoutToUse.addWidget(self.displayButton, currentRow, 3, 1, 2)
        
        return layoutToUse
    
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
        filePath.setToolTip("The currently loaded file.")
        filePath.setReadOnly(True) # this is mostly for displaying the file path, the load button selects it
        self.widgetInfo[file_prefix]['path'] = filePath
        grid_layout.addWidget(filePath, currentRow, 1, 1, 3)
        loadButton = QtGui.QPushButton("Load")
        # set some tooltip text
        loadButton.setToolTip("Load a file: glance can handle NetCDF, HDF4, HDF5, and AERI files")
        # connect the button to an action
        loadButton.clicked.connect(partial(self.selectFileToLoad, file_prefix=file_prefix))
        self.widgetInfo[file_prefix]['load'] = loadButton
        grid_layout.addWidget(loadButton, currentRow, 4)
        
        currentRow = currentRow + 1
        
        # set up the drop down for the variable select
        grid_layout.addWidget(QtGui.QLabel("variable name:"), currentRow, 1)
        variableSelection = QtGui.QComboBox()
        variableSelection.setDisabled(True)
        variableSelection.activated.connect(partial(self.reportVariableSelected, file_prefix=file_prefix))
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
        overrideFillButton.setToolTip("Check to override the default fill value.")
        overrideFillButton.setDisabled(True)
        overrideFillButton.stateChanged.connect(partial(self.reportOverrideChange, file_prefix=file_prefix))
        self.widgetInfo[file_prefix]['override'] = overrideFillButton
        grid_layout.addWidget(overrideFillButton, currentRow, 1)
        
        # now set up the input of the fill value that will be used
        grid_layout.addWidget(QtGui.QLabel("fill value:"), currentRow+1, 1)
        fillValue = QtGui.QLineEdit()
        fillValue.setToolTip("The fill value that will be used.")
        tempValidator = QtGui.QDoubleValidator(fillValue)
        tempValidator.setNotation(QtGui.QDoubleValidator.StandardNotation)
        fillValue.setValidator(tempValidator)
        fillValue.setDisabled(True)
        fillValue.editingFinished.connect(partial(self.fillValueChanged, file_prefix=file_prefix))
        self.widgetInfo[file_prefix]['fillValue'] = fillValue
        grid_layout.addWidget(fillValue, currentRow+1, 2, 1, 3)
        
        currentRow = currentRow + 2
        
        return currentRow
    
    def _build_settings_tab (self) :
        """
        built the basic qt controls for the settings tab and lay them out in a grid layout
        """
        
        # create the layout and set up some of the overall record keeping
        layoutToUse = QtGui.QGridLayout()
        currentRow = 0
        
        # add the filter entry areas
        currentRow = self._add_filter_controls("A", layoutToUse, currentRow)
        currentRow = self._add_filter_controls("B", layoutToUse, currentRow)
        
        # add the drop down for selecting the data display form
        layoutToUse.addWidget(QtGui.QLabel("Data Form:"), currentRow, 0)
        self.dataDisplayFormDropDown = QtGui.QComboBox()
        self.dataDisplayFormDropDown.activated.connect(self.reportDataFormSelected)
        layoutToUse.addWidget(self.dataDisplayFormDropDown, currentRow, 1, 1, 2)
        
        currentRow = currentRow + 1
        
        # add a check box to constrain originals to the same range
        showOriginalsInSameRange = QtGui.QCheckBox("show original plots in same range")
        showOriginalsInSameRange.setToolTip("Check to constrain the colorbar range of the Original A and " +
                                            "Original B data plots to be the same.\nThe range will include all data in both A and B.")
        showOriginalsInSameRange.setDisabled(False)
        #showOriginalsInSameRange.stateChanged.connect(TODO)
        self.useSameRangeWidget = showOriginalsInSameRange
        layoutToUse.addWidget(showOriginalsInSameRange, currentRow, 1, 1, 2)
        
        currentRow = currentRow + 1
        
        # add lon/lat controls
        
        # add the lon/lat controls that are separated by file
        currentRow = self._add_lon_lat_controls("A", layoutToUse, currentRow)
        currentRow = self._add_lon_lat_controls("B", layoutToUse, currentRow)
        
        # add box to enter lon/lat epsilon
        layoutToUse.addWidget(QtGui.QLabel("lon/lat epsilon:"), currentRow, 0)
        llepsilonWidget = QtGui.QLineEdit()
        self.llepsilonWidget = llepsilonWidget
        llepsilonWidget.setToolTip("Maximum acceptible difference between the longitudes or latitudes in the two files.")
        tempValidator = QtGui.QDoubleValidator(llepsilonWidget)
        tempValidator.setBottom(0.0) # only accept positive epsilons
        tempValidator.setNotation(QtGui.QDoubleValidator.StandardNotation)
        llepsilonWidget.setValidator(tempValidator)
        llepsilonWidget.editingFinished.connect(self.reportLLEpsilonChanged)
        layoutToUse.addWidget(llepsilonWidget, currentRow, 1, 1, 2)
        
        # TODO, possible enhancements for the future
        #       add filters for lon/lat?
        #       allow lon/lat to be loaded from a different file?
        #       add buttons to let you plot lon/lat specific errors?
        
        return layoutToUse
    
    def _add_filter_controls (self, file_prefix, grid_layout, current_row) :
        """
        Add controls for the user to enter the filter function and any adjunct
        support code needed to set up that function using the given grid_layout
        
        return the next free current_row number when finished adding widgets
        """
        
        # add the entry to name the specific function
        grid_layout.addWidget(QtGui.QLabel(str(file_prefix) + " filter function:"), current_row, 0)
        functionEntry = QtGui.QLineEdit()
        functionEntry.setToolTip("Enter python function to use for filtering " + str(file_prefix) + " data here."
                                 + "\nThis function should work in the form function_name(data_to_filter)"
                                 + "\n and you should enter only the function_name here.")
        #functionEntry.editingFinished.connect(TODO)
        grid_layout.addWidget(functionEntry, current_row, 1, 1, 3)
        self.widgetInfo[file_prefix]["filter_function"] = functionEntry
        
        current_row = current_row + 1
        
        grid_layout.addWidget(QtGui.QLabel(str(file_prefix) + " filter setup: "), current_row, 0)
        setupCodeEntry = QtGui.QTextEdit()
        setupCodeEntry.setToolTip("Enter additional python code needed to set up the filter function for " + str(file_prefix) + ".")
        setupCodeEntry.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Preferred) # TODO, this is not working!
        #setupCodeEntry.editingFinished.connect(TODO)
        self.widgetInfo[file_prefix]["filter_setup"] = setupCodeEntry
        grid_layout.addWidget(setupCodeEntry, current_row, 1, 1, 3)
        
        current_row = current_row + 1
        
        return current_row
    
    def _add_lon_lat_controls (self, file_prefix, grid_layout, current_row) :
        """
        Add the longitude and latitude controls for the given file_prefix to
        the given grid_layout starting on current_row
        
        return the next free current_row number when finished adding widgets
        """
        
        # add the label so we know which file this is for
        grid_layout.addWidget(QtGui.QLabel(file_prefix + " Navigation:"), current_row, 0)
        
        current_row = current_row + 1
        
        # add drop down to select latitude
        grid_layout.addWidget(QtGui.QLabel("Latitude:"), current_row, 1)
        latNameDropDown = QtGui.QComboBox()
        latNameDropDown.activated.connect(partial(self.reportLatitudeSelected, file_prefix=file_prefix))
        self.widgetInfo[file_prefix]["latName"] = latNameDropDown
        grid_layout.addWidget(latNameDropDown, current_row, 2, 1, 2)
        
        current_row = current_row + 1
        
        # add drop down to select longitude
        grid_layout.addWidget(QtGui.QLabel("Longitude:"), current_row, 1)
        lonNameDropDown = QtGui.QComboBox()
        lonNameDropDown.activated.connect(partial(self.reportLongitudeSelected, file_prefix=file_prefix))
        self.widgetInfo[file_prefix]["lonName"] = lonNameDropDown
        grid_layout.addWidget(lonNameDropDown, current_row, 2, 1, 2)
        
        current_row = current_row + 1
        
        return current_row
    
    ################# start methods related to user input #################
    
    def selectFileToLoad (self, file_prefix=None) :
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
    
    def reportVariableSelected (self, file_prefix=None) :
        """
        when a variable is selected for one of the files, report it to any user update listeners
        """
        
        selectionText = self.widgetInfo[file_prefix]['variable'].currentText()
        
        # let our listeners know the user selected a variable
        for listener in self.userUpdateListeners :
            listener.userSelectedVariable(file_prefix, selectionText)
    
    def reportOverrideChange (self, file_prefix=None) :
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
    
    def fillValueChanged (self, file_prefix=None) :
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
    
    def reportLLEpsilonChanged (self) :
        """
        when the lon/lat epsilon changes, report it to user update listeners
        """
        
        newLLEpsilon = self.llepsilonWidget.text()
        # it's still possible for epsilon to not be a number, so fix that
        newLLEpsilon = self._extra_number_validation(newLLEpsilon)
        self.llepsilonWidget.setText(str(newLLEpsilon))
        
        # let our user update listeners know the llepsilon changed
        for listener in self.userUpdateListeners :
            listener.userChangedLLEpsilon(newLLEpsilon)
    
    
    def reportEpsilonPercentChanged (self) :
        """
        when the epsilon percent changes, report it to user update listeners
        """
        
        newEpsilonPer = self.epsilonPerWidget.text()
        # it's still possible for epsilon percent to not be a number, so fix that
        newEpsilonPer = self._extra_number_validation(newEpsilonPer)
        self.epsilonPerWidget.setText(str(newEpsilonPer))
        
        # let our user update listeners know the epsilon percent changed
        for listener in self.userUpdateListeners :
            listener.userChangedEpsilonPercent(newEpsilonPer)
    
    def reportLongitudeSelected (self, file_prefix=None) :
        """
        when a longitude variable is selected for one of the files, report it to any user update listeners
        """
        
        selectionText = self.widgetInfo[file_prefix]['lonName'].currentText()
        
        # let our listeners know the user selected a longitude
        for listener in self.userUpdateListeners :
            listener.userSelectedLongitude(file_prefix, selectionText)
    
    def reportLatitudeSelected (self, file_prefix=None) :
        """
        when a latitude variable is selected for one of the files, report it to any user update listeners
        """
        
        selectionText = self.widgetInfo[file_prefix]['latName'].currentText()
        
        # let our listeners know the user selected a latitude
        for listener in self.userUpdateListeners :
            listener.userSelectedLatitude(file_prefix, selectionText)
    
    def reportImageTypeSelected (self) :
        """
        the user selected a new image type, so let our user update listeners know that
        """
        
        newImageType = self.imageSelectionDropDown.currentText()
        
        # report the new image type to our user update listeners
        for listener in self.userUpdateListeners :
            listener.userSelectedImageType(newImageType)
    
    def reportDataFormSelected (self) :
        """
        the user selected a new data form, so let our user update listeners know that
        """
        
        newDataForm = self.dataDisplayFormDropDown.currentText()
        
        # report the new data form to our user update listeners
        for listener in self.userUpdateListeners :
            listener.userSelectedDataForm(newDataForm)
    
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
                                                            statsTextToDisplay=statsAnalysis, stored_in=self.statsWindows)
    
    def fileDataUpdate (self, file_prefix, file_path, selected_variable, use_fill_override, new_fill_value, variable_dimensions,
                        variable_list=None, attribute_list=None) :
        """
        The file data for one of the two files has changed. Update the GUI to reflect the change in data.
        """
        
        # set the path
        self.widgetInfo[file_prefix]['path'].setText(file_path)
        
        # if we got a new variable list, set up the appropriate drop down list
        if variable_list is not None :
            # set up the list of selectable variables for analysis
            self.widgetInfo[file_prefix]['variable'].clear()
            self.widgetInfo[file_prefix]['variable'].addItems(variable_list)
        
        # set the selected variable
        tempPosition = self.widgetInfo[file_prefix]['variable'].findText(selected_variable)
        self.widgetInfo[file_prefix]['variable'].setCurrentIndex(tempPosition)
        
        # set the override
        self.widgetInfo[file_prefix]['override'].setChecked(use_fill_override)
        
        # set the fill value that's going to be used
        self.widgetInfo[file_prefix]['fillValue'].setText(str(new_fill_value))
        
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
    
    def updateSelectedLatLon(self, filePrefix, newLatitude, newLongitude, lonlatList=None) :
        """
        Update the latitude and longitude names that are selected in the drop down,
        if a list is given, then replace the list of options that are being displayed for that file.
        """
        
        """ TODO, uncomment once that tab is set up
        
        # if we got a new list, set up the appropriate drop down lists
        if lonlatList is not None :
            
            # set up the longitude and latitude selectors
            self.widgetInfo[filePrefix]['latName'].clear()
            self.widgetInfo[filePrefix]['latName'].addItems(lonlatList)
            self.widgetInfo[filePrefix]['lonName'].clear()
            self.widgetInfo[filePrefix]['lonName'].addItems(lonlatList)
        
        # set the selected latitude
        tempPosition = self.widgetInfo[filePrefix]['latName'].findText(newLatitude)
        self.widgetInfo[filePrefix]['latName'].setCurrentIndex(tempPosition)
        
        # set the selected longitude
        tempPosition = self.widgetInfo[filePrefix]['lonName'].findText(newLongitude)
        self.widgetInfo[filePrefix]['lonName'].setCurrentIndex(tempPosition)
        
        """
    
    def updateEpsilon (self, epsilon) :
        """
        update the comparison epsilon displayed to the user
        """
        
        self.epsilonWidget.setText(str(epsilon))
    
    def updateEpsilonPercent (self, epsilonPercent) :
        """
        update the epsilon percent displayed to the user
        """
        
        self.epsilonPerWidget.setText(str(epsilonPercent))
    
    def updateLLEpsilon (self, newLonLatEpsilon) :
        """
        update the epsilon for longitude and latitude displayed to the user
        """
        
        #self.llepsilonWidget.setText(str(newLonLatEpsilon)) TODO, uncomment once that tab is set up
    
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
    
    def updateDataForms(self, dataForm, list=None) :
        """
        update the data form that's selected,
        if the list is given, clear and reset the list of possible data forms
        """
        
        """ TODO, uncomment once that tab is set up
        # replace the list if needed
        if list is not None :
            self.dataDisplayFormDropDown.clear()
            self.dataDisplayFormDropDown.addItems(list)
        
        # change the currently selected data form
        tempPosition = self.dataDisplayFormDropDown.findText(dataForm)
        self.dataDisplayFormDropDown.setCurrentIndex(tempPosition)
        """
    
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
        self.statsText = QtGui.QTextEdit()
        self.statsText.setHtml(statsTextToDisplay)
        self.statsText.setReadOnly(True)
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
