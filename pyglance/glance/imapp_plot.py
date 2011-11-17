#!/usr/bin/env python
# encoding: utf-8
"""
Plot IMAPP IDEA data.

Created by evas Oct 2011.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
"""

import sys, os, logging

# these first two lines must stay before the pylab import
import matplotlib
matplotlib.use('Agg') # use the Anti-Grain Geometry rendering engine

from pylab import *

import matplotlib.cm     as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap

import numpy as np

import glance.data as dataobj

LOG = logging.getLogger(__name__)

defaultValues   = {
                    'longitudeVar': 'xtraj',
                    'latitudeVar':  'ytraj',
                    'initAODVar':   'aod_traj',
                    'trajPressVar': 'ptraj',
                    'timeVar':      'time',
                    'nePiece':      'NE_',
                    'swPiece':      'SW_',
                    'latPiece':     'LAT',
                    'lonPiece':     'LON',
                    'figureName':   'frame.png',
                    'figureDPI':    200
                  }

# a custom colormap or the Trajectory Pressures
color_data = {
    'red'   : ( (0.0,  1.0,  1.0),
                (0.75, 1.0,  1.0),
                (0.9,  1.0,  1.0),
                (1.0,  0.0,  0.0) ),
    'green' : ( (0.0,  1.0,  1.0),
                (0.75, 1.0,  1.0),
                (0.9,  0.08, 0.08),
                (1.0,  0.0,  0.0) ),
    'blue'  : ( (0.0,  1.0,  1.0),
                (0.75, 1.0,  1.0),
                (0.9,  0.58, 0.58),
                (1.0,  0.0,  0.0) )
              }
dark_trajectory_pressure_color_map = matplotlib.colors.LinearSegmentedColormap('darkTrajPressCM', color_data, 256)
color_data = {
    'red'   : ( (0.0,  0.0,  0.0),
                (0.75, 0.0,  0.0),
                (0.9,  1.0,  1.0),
                (1.0,  1.0,  1.0) ),
    'green' : ( (0.0,  0.0,  0.0),
                (0.75, 0.0,  0.0),
                (0.9,  0.08, 0.08),
                (1.0,  1.0,  1.0) ),
    'blue'  : ( (0.0,  0.0,  0.0),
                (0.75, 0.0,  0.0),
                (0.9,  0.58, 0.58),
                (1.0,  1.0,  1.0) )
              }
light_trajectory_pressure_color_map = matplotlib.colors.LinearSegmentedColormap('lightTrajPressCM', color_data, 256)

def _create_imapp_figure (initAODdata,       initLongitudeData,       initLatitudeData,
                          pressureData=None, pressLongitudeData=None, pressLatitudeData=None,
                          baseMapInstance=None, figureTitle="MODIS AOD & AOD Trajectories",
                          useDarkBackground=True, latRange=None, lonRange=None) :
    """
    TODO, I'm still deciding how this funciton works so it will probably change drastically
    """
    
    # build the plot
    figure = plt.figure()
    tempColor = 'k' if useDarkBackground else 'w' # either use a black or white background
    axes = figure.add_subplot(111, axisbg=tempColor)
    
    # build extra info to go to the map plotting function
    kwargs = { } 
    
    # draw the basic physical and geopolitical features
    tempColor='w' if useDarkBackground else 'k' # either draw black or white lines
    if baseMapInstance is not None :
        baseMapInstance.drawcoastlines( color=tempColor, linewidth=0.5)
        baseMapInstance.drawcountries(  color=tempColor, linewidth=0.5)
        baseMapInstance.drawstates(     color=tempColor, linewidth=0.5)
        baseMapInstance.drawmapboundary(color=tempColor, linewidth=0.5)
    # draw the parallels and meridians
    if latRange is not None :
        parallels = arange(-80., 90.,  latRange / 4.0)
        baseMapInstance.drawparallels(parallels,labels=[1,0,0,1], color=tempColor, linewidth=0.5)
    if lonRange is not None :
        meridians = arange(0.,   360., lonRange / 4.0)
        baseMapInstance.drawmeridians(meridians,labels=[1,0,0,1], color=tempColor, linewidth=0.5)
    
    # translate the longitude and latitude data sets into map coordinates
    pressX, pressY = None, None
    if pressureData is not None :
        pressX, pressY =  baseMapInstance(pressLongitudeData, pressLatitudeData)
    initX,  initY  = baseMapInstance(initLongitudeData,  initLatitudeData)
    
    
    color_map_to_use = dark_trajectory_pressure_color_map if useDarkBackground else light_trajectory_pressure_color_map
    # if we have pressure data plot that
    if pressureData is not None :
        # I'm taking advantage of the fact that I can remove the black edge line with lw=0 (line width = 0) and scale the marker down with s=0.5
        baseMapInstance.scatter(pressX, pressY, s=0.5, c=pressureData, marker='o', cmap=color_map_to_use, vmin=0, vmax=1000, lw=0)
        
        # create the pressure colorbar
        cbar1 = colorbar(format='%.3g', orientation='horizontal', shrink=0.25)
        cbar1.ax.set_position([0.1, -0.16, 0.25, 0.25])
        for tempText in cbar1.ax.get_xticklabels():
            #print(tempText)
            tempText.set_fontsize(5)
        cbar1.set_label("Trajectory Pressure (mb)")
    
    # plot the origin points after the pressure so they'll appear on top, I'm assuming we will always have this data
    baseMapInstance.scatter(initX,  initY, s=10,  c=initAODdata,  marker='o', cmap=cm.jet, vmin=0.0, vmax=1.0, lw=0.5)
    
    # make a color bar
    cbar2 = colorbar(format='%.3g', orientation='horizontal', shrink=0.25)
    cbar2.ax.set_position([0.4, -0.16, 0.25, 0.25])
    for tempText in cbar2.ax.get_xticklabels():
        #print(tempText)
        tempText.set_fontsize(5)
    cbar2.set_label("MODIS AOD") # TODO, how do I get a second Trajectory Pressure colorbar in the right place?
    
    # now that we've moved everything around, make sure our main image is in the right place
    axes.set_position([0.1, 0.15, 0.8, 0.8]) # why was this method so hard to find?
    
    # set up the figure title
    # TODO compose the figure title with time/date info?
    axes.set_title(figureTitle)
    
    return figure

def main():
    import optparse
    usage = """
%prog [options] 
run "%prog help" to list commands
examples:

python -m glance.imapp_plot plot A.nc

"""
    # the following represent options available to the user on the command line:
    
    parser = optparse.OptionParser(usage)
    
    # logging output options
    parser.add_option('-q', '--quiet', dest="quiet",
                    action="store_true", default=False, help="only error output")
    parser.add_option('-v', '--verbose', dest="verbose",
                    action="store_true", default=False, help="enable more informational output")   
    parser.add_option('-w', '--debug', dest="debug",
                    action="store_true", default=False, help="enable debug output")   

    # file generation related options TODO, make this work
    parser.add_option('-p', '--outpath', dest="outpath", type='string', default='./',
                    help="set path to output directory")
    
    # file variable settings
    parser.add_option('-o', '--longitude', dest="longitudeVar", type='string',
                    help="set name of longitude variable")
    parser.add_option('-a', '--latitude', dest="latitudeVar", type='string',
                    help="set name of latitude variable")
    
    # time related settings
    parser.add_option('-s', '--startTime', dest="startTime", type='int',
                    default=0, help="set first time to process")
    parser.add_option('-e', '--endTime', dest="endTime", type='int',
                    help="set last time to process")
    parser.add_option('-f', '--futureWindow', dest="futureWindow", type='int',
                    default=6, help="set number of hours of future pressures to show")
    
    parser.add_option('-t', '--test', dest="self_test",
                action="store_true", default=False, help="run internal unit tests")
    
    # TODO will add this once we're out of alpha
    #parser.add_option('-n', '--version', dest='version',
    #                  action="store_true", default=False, help="view the glance version")
    
    
    # parse the uers options from the command line
    options, args = parser.parse_args()
    if options.self_test:
        import doctest
        doctest.testmod()
        sys.exit(2)
    
    # set up the logging level based on the options the user selected on the command line
    lvl = logging.WARNING
    if options.debug: lvl = logging.DEBUG
    elif options.verbose: lvl = logging.INFO
    elif options.quiet: lvl = logging.ERROR
    logging.basicConfig(level = lvl)
    
    # TODO display the version
    #if options.version :
    #    print (_get_version_string() + '\n')

    commands = {}
    prior = None
    prior = dict(locals())
    
    """
    The following functions represent available menu selections.
    """
    
    def plot(*args):
        """plot trajectory frames
        Given a file with trajectory possitions and pressures over time, plot out
        images of these trajectories on the Earth and save it to disk.
        """
        
        LOG.debug("startTime:    " + str(options.startTime))
        LOG.debug("endTime:      " + str(options.endTime))
        LOG.debug("futureWindow: " + str(options.futureWindow))
        
        # setup the output directory now
        if not (os.path.isdir(options.outpath)) :
            LOG.info("Specified output directory (" + options.outpath + ") does not exist.")
            LOG.info("Creating output directory.")
            os.makedirs(options.outpath)
        
        # open the file
        LOG.info("Opening trajectory data file.")
        trajectoryFilePath = args[0]
        trajectoryFileObject = dataobj.FileInfo(trajectoryFilePath)
        if trajectoryFileObject is None:
            LOG.warn("Trajectory file (" + trajectoryFilePath + ") could not be opened.")
            LOG.warn("Aborting attempt to plot trajectories.")
            sys.exit(1)
        
        # load the required variables
        # TODO, allow the user control over the names?
        LOG.info("Loading variable data from trajectory data file.")
        initialAODdata         = trajectoryFileObject.file_object[defaultValues['initAODVar']]
        trajectoryPressureData = trajectoryFileObject.file_object[defaultValues['trajPressVar']]
        latitudeData           = trajectoryFileObject.file_object[defaultValues['latitudeVar']]
        longitudeData          = trajectoryFileObject.file_object[defaultValues['longitudeVar']]
        trajectoryTimeData     = trajectoryFileObject.file_object[defaultValues['timeVar']]
        
        # get information on where we should display the data
        northeastLon = trajectoryFileObject.file_object.get_global_attribute( defaultValues['nePiece'] + defaultValues['lonPiece'] )
        northeastLat = trajectoryFileObject.file_object.get_global_attribute( defaultValues['nePiece'] + defaultValues['latPiece'] )
        southwestLon = trajectoryFileObject.file_object.get_global_attribute( defaultValues['swPiece'] + defaultValues['lonPiece'] )
        southwestLat = trajectoryFileObject.file_object.get_global_attribute( defaultValues['swPiece'] + defaultValues['latPiece'] )
        latRange     = abs(northeastLat - southwestLat)
        lonRange     = abs(southwestLon - northeastLon)
        
        # build a basemap
        LOG.info("Building basemap object.")
        projectionName = 'merc' # use the Mercator Projection
        basemapObject  = Basemap (projection=projectionName,llcrnrlat=southwestLat,urcrnrlat=northeastLat,
                                  llcrnrlon=southwestLon, urcrnrlon=northeastLon, lat_ts=20, resolution='l') # TODO do I need to use lat_ts=20, ?
        
        # sort out the times we're using
        futureWindow = options.futureWindow
        startTime    = options.startTime
        endTime      = options.endTime if options.endTime is not None else trajectoryTimeData[-1]
        # as a sanity check, it's not really productive to make frames for times after our data ends
        if endTime > trajectoryTimeData[-1] :
            endTime = trajectoryTimeData[-1]
        
        # loop over time to create each frame
        initAODFlat =   initialAODdata.ravel()
        initLonData = longitudeData[0].ravel()
        initLatData =  latitudeData[0].ravel()
        for currentTime in range (startTime, endTime + 1) :
            
            # select only the window of trajectory data we need
            alreadyHappenedTimeIndex = 0
            tooFarInFutureTimeIndex  = trajectoryTimeData.size
            for tempIndex in range (0, trajectoryTimeData.size) :
                
                # first find the time that corresponds to the "most recently happened" index
                if trajectoryTimeData[tempIndex] <= currentTime :
                # TODO if time doesn't always increase, this needs another check
                    alreadyHappenedTimeIndex = tempIndex
                
                # then figure out how much further we can go before we pass the futureWindow's edge
                if (trajectoryTimeData[tempIndex] > (currentTime + futureWindow)) :
                # TODO, if time data doesn't always increase I also need & (trajectoryTimeData[tempIndex] < trajectoryTimeData[tooFarInFutureTimeIndex])
                    tooFarInFutureTimeIndex = tempIndex
                    break # TODO, this assumes time data always increases; also breaks suck so eventually take this out
            
            # only get data to plot the trajectory points if there is some available
            thisFramePressures = None
            thisFramePressLon  = None
            thisFramePressLat  = None
            if alreadyHappenedTimeIndex + 1 < tooFarInFutureTimeIndex :
                
                LOG.debug("Already happened index: " + str(alreadyHappenedTimeIndex))
                LOG.debug("Too far future index:   " + str(tooFarInFutureTimeIndex))
                
                # now pull out the pressure data and related lon/lat info for plotting
                thisFramePressures = trajectoryPressureData[alreadyHappenedTimeIndex + 1 : tooFarInFutureTimeIndex].ravel()
                thisFramePressLon  =          longitudeData[alreadyHappenedTimeIndex + 1 : tooFarInFutureTimeIndex].ravel()
                thisFramePressLat  =           latitudeData[alreadyHappenedTimeIndex + 1 : tooFarInFutureTimeIndex].ravel()
                
            # make the plot
            LOG.info("Creating trajectory plot for time " + str(float(currentTime)) + ".")
            tempFigure = _create_imapp_figure (initAODFlat,                     initLonData,                          initLatData,
                                               pressureData=thisFramePressures, pressLongitudeData=thisFramePressLon, pressLatitudeData=thisFramePressLat,
                                               baseMapInstance=basemapObject, latRange=latRange, lonRange=lonRange)
            
            # save the plot to disk
            LOG.info("Saving plot to disk.")
            tempFigure.savefig(os.path.join(options.outpath, str(currentTime) + defaultValues['figureName']), dpi=defaultValues['figureDPI'])
            
            # get rid of the figure 
            plt.close(tempFigure)
            del(tempFigure)
    
    def help(command=None):
        """print help for a specific command or list of commands
        e.g. help stats
        """
        if command is None: 
            # print first line of docstring
            for cmd in commands:
                ds = commands[cmd].__doc__.split('\n')[0]
                print "%-16s %s" % (cmd,ds)
        else:
            print commands[command].__doc__
            
    # def test():
    #     "run tests"
    #     test1()
    #
    
    # all the local public functions are considered part of this program, collect them up
    commands.update(dict(x for x in locals().items() if x[0] not in prior))    
    
    # if what the user asked for is not one of our existing functions, print the help
    if (not args) or (args[0] not in commands): 
        parser.print_help()
        help()
        return 9
    else:
        # call the function the user named, given the arguments from the command line  
        locals()[args[0]](*args[1:])

    return 0


if __name__=='__main__':
    sys.exit(main())