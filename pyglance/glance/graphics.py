#!/usr/bin/env python
# encoding: utf-8
"""
$Id: graphics.py,v 1.3 2008/05/08 16:46:02 rayg Exp $
Library of routines for plotting HSR data on a map

This file is based on the Keoni.map.graphics but has a significantly changed
interface to support more modular interactions with the basemap object.
~ Eva Schiffer, August 20th, 2009
"""

from mpl_toolkits.basemap import Basemap, shiftgrid
from numpy import arange, array, reshape, concatenate, nan

# the value that will denote "bad" longitudes and latitudes
badLonLat = 1.0E30

def create_basemap (lon, lat=None, axis=None, projection='lcc', resolution='i') :
    """
    Create an instance of basemap using either the specified axis info or the
    specified lon and lat info to pick the viewing area.
    the format of the axis is [lon min, lon max, lat min, lat max]
    where the min and max are for the entire area that you wish to show.
    
    Note: There are known viewing area problems with conic projections that
    may cause "rectangular" data to be clipped.
    """
    
    if lat is None:
        # then assume (lon,lat) pairs in a generator
        perimeter = list(lon)
        lon = array([x[0] for x in perimeter])
        lat = array([x[1] for x in perimeter])
    
    # make sure the axis is the correct shape
    # and fill in any missing values
    if axis is None: 
        axis = [None,None,None,None]
    if  axis[0] is None:
        axis[0] = lon.min()
    if  axis[1] is None:
        axis[1] = lon.max()
    if  axis[2] is None:
        axis[2] = lat.min()
    if  axis[3] is None:
        axis[3] = lat.max()
    
    # pull out the longitude/latitude info
    lon_left   = axis[0] 
    lat_bottom = axis[2] 
    lon_right  = axis[1] 
    lat_top    = axis[3] 
    lon_mid    = (lon_left + lon_right ) / 2.
    lat_mid    = (lat_top  + lat_bottom) / 2.
    
    # make our basemap
    m = None
    if projection is 'ortho' :
        # orthographic projections require this call
        m = Basemap(resolution=resolution, area_thresh=10000., projection=projection,
                    lat_0=lat_mid, lon_0=lon_mid)
    else :
        # most of the other projections use this call
        m = Basemap(llcrnrlon=lon_left,llcrnrlat=lat_bottom,urcrnrlon=lon_right,urcrnrlat=lat_top,
                    resolution=resolution, area_thresh=10000., projection=projection,
                    lat_1=lat_mid,lon_0=lon_mid)
    
    return m, axis

def draw_basic_features(baseMapInstance, axis) :
    """
    Draw the basic outlines of the earth's features.
    """
    # draw the basic physical and geopolitical features
    baseMapInstance.drawcoastlines()
    baseMapInstance.drawcountries()
    baseMapInstance.drawstates()
    baseMapInstance.drawmapboundary()
    
    # pull out the longitude/latitude info
    lon_left   = axis[0] 
    lat_bottom = axis[2] 
    lon_right  = axis[1] 
    lat_top    = axis[3] 
    
    # draw the parallels and meridians
    parallels = arange(-80.,90.,abs(lat_top - lat_bottom) / 4.0)
    baseMapInstance.drawparallels(parallels,labels=[1,0,0,1])
    meridians = arange(0., 360.,abs(lon_left - lon_right) / 4.0)
    baseMapInstance.drawmeridians(meridians,labels=[1,0,0,1])    
    
    return

def show_lon_lat_data(lon, lat, baseMapInstance, data=None, levelsToUse=None, **kwargs) :
    """
    Show data corresponding to the longitude and latitude set provided on the earth using the provided basemap.
    levelsToUse is a list of numbers representing data ranges that will be used
    """
    x, y = baseMapInstance(lon, lat) # translate into the coordinate system of the basemap
    
    return show_x_y_data(x, y, baseMapInstance, data, levelsToUse, **kwargs)

def show_x_y_data(x, y, baseMapInstance, data=None, levelsToUse=None, **kwargs) :
    """
    Show data corresponding to a given x, y using the provided basemap.
    levelsToUse is a list of numbers representing data ranges that will be used
    """
    artistsAdded = [ ]
    
    # only try to plot the data if there is some
    if data is not None:
        
        newX    = make_2D_version_if_needed(x,    badLonLat)
        newY    = make_2D_version_if_needed(y,    badLonLat)
        newData = make_2D_version_if_needed(data, nan)
        
        if levelsToUse is not None :
            p = baseMapInstance.contourf(newX, newY, newData, levelsToUse, **kwargs)
        else :
            p = baseMapInstance.contourf(newX, newY, newData, **kwargs)
    
    # return the original x and y so the caller can match any external data in shape
    return baseMapInstance, x, y

def show_quiver_plot (lon, lat, baseMapInstance, (uData, vData)=(None,None), colordata=None, **kwargs) :
    """
    Show a quiver plot of the given vector data at the given longitude and latitude
    """
    
    x, y = baseMapInstance(lon, lat) # translate into the coordinate system of the basemap
    
    # show the quiver plot if there is data
    
    if (uData is not None) and (vData is not None) :
        if colordata is None:
            p = baseMapInstance.quiver(x, y, uData, vData, **kwargs)
        else :
            p = baseMapInstance.quiver(x, y, uData, vData, colordata, **kwargs)
    
    # return the original x and y so the caller can match any external data in shape
    return baseMapInstance, x, y

def make_2D_version_if_needed(originalArray, fillValue) :
    """
    if the array isn't already 2D, add some fillValued points and reshape it
    to make it so
    """
    arrayToChange = originalArray
    
    # if the array has too few dimentions, we'll need to reshape it
    if len(arrayToChange.shape) < 2 :
            
            # if there are an odd number of points, add a null point
            if (arrayToChange.size % 2) > 0 :
                arrayToChange = concatenate((arrayToChange, [fillValue]))
            
            # now fold the array in half
            newDim = arrayToChange.size / 2
            arrayToChange = reshape(arrayToChange, (newDim, 2))
    
    return arrayToChange

def mapshow(lon, lat=None, data=None, axis=None, projection='lcc',resolution='i', levelsToUse=None, **kwargs):
    """Like matshow(), but requiring longitude and latitude
    optional axis is in the form of [left-longitude, right-longitude, bottom-latitude, top-latitude]
    axis values left as None will be auto-maximized
    returns basemap object, x-transform, y-transform
    if lat is None, then assume that lon is a sequence of (lon,lat) points
    """
    # make the basemap object
    baseMapInstance, axis = create_basemap(lon, lat, axis, projection, resolution)
    
    # draw the minimal map features
    draw_basic_features(baseMapInstance, axis)
    
    # plot the actual data on the map
    baseMapInstance, x, y = show_lon_lat_data(lon, lat, baseMapInstance, data, levelsToUse, **kwargs)
    
    return baseMapInstance, x, y

