#!/usr/bin/env python
# encoding: utf-8
"""
This module handles physical collocation of two data sets. The code
in this module is based on previous versions of delta.py.

Created by evas Apr 2010.
Copyright (c) 2010 University of Wisconsin SSEC. All rights reserved.
"""

import logging
import numpy as np

LOG = logging.getLogger(__name__)

# TODO, move to using this class (not yet finished)
class CollocationMapping :
    """
    This class represents a mapping that collocates points in two
    sets of data within certain tolerances. These points may be
    filtered based on secondary criteria beyond spatial tolerances,
    but once the inital mapping is defined, an instance of this class
    can only be used to convert size-compatable data from the two
    original sets into the final, collocated order.
    
    Note: the mapping will include all pairs of points that match
    within the given epsilon and additional data criteria given,
    this means an individual a or b point may be repeated if it matches
    multiple points in the other data set satisfactorily
    
    WARNING: The algorithm used may fail to find all possible matches
    if a longitude/latitude epsilon of greater than or equal to one
    degree is given TODO
    """
    
    """
    raw_lon_lat_matches  - the number of points matched, based on the longitude and latitude alone
    raw_unmatchable_in_a - the number of unmacthed pointes in a, based on longitude and latitude alone
    raw_unmatchable_in_b - the number of unmacthed pointes in b, based on longitude and latitude alone
    
    a_point_mapping      - mapping of points in a to points in b; TODO give form
    b_point_mapping      - mapping of points in b to points in a; TODO give form
    
    matched_longitude    - the list of matched longitude values
    matched_latitude     - the list of matched latitude  values
    
    unmatchable_longitude_a - the list of longitude values in a that could not be matched
    unmatchable_latitude_a  - the list of latitude  values in a that could not be matched
    unmatchable_longitude_b - the list of longitude values in b that could not be matched
    unmatchable_latitude_b  - the list of latitude  values in b that could not be matched
    """
    
    def __init__(self,
                 (a_longitude, a_latitude),
                 (b_longitude, b_latitude),
                 lon_lat_epsilon,
                 valid_in_a_mask=None, valid_in_b_mask=None,
                 additional_data_criteria=[ ],
                 additional_filter_functions=[ ]) :
        """
        Build collocation mapping data based on two sets of
        longitude and latitude, as well as an epsilon that
        describes acceptable differences in degrees.
        
        Note: additional "criteria" variable data may be inculuded.
        additional_filter_functions will be called in the form:
        
            additional_filter_functions[i](aRow, aCol, bRow, bCol, additional_data_criteria[i])
        
        a return of True or False is expected, to indicate whether or not
        the match should be accepted
        """
        pass
    
    def _create_basic_mapping_from_lon_lat((alongitude, alatitude),
                                           (blongitude, blatitude),
                                           lonlatEpsilon,
                                           invalidAMask=None, invalidBMask=None) :
        """
        match points together based on their longitude and latitude values
        to match points must be within lonlatEpsilon degrees in both longitude and latitude
        
        if the longitude and latitude variables contain invalid data the invalidAMask and
        invalidBMask should be passed with the appropriate masking to remove the invalid values
        
        the return will be in the form of two dictionaries of points, one from a and one from b,
        indexed on the index number in the A or B data where they can be found. Each entry will
        consist of a list of:
            [longitudeValue, latitudeValue, indexNumber, [list of matching indexes in the other set]]
        
        Note: the return will include all pairs of points that match,
        this means an individual a or b point may be repeated if it matches
        multiple points within the lonlatEpsilon provided
        
        Warning: This algorithm will fail to find all matching points if the lonlatEpsilon is set to a
        value greater than or equal to 1.0 degrees. This is related to the bin size used for searching
        thoretically the bin size could be corrected to scale with the lonlatEpsilon in the future. TODO
        """

def create_colocation_mapping_within_epsilon((alongitude, alatitude),
                                             (blongitude, blatitude),
                                             lonlatEpsilon,
                                             invalidAMask=None, invalidBMask=None):
    """
    match points together based on their longitude and latitude values
    to match points must be within lonlatEpsilon degrees in both longitude and latitude
    
    if the longitude and latitude variables contain invalid data the invalidAMask and
    invalidBMask should be passed with the appropriate masking to remove the invalid values
    
    the return will be in the form of two dictionaries of points, one from a and one from b,
    indexed on the index number in the A or B data where they can be found. Each entry will
    consist of a list of:
        [longitudeValue, latitudeValue, indexNumber, [list of matching indexes in the other set]]
    
    Note: the return will include all pairs of points that match,
    this means an individual a or b point may be repeated if it matches
    multiple points within the lonlatEpsilon provided
    
    Warning: This algorithm will fail to find all matching points if the lonlatEpsilon is set to a
    value greater than or equal to 1.0 degrees. This is related to the bin size used for searching
    thoretically the bin size could be corrected to scale with the lonlatEpsilon in the future. TODO
    """
    assert(alongitude.shape == alatitude.shape)
    assert(blongitude.shape == blatitude.shape)
    assert(lonlatEpsilon >= 0.0)
    
    LOG.debug("Preparing to colocate longitude and latitude points (acceptable epsilon: " + str(lonlatEpsilon) + " degrees)")
    LOG.debug("size of A: " + str(alongitude.shape))
    LOG.debug("size of B: " + str(blongitude.shape))
    
    # make blank invalid masks if none were passed in
    if invalidAMask is None :
        invalidAMask = np.zeros(alongitude.shape, dtype=bool)
    if invalidBMask is None :
        invalidBMask = np.zeros(blongitude.shape, dtype=bool)
    
    # make flat versions of our longitude and latitude
    # so that our index correlations will be simple
    flatALatitude  =  alatitude.ravel()
    flatALongitude = alongitude.ravel()
    flatBLatitude  =  blatitude.ravel()
    flatBLongitude = blongitude.ravel()
    
    # find the ranges of the longitude and latitude
    minLatitude  = np.min(np.min(flatALatitude),  np.min(flatBLatitude))
    maxLatitude  = np.max(np.max(flatALatitude),  np.max(flatBLatitude))
    minLongitude = np.min(np.min(flatALongitude), np.min(flatBLongitude))
    maxLongitude = np.max(np.max(flatALongitude), np.max(flatBLongitude))
    
    # make the bins for the data in longitude/latitude space
    aBins = { }
    bBins = { }
    allAPts = { }
    allBPts = { }
    
    # loop to put all the aData in the bins
    for index in range(flatALatitude.size) :
        filingLat = int( flatALatitude[index])
        filingLon = int(flatALongitude[index])
        
        if (filingLat, filingLon) not in aBins :
            aBins[(filingLat, filingLon)] = [ ]
        
        # create the simple list holding that point in the form:
        # the lon/lat values (for ease of comparison), the index number in A, and the list of matches
        aPoint = [flatALatitude[index], flatALongitude[index], index, [ ]]
        
        # put the point in the list and bin
        allAPts[index] = aPoint
        aBins[(filingLat, filingLon)].append(aPoint)
    
    # loop to put all the bData in the bins
    for index in range(flatBLatitude.size) :
        filingLat = int( flatBLatitude[index])
        filingLon = int(flatBLongitude[index])
        
        if (filingLat, filingLon) not in bBins :
            bBins[(filingLat, filingLon)] = [ ]
        
        # create the simple list holding that point in the form:
        # the lon/lat values (for ease of comparison), the index number in A, and the list of matches
        bPoint = [flatBLatitude[index], flatBLongitude[index], index, [ ]]
        
        # put the point in the list and bin
        allBPts[index] = bPoint
        bBins[(filingLat, filingLon)].append(bPoint)
    
    # for debugging purposes
    totalMatches = 0
    
    # look in all the aData bins and select point pairs that match within epsilon
    for binLatitude, binLongitude in aBins.keys() :
        
        # figure out the lat/lon of the 9 bins that are "near" this one
        toSearch = [ ]
        for latValue in range(binLatitude - 1, binLatitude + 1) :
            for lonValue in range(binLongitude -1, binLongitude + 1) :
                toSearch.append((latValue, lonValue))
        
        # for each A pt in this bin
        for aLat, aLon, aIndex, aMatches in aBins[(binLatitude, binLongitude)] :
            
            # look through my nearby B bins and find any
            # "matching" points that fall within epsilon
            for latValue, lonValue in toSearch :
                
                # if there's anything in that B bin
                if ((latValue, lonValue) in bBins) and (bBins[(latValue, lonValue)] is not None) :
                    
                    # for each data point in the B bin, check if it matches our current A point
                    for bLat, bLon, bIndex, bMatches in bBins[(latValue, lonValue)] :
                        
                        if (abs(bLat - aLat) < lonlatEpsilon) and (abs(aLon - bLon) < lonlatEpsilon) :
                            totalMatches = totalMatches + 1
                            
                            # put the point on our matched lists
                            aMatches.append(bIndex)
                            bMatches.append(aIndex)
    
    LOG.debug('Found ' + str(totalMatches) + ' matched pairs.')
    
    return allAPts, allBPts, totalMatches

def create_colocated_lonlat_with_lon_lat_colocation(listOfColocatedALonLat, listOfColocatedBLonLat,
                                                    totalMatches,
                                                    aLongitude, aLatitude,
                                                    bLongitude, bLatitude) :
    """
    given a pre colocated list of A and B lon/lat info from create_colocation_mapping_within_epsilon,
    match up the longitude and latitude and return the colocated sets
    """
    
    # some general statistics
    multipleMatchesInA      = 0
    multipleMatchesInB      = 0
    totalValidMatchedPairs  = 0
    
    # our final data sets
    matchedLongitude   = np.zeros(totalMatches, dtype=aLongitude.dtype) 
    matchedLatitide    = np.zeros(totalMatches, dtype=aLatitude.dtype) 
    
    # we don't know how many unmatched points we may have
    unmatchedALongitude = [ ]
    unmatchedALatitude  = [ ]
    unmatchedBLongitude = [ ]
    unmatchedBLatitude  = [ ]
    
    # go through the A list and collect all the matches
    currentIndex = 0
    # look through all the A points
    for aIndex in sorted(listOfColocatedALonLat.keys()) :
        
        [aLon, aLat, aIndex, aMatches] = listOfColocatedALonLat[aIndex]
        tempMatches = 0
        
        # for each match you found on a given a point
        for matchIndex in sorted(aMatches) :
            
            [bLon, bLat, bIndex, bMatches] = listOfColocatedBLonLat[matchIndex]
            
            # copy the lon/lat info 
            matchedLongitude[currentIndex] = (aLon + bLon) / 2
            matchedLatitide[currentIndex]  = (aLat + bLat) / 2
            
            currentIndex = currentIndex + 1
            tempMatches  = tempMatches  + 1
        
        # update statistics based on the number of matches
        totalValidMatchedPairs = totalValidMatchedPairs + tempMatches
        if tempMatches > 1 :
            multipleMatchesInA = multipleMatchesInA + tempMatches
        elif tempMatches <= 0 :
            unmatchedALatitude.append(aLat)
            unmatchedALongitude.append(aLon)
    
    # gather some additional statistics from the B list
    # go through each b point
    for bIndex in sorted(listOfColocatedBLonLat) :
        
        [bLon, bLat, bIndex, bMatches] = listOfColocatedBLonLat[bIndex]
        tempMatches = len(bMatches)
        
        # update some statistics based on the number of matches
        if tempMatches > 1 :
            multipleMatchesInB = multipleMatchesInB + tempMatches
        elif tempMatches <= 0 :
            unmatchedBLatitude.append(bLat)
            unmatchedBLongitude.append(bLon)
    
    # make the unmatched lists into proper numpy arrays
    unmatchedALatitude  = np.array(unmatchedALatitude,  dtype=aLatitude.dtype)
    unmatchedALongitude = np.array(unmatchedALongitude, dtype=aLongitude.dtype)
    unmatchedBLatitude  = np.array(unmatchedBLatitude,  dtype=bLatitude.dtype)
    unmatchedBLongitude = np.array(unmatchedBLongitude, dtype=bLongitude.dtype)
    
    LOG.debug("Total matched pairs of longitude/latitide: " + str(totalValidMatchedPairs))
    
    return (matchedLongitude,    matchedLatitide, (multipleMatchesInA, multipleMatchesInB)), \
           (unmatchedALongitude, unmatchedALatitude), \
           (unmatchedBLongitude, unmatchedBLatitude)

def create_colocated_data_with_lon_lat_colocation(listOfColocatedALonLat, listOfColocatedBLonLat,
                                                  colocatedLongitude, colocatedLatitude,
                                                  aData, bData,
                                                  missingData=None, altMissingDataInB=None,
                                                  invalidAMask=None, invalidBMask=None) :
    """
    given a pre colocated list of A and B lon/lat info from create_colocation_mapping_within_epsilon,
    match up the valid data in two data sets and return the list of valid data, padded with missing
    values so that it will match the original longitude and latitude
    """
    
    if altMissingDataInB is None :
        altMissingDataInB = missingData
    
    # some general statistics
    multipleMatchesInA      = 0
    multipleMatchesInB      = 0
    totalValidMatchedPairs  = 0
    
    # our final data sets
    matchedAPoints     = np.ones(colocatedLatitude.shape, dtype=aData.dtype) * missingData
    matchedBPoints     = np.ones(colocatedLatitude.shape, dtype=bData.dtype) * altMissingDataInB
    
    # we don't know how many unmatched points we may have
    unmatchedAPoints    = [ ]
    unmatchedBPoints    = [ ]
    unmatchedALongitude = [ ]
    unmatchedALatitude  = [ ]
    unmatchedBLongitude = [ ]
    unmatchedBLatitude  = [ ]
    
    # go through the A list and sort all the valid matches
    currentIndex = 0
    # go through all the a points
    for aIndex in sorted(listOfColocatedALonLat.keys()) :
        
        [aLon, aLat, aIndex, aMatches] = listOfColocatedALonLat[aIndex]
        tempMatches = 0
        isInvalidA = invalidAMask[aIndex]
        
        # for each point that matched to a given a point
        for matchIndex in sorted(aMatches) :
            
            [bLon, bLat, bIndex, bMatches] = listOfColocatedBLonLat[matchIndex]
            isInvalidB = invalidBMask[bIndex]
            
            # if either of our data points is invalid, then the data doesn't match
            if isInvalidA or isInvalidB :
                # fill in missing data in the matches
                matchedAPoints[currentIndex] = missingData
                matchedBPoints[currentIndex] = altMissingDataInB
            else: # we have a valid match!
                tempMatches = tempMatches + 1
                matchedAPoints[currentIndex] = aData[aIndex]
                matchedBPoints[currentIndex] = bData[bIndex]
            
            currentIndex = currentIndex + 1
        
        totalValidMatchedPairs = totalValidMatchedPairs + tempMatches
        if tempMatches > 1 :
            multipleMatchesInA = multipleMatchesInA + tempMatches
        elif (tempMatches <= 0) and (not isInvalidA) :
            unmatchedAPoints.append(aData[aIndex])
            unmatchedALongitude.append(aLon)
            unmatchedALatitude.append(aLat)
    
    # gather some additional statistics from the B list
    # go through all the b points
    for bIndex in sorted(listOfColocatedBLonLat.keys()) :
        
        [bLon, bLat, bIndex, bMatches] = listOfColocatedBLonLat[bIndex]
        tempMatches = 0
        isInvalidB = invalidBMask[bIndex]
        
        # for each point that matched to a given b point
        for matchIndex in sorted(bMatches) :
            
            [aLon, aLat, aIndex, aMatches] = listOfColocatedALonLat[matchIndex]
            isInvalidA = invalidAMask[aIndex]
            
            # if either of our data points is invalid, then the data doesn't match
            if isInvalidB or isInvalidA :
                # we've already built our matched data, so no need to missing it out
                pass
            else: # we have a valid match!
                tempMatches = tempMatches + 1
        
        if tempMatches > 1 :
            multipleMatchesInB = multipleMatchesInB + tempMatches
        elif (tempMatches <= 0) and (not isInvalidB) :
            unmatchedBPoints.append(bData[bIndex])
            unmatchedBLongitude.append(bLon)
            unmatchedBLatitude.append(bLat)
    
    # make the unmatched lists into proper numpy arrays
    unmatchedAPoints    = np.array(unmatchedAPoints,    dtype=aData.dtype)
    unmatchedBPoints    = np.array(unmatchedBPoints,    dtype=bData.dtype)
    unmatchedALongitude = np.array(unmatchedALongitude, dtype=colocatedLongitude.dtype)
    unmatchedALatitude  = np.array(unmatchedALatitude,  dtype=colocatedLatitude.dtype)
    unmatchedBLongitude = np.array(unmatchedBLongitude, dtype=colocatedLongitude.dtype)
    unmatchedBLatitude  = np.array(unmatchedBLatitude,  dtype=colocatedLatitude.dtype)
    
    LOG.debug("Total matched data point pairs found: " + str(totalValidMatchedPairs))
    
    return (matchedAPoints, matchedBPoints, (multipleMatchesInA, multipleMatchesInB)), \
           (unmatchedAPoints, unmatchedALongitude, unmatchedALatitude), \
           (unmatchedBPoints, unmatchedBLongitude, unmatchedBLatitude)

if __name__=='__main__':
    import doctest
    doctest.testmod()
