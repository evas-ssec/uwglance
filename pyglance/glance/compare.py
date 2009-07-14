#!/usr/bin/env python
# encoding: utf-8
"""

Top-level routines to compare two files.


Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging, re, subprocess, datetime
import imp as imp
from pprint import pprint, pformat
import numpy as np

import glance.io as io
import glance.delta as delta
import glance.plot as plot
import glance.report as report

LOG = logging.getLogger(__name__)

glance_default_longitude_name = 'pixel_longitude' 
glance_default_latitude_name = 'pixel_latitude'

# these are the built in default settings
glance_analysis_defaults = {'epsilon': 0.0,
                            'missing_value': None,
                            'epsilon_failure_tolerance': None, # this means don't do the check
                            'nonfinite_data_tolerance': None # this means don't do the check
                            }

def _cvt_names(namelist, epsilon, missing):
    """"if variable names are of the format name:epsilon, yield name,epsilon, missing
        otherwise yield name,default-epsilon,default-missing
    """
    for name in namelist:
        if ':' not in name:
            yield name, epsilon
        else:
            n,e,m = name.split(':')
            if not e: e = epsilon
            else: e = float(e)
            if not m: m = missing
            else: m = float(m)
            yield n, e, m

def _parse_varnames(names, terms, epsilon=0.0, missing=None):
    """filter variable names and substitute default epsilon and missing settings if none provided
    returns name,epsilon,missing triples
    >>> _parse_varnames( ['foo','bar', 'baz', 'zoom', 'cat'], ['f..:0.5:-999', 'ba.*:0.001', 'c.t::-9999'], 1e-7 )
    set([('foo', 0.5, -999.0), ('cat', 9.9999999999999995e-08, -9999.0), ('bar', 0.001, None), ('baz', 0.001, None)])
    """
    terms = [x.split(':') for x in terms]
    terms = [(re.compile(x[0]).match,x[1:]) for x in terms]
    def _cvt_em(eps=None, mis=None):
        eps = float(eps) if eps else epsilon
        mis = float(mis) if mis else missing
        return eps, mis
    sel = [ ((x,)+_cvt_em(*em)) for x in names for (t,em) in terms if t(x) ]
    return set(sel)

def _setup_file(fileNameAndPath, prefexText='') :
    '''
    open the provided file name/path and extract information on the md5sum and last modification time
    optional prefext text may be passed in for informational output formatting
    '''
    # some info to return
    fileInfo = {'path': fileNameAndPath}
    
    # open the file
    LOG.info(prefexText + "opening " + fileNameAndPath)
    fileObject = io.open(fileNameAndPath)
    
    # get the file md5sum
    tempSubProcess = subprocess.Popen("md5sum " + fileNameAndPath, shell=True, stdout=subprocess.PIPE)
    fileInfo['md5sum'] = tempSubProcess.communicate()[0].split()[0]
    LOG.info(prefexText + "file md5sum: " + str(fileInfo['md5sum']))
    
    # get the last modified stamp
    statsForFile = os.stat(fileNameAndPath)
    fileInfo['lastModifiedTime'] = datetime.datetime.fromtimestamp(statsForFile.st_mtime).ctime() # should time zone be forced?
    LOG.info (prefexText + "file was last modified: " + fileInfo['lastModifiedTime'])
    
    return fileObject, fileInfo

def _check_file_names(fileAObject, fileBObject) :
    """
    get information about the names in the two files and how they compare to each other
    """
    # get information about the variables stored in the files
    aNames = set(fileAObject())
    bNames = set(fileBObject())
    
    # get the variable names they have in common
    commonNames = aNames.intersection(bNames)
    # which names are unique to only one of the two files?
    uniqueToANames = aNames - commonNames
    uniqueToBNames = bNames - commonNames
    
    return {'sharedVars': commonNames,  'uniqueToAVars': uniqueToANames, 'uniqueToBVars': uniqueToBNames}

def _resolve_names(fileAObject, fileBObject, defaultValues,
                   requestedNames, usingConfigFileFormat=False) :
    """
    figure out which names the two files share and which are unique to each file, as well as which names
    were requested and are in both sets
    
    usingConfigFileFormat signals whether the requestedNames parameter will be in the form of the inputed
    names from the command line or a more complex dictionary holding information about the names read in
    from a configuration file
    
    Note: if we ever need a variable with different names in file A and B to be comparable, this logic
    will need to be changed.
    """
    # look at the names present in the two files and compare them
    nameComparison = _check_file_names(fileAObject, fileBObject)
    
    # figure out which set should be selected based on the user requested names
    fileCommonNames = nameComparison['sharedVars']
    finalNames = {}
    if (usingConfigFileFormat) :
        
        # if the user didn't ask for any, try everything
        if (requestedNames == {}) :
            finalFromCommandLine = _parse_varnames(fileCommonNames, ['.*'],
                                                   defaultValues['epsilon'], defaultValues['missing_value'])
            for name, epsilon, missing in finalFromCommandLine :
                # note that we'll use the variable's name as the display name for the time being
                finalNames[name] = {'variable_name': name,
                                    'epsilon': defaultValues['epsilon'],
                                    'missing_value': defaultValues['missing_value']
                                    }
        # otherwise just do the ones the user asked for
        else : 
            # TODO, how could I do this without a loop?
            for name in fileCommonNames :
                if requestedNames.has_key(name) :
                    finalNames[name] = defaultValues.copy()
                    finalNames[name]['variable_name'] = name
                    finalNames[name].update(requestedNames[name])
    else:
        # format this similarly to the stuff from the config file
        finalFromCommandLine = _parse_varnames(fileCommonNames, requestedNames,
                                               defaultValues['epsilon'], defaultValues['missing_value'])
        for name, epsilon, missing in finalFromCommandLine :
            # note that we'll use the variable's name as the display name for the time being
            finalNames[name] = {'variable_name': name,
                                'epsilon': defaultValues['epsilon'],
                                'missing_value': defaultValues['missing_value']
                                }
    
    LOG.debug("Final selected set of variables to analyze:")
    LOG.debug(str(finalNames))
    
    return finalNames, nameComparison

def _load_config_or_options(optionsSet, originalArgs) :
    """
    load information on how the user wants to run the command either from the command line options or
    from a configuration file
    """
    
    # basic defaults for stuff we will need to return
    runInfo = {}
    runInfo['shouldIncludeReport'] = True
    runInfo['shouldIncludeImages'] = False
    runInfo['latitude'] = glance_default_latitude_name
    runInfo['longitude'] = glance_default_longitude_name
    # by default, we don't have any particular variables to analyze
    desiredVariables = {}
    # use the built in default values, to start with
    defaultsToUse = glance_analysis_defaults.copy()
    
    requestedNames = None
    
    # set up the paths, they can only come from the command line
    paths = {}
    paths['a'], paths['b'] = originalArgs[:2] # todo, let caller control # of paths expected?
    paths['out'] = optionsSet.outputpath
    
    # check to see if the user wants to use a config file and if the path exists
    requestedConfigFile = optionsSet.configFile
    usedConfigFile = False
    if (not (requestedConfigFile is None)) and os.path.exists(requestedConfigFile):
        
        LOG.info ("Using Config File Settings")
        
        # this will handle relative paths, but not '~'?
        requestedConfigFile = os.path.abspath(os.path.expanduser(requestedConfigFile))
        
        # split out the file base name and the file path
        (filePath, fileName) = os.path.split(requestedConfigFile)
        splitFileName = fileName.split('.')
        fileBaseName = fileName[:-3] # remove the '.py' from the end
        
        # load the file
        print('loading config file: ' + str(requestedConfigFile))
        glanceRunConfig = imp.load_module(fileBaseName, file(requestedConfigFile, 'U'),
                                          filePath, ('.py' , 'U', 1))
        
        # get everything from the config file
        runInfo['shouldIncludeImages'] = glanceRunConfig.shouldIncludeImages
        runInfo['latitude'] = glanceRunConfig.latitudeVar or runInfo['latitude']
        runInfo['longitude'] = glanceRunConfig.longitudeVar or runInfo['longitude']
        
        # get any requested names
        requestedNames = glanceRunConfig.setOfVariables.copy()
        
        # user selected defaults, if they omit any we'll still be using the program defaults
        defaultsToUse.update(glanceRunConfig.defaultValues)
        
        # this is an exception, since it is not advertised to the user we don't expect it to be in the file
        # (at least not at the moment, it could be added later)
        runInfo['shouldIncludeReport'] = not optionsSet.imagesOnly 
        
        usedConfigFile = True
    
    # if we didn't get the info from the config file for some reason
    # (the user didn't want to, we couldn't, etc...) get it from the command line options
    if not usedConfigFile:
        
        LOG.info ('Using Command Line Settings')
        
        # so get everything from the options directly
        runInfo['shouldIncludeReport'] = not optionsSet.imagesOnly
        runInfo['shouldIncludeImages'] = not optionsSet.htmlOnly
        runInfo['latitude'] = optionsSet.latitudeVar or runInfo['latitude']
        runInfo['longitude'] = optionsSet.longitudeVar or runInfo['longitude']
        
        # get any requested names from the command line
        requestedNames = originalArgs[2:] or ['.*']
        
        # user selected defaults
        defaultsToUse['epsilon'] = optionsSet.epsilon
        defaultsToUse['missing_value'] = optionsSet.missing
        # there is no way to set the tolerances from the command line at the moment
    
    return paths, runInfo, defaultsToUse, requestedNames, usedConfigFile

def _get_and_analyze_lon_lat (fileObject, latitudeVariableName, longitudeVariableName) :
    """
    get the longitude and latitude data from the given file, assuming they are in the given variable names
    and analyze them to identify spacially invalid data (ie. data that would fall off the earth)
    """
    # get the data from the file
    longitudeData = np.array(fileObject[longitudeVariableName][:], dtype=np.float)
    latitudeData  = np.array(fileObject[latitudeVariableName][:],  dtype=np.float)
    
    # build a mask of our spacially invalid data
    invalidLatitude = (latitudeData < -90) | (latitudeData > 90)
    invalidLongitude = (longitudeData < -180)   | (longitudeData > 180)
    spaciallyInvalidMask = invalidLatitude | invalidLongitude
    
    # analyze our spacially invalid data
    percentageOfSpaciallyInvalidPts, numberOfSpaciallyInvalidPts = _get_percentage_from_mask(spaciallyInvalidMask)
    
    return longitudeData, latitudeData, spaciallyInvalidMask, numberOfSpaciallyInvalidPts, percentageOfSpaciallyInvalidPts

def _get_percentage_from_mask(dataMask) :
    """
    given a mask that marks the elements we want the percentage of as True (and is the size of our original data),
    figure out what percentage of the whole they are
    """
    numMarkedDataPts = len(dataMask[dataMask].ravel())
    dataShape = dataMask.shape
    totalDataPts = dataShape[0] * dataShape[1]
    percentage = 100.0 * float(numMarkedDataPts) / float(totalDataPts)
    
    return percentage, numMarkedDataPts

def main():
    import optparse
    usage = """
%prog [options] 
run "%prog help" to list commands
examples:

python -m glance.compare info A.hdf
python -m glance.compare stats A.hdf B.hdf '.*_prof_retr_.*:1e-4' 'nwp_._index:0'
python -m glance.compare plotDiffs A.hdf B.hdf [optional output path]

"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--test', dest="self_test",
                    action="store_true", default=False, help="run internal unit tests")            
    parser.add_option('-q', '--quiet', dest="quiet",
                    action="store_true", default=False, help="only error output")
    parser.add_option('-v', '--verbose', dest="verbose",
                    action="store_true", default=False, help="enable more informational output")   
    parser.add_option('-w', '--debug', dest="debug",
                    action="store_true", default=False, help="enable debug output")   
    parser.add_option('-e', '--epsilon', dest="epsilon", type='float', default=0.0,
                    help="set default epsilon value for comparison threshold")   
    parser.add_option('-m', '--missing', dest="missing", type='float', default=None,
                    help="set default missing-value")
    #report generation related options
    parser.add_option('-p', '--outputpath', dest="outputpath", type='string', default='./',
                    help="set path to output directory")
    parser.add_option('-o', '--longitude', dest="longitudeVar", type='string',
                    help="set name of longitude variable")
    parser.add_option('-a', '--latitude', dest="latitudeVar", type='string',
                    help="set name of latitude variable")
    parser.add_option('-i', '--imagesonly', dest="imagesOnly", 
                      action="store_true", default=False,
                      help="generate only image files (no html report)")
    parser.add_option('-r', '--reportonly', dest="htmlOnly", 
                      action="store_true", default=False,
                      help="generate only html report files (no images)")
    parser.add_option('-c', '--configfile', dest="configFile", type='string', default=None,
                      help="set optional configuration file")
    
                    
    options, args = parser.parse_args()
    if options.self_test:
        import doctest
        doctest.testmod()
        sys.exit(2)

    lvl = logging.WARNING
    if options.debug: lvl = logging.DEBUG
    elif options.verbose: lvl = logging.INFO
    elif options.quiet: lvl = logging.ERROR
    logging.basicConfig(level = lvl)

    commands = {}
    prior = None
    prior = dict(locals())
    
    def info(*args):
        """list information about a list of files
        List available variables for comparison.
        """
        for fn in args:
            lal = list(io.open(fn)())
            lal.sort()
            print fn + ': ' + ('\n  ' + ' '*len(fn)).join(lal)
    
    def sdr_cris(*args):
        """compare sdr_cris output
        parameters are variable name followed by detector number
        sdr_cris desired.h5 actual.h5 ESRealLW 0
        """
        afn,bfn = args[:2]
        LOG.info("opening %s" % afn)
        a = io.open(afn)
        LOG.info("opening %s" % bfn)
        b = io.open(bfn)
        # shape is [scanline, field, detector, wnum]
        vname = '/All_Data/CrIS-SDR_All/' + args[2]
        det_idx = int(args[3])
        def get(f):
            spc = f[vname][:,:,det_idx,:]
            nsl,nfor,nwn = spc.shape
            return spc.reshape( (nsl*nfor,nwn) )
        aspc = get(a)
        bspc = get(b)
        plot.compare_spectra(bspc,aspc)
        plot.show()
    
    def noisecheck(*args):
        """gives statistics for dataset comparisons against truth with and without noise
        usage: noisecheck truth-file noise-file actual-file variable1{:epsilon{:missing}} {variable2...}
        glance noisecheck /Volumes/snaapy/data/justins/abi_graffir/coreg/pure/l2_data/geocatL2.GOES-R.2005155.220000.hdf.gz /Volumes/snaapy/data/justins/abi_graffir/noise/noise1x/l2_data/geocatL2.GOES-R.2005155.220000.hdf 
        """
        afn,noizfn,bfn = args[:3]
        LOG.info("opening truth file %s" % afn)
        a = io.open(afn)
        LOG.info("opening actual file %s" % noizfn)
        noiz = io.open(noizfn)
        LOG.info("opening noise file %s" % bfn)
        b = io.open(bfn)
        
        anames = set(a())
        bnames = set(b()) 
        cnames = anames.intersection(bnames) # common names
        pats = args[3:] or ['.*']
        names = _parse_varnames( cnames, pats, options.epsilon, options.missing )
        for name,epsilon,missing in names:
            avar = a[name]
            bvar = b[name]
            nvar = noiz[name]
            if missing is None:
                amiss = a.missing_value(name)
                bmiss = b.missing_value(name)
            else:
                amiss,bmiss = missing,missing
            x = avar[:]
            y = bvar[:]
            z = nvar[:]
            def scat(x,xn,y):
                from pylab import plot,show,scatter
                scatter(x,y)
                show()
            nfo = delta.rms_corr_withnoise(x,y,z,epsilon,(amiss,bmiss),plot=scat)
            print '-'*32
            print name
            for kv in sorted(nfo.items()):
                print '  %s: %s' % kv
    
    def stats(*args):
        """create statistics summary of variables
        Summarize difference statistics between listed variables.
        If no variable names are given, summarize all common variables.
        Variable names can be of the form varname:epsilon:missing to use non-default epsilon or missing value.
        Variable names can be regular expressions, e.g. 'image.*' or '.*prof_retr.*::-999'
        Either epsilon or missing can be empty to stay with default.
        If _FillValue is an attribute of a variable, that will be used to find missing values where no value is given.
        Run with -v to get more detailed information on statistics.
        Examples:
         python -m glance.compare stats hdffile1 hdffile2
         python -m glance.compare stats --epsilon=0.00001 A.hdf B.hdf baseline_cmask_seviri_cloud_mask:0.002:
         python -m glance.compare -w stats --epsilon=0.00001 A.hdf A.hdf imager_prof_retr_abi_total_precipitable_water_low::-999
        """ 
        afn,bfn = args[:2]
        LOG.info("opening %s" % afn)
        a = io.open(afn)
        LOG.info("opening %s" % bfn)
        b = io.open(bfn)
        anames = set(a())
        bnames = set(b()) 
        cnames = anames.intersection(bnames) # common names
        pats = args[2:] or ['.*']
        names = _parse_varnames( cnames, pats, options.epsilon, options.missing )
        LOG.debug(str(names))
        doc_each = (options.verbose or options.debug) and len(names)==1
        doc_atend = (options.verbose or options.debug) and len(names)!=1
        for name,epsilon,missing in names:
            avar = a[name]
            bvar = b[name]
            if missing is None:
                amiss = a.missing_value(name)
                bmiss = b.missing_value(name)
            else:
                amiss,bmiss = missing,missing
            LOG.debug('comparing %s with epsilon %s and missing %s,%s' % (name,epsilon,amiss,bmiss))
            aval = avar[:]
            bval = bvar[:]
            print '-'*32
            print name
            print 
            lal = list(delta.summarize(aval,bval,epsilon,(amiss,bmiss)).items()) 
            # lal = list(delta.stats(*delta.diff(aval,bval,epsilon,(amiss,bmiss))).items())
            lal.sort()
            for dictionary_title, dict_data in lal:
                print '%s' %  dictionary_title
                dict_data
                for each_stat in sorted(list(dict_data)):
                    print '  %s: %s' % (each_stat, dict_data[each_stat])
                    if doc_each: print('    ' + delta.STATISTICS_DOC[each_stat])
                print 
        if doc_atend:
            print('\n\n' + delta.STATISTICS_DOC_STR)

    def plotDiffs(*args) :
        """generate a set of images comparing two files
        This option creates a set of graphical comparisons of variables in the two given hdf files.
        The images detailing the differences between variables in the two hdf files will be
        generated and saved to disk. 
        Variables to be compared may be specified after the names of the two input files. If no variables
        are specified, all variables that match the shape of the longitude and latitude will be compared.
        Specified variables that do not exist, do not match the correct data shape, or are the longitude/latitude
        variables will be ignored.
        The user may also use the notation variable_name:epsilon:missing_value to specify the acceptible epsilon
        for comparison and the missing_value which indicates missing data. If one or both of these values is absent
        (in the case of variable_name:epsilon: variable_name::missing_value or just variable_name) the default value
        of 0.0 will be used for epsilon and no missing values will be analyzed. 
        The created images will be stored in the provided path, or if no path is provided, they will be stored in
        the current directory.
        The longitude and latitude variables may be specified with --longitude and --latitude
        If no longitude or latitude are specified the pixel_latitude and pixel_longitude variables will be used.
        Examples:
         python -m glance.compare plotDiffs A.hdf B.hdf variable_name_1:epsilon1: variable_name_2 variable_name_3:epsilon3:missing3 variable_name_4::missing4
         python -m glance.compare --outputpath=/path/where/output/will/be/placed/ plotDiffs A.hdf B.hdf
         python -m glance.compare plotDiffs --longitude=lon_variable_name --latitude=lat_variable_name A.hdf B.hdf variable_name
        """
        # set the options so that a report will not be generated
        options.imagesOnly = True
        
        # make the images
        reportGen(*args)
        
        return

    def reportGen(*args) :
        """generate a report comparing two files
        This option creates a report comparing variables in the two given hdf files.
        An html report and images detailing the differences between variables in the two hdf files will be
        generated and saved to disk. The images will be embedded in the report or visible as separate .png files.
        Variables to be compared may be specified after the names of the two input files. If no variables
        are specified, all variables that match the shape of the longitude and latitude will be compared.
        Specified variables that do not exist, do not match the correct data shape, or are the longitude/latitude
        variables will be ignored.
        The user may also use the notation variable_name:epsilon:missing_value to specify the acceptible epsilon
        for comparison and the missing_value which indicates missing data. If one or both of these values is absent
        (in the case of variable_name:epsilon: variable_name::missing_value or just variable_name) the default value
        of 0.0 will be used for epsilon and no missing values will be analyzed. 
        The html report page(s) and any created images will be stored in the provided path, or if no path is provided,
        they will be stored in the current directory.
        If for some reason you would prefer to generate the report without images, use the --reportonly option. This
        option will generate the html report but omit the images. This may be significantly faster, depending on
        your system, but the differences between the files may be quite a bit more difficult to interpret.
        The longitude and latitude variables may be specified with --longitude and --latitude
        If no longitude or latitude are specified the pixel_latitude and pixel_longitude variables will be used.
        Examples:
         python -m glance.compare reportGen A.hdf B.hdf variable_name_1:epsilon1: variable_name_2 variable_name_3:epsilon3:missing3 variable_name_4::missing4
         python -m glance.compare --outputpath=/path/where/output/will/be/placed/ reportGen A.hdf B.hdf
         python -m glance.compare reportGen --longitude=lon_variable_name --latitude=lat_variable_name A.hdf B.hdf variable_name
         python -m glance.compare reportGen --imagesonly A.hdf B.hdf
        """
        
        # load the user settings from either the command line or a user defined config file
        pathsTemp, runInfo, defaultValues, requestedNames, usedConfigFile = _load_config_or_options(options, args)
        
        # note some of this information for debugging purposes
        LOG.debug('paths: ' + str(pathsTemp))
        LOG.debug('defaults: ' + str(defaultValues))
        LOG.debug(str("longitude variable: " + runInfo['longitude']))
        LOG.debug(str("latitude variable: " + runInfo['latitude']))
        
        # if we wouldn't generate anything, just stop now
        if (not runInfo['shouldIncludeImages']) and (not runInfo['shouldIncludeReport']) :
            LOG.warn("User selection of no image generation and no report generation will result in no " +
                     "content being generated. Aborting report generation function.")
            return
        
        # get info on who's doing the run and where
        runInfo['machine'] = os.uname()[1] # the name of the machine running the report
        runInfo['user'] = os.getlogin() # the name of the user running the report
        
        # get information about the input and output files
        outputPath = pathsTemp['out']
        LOG.debug("output path: " + str(outputPath))
        # open the files
        files = {}
        LOG.info("Processing File A:")
        aFile, files['file A'] = _setup_file(pathsTemp['a'], "\t")
        LOG.info("Processing File B:")
        bFile, files['file B'] = _setup_file(pathsTemp['b'], "\t")
        
        # get information about the names the user requested
        hadUserRequest = (not (requestedNames is None)) and (len(requestedNames) > 0)
        finalNames, nameStats = _resolve_names(aFile, bFile,
                                               defaultValues,
                                               requestedNames, usedConfigFile)
        
        # get and analyze our longitude and latitude data
        spatialInfo = {'file A': {}, 'file B':{}}
        longitudeA, latitudeA, spaciallyInvalidMaskA, \
            spatialInfo['file A']['totNumInvPts'], spatialInfo['file A']['perInvPts'] = \
            _get_and_analyze_lon_lat (aFile, runInfo['latitude'], runInfo['longitude'])
        longitudeB, latitudeB, spaciallyInvalidMaskB, \
            spatialInfo['file B']['totNumInvPts'], spatialInfo['file B']['perInvPts'] = \
            _get_and_analyze_lon_lat (bFile, runInfo['latitude'], runInfo['longitude'])
        
        # test the "valid" values in our lon/lat
        longitude = longitudeA
        latitude = latitudeA
        if not all(longitudeA.ravel() == longitudeB.ravel()) and all(latitudeA.ravel() == latitudeB.ravel()) :
            LOG.warn("Possible mismatch in values stored in file a and file b longitude and latitude values."
                     + " Depending on the degree of mismatch, some data value comparisons may be "
                     + "distorted or spacially nonsensical.")
            # TODO, this should cause a warning to be put in the report
        
        # compare our spacialy invalid info
        spaciallyInvalidMask = spaciallyInvalidMaskA | spaciallyInvalidMaskB
        spatialInfo['perInvPtsInBoth'] = spatialInfo['file A']['perInvPts']
                # a default that will hold if the two files have the same spatially invalid pts
        
        if not all(spaciallyInvalidMaskA.ravel() == spaciallyInvalidMaskB.ravel()) : 
            LOG.info("Mismatch in number of spatially invalid points. " +
                     "Files may not have corresponding data where expected.")
            
            # figure out which points are only valid in one of the two files
            validOnlyInAMask = (~spaciallyInvalidMaskA) & spaciallyInvalidMaskB
            spatialInfo['file A']['numInvPts'] = len(validOnlyInAMask[validOnlyInAMask])
            validOnlyInBMask = (~spaciallyInvalidMaskB) & spaciallyInvalidMaskA
            spatialInfo['file B']['numInvPts'] = len(validOnlyInBMask[validOnlyInBMask])
            
            # so how many do they have together?
            spatialInfo['perInvPtsInBoth'], totalNumSpaciallyInvPts = _get_percentage_from_mask(spaciallyInvalidMask)
            # make a "clean" version of the lon/lat
            longitude[validOnlyInAMask] = longitudeA[validOnlyInAMask]
            longitude[validOnlyInBMask] = longitudeB[validOnlyInBMask]
            latitude [validOnlyInAMask] = latitudeA [validOnlyInAMask]
            latitude [validOnlyInBMask] = latitudeB [validOnlyInBMask]
            
            # plot the points that are only valid one file and not the other
            if (spatialInfo['file A']['numInvPts'] > 0) and (runInfo['shouldIncludeImages']) :
                plot.plot_and_save_spacial_trouble(longitude, latitude,
                                                   validOnlyInAMask,
                                                   spaciallyInvalidMaskA,
                                                   "A", outputPath, True)
            if (spatialInfo['file B']['numInvPts'] > 0) and (runInfo['shouldIncludeImages']) :
                plot.plot_and_save_spacial_trouble(longitude, latitude,
                                                   validOnlyInBMask,
                                                   spaciallyInvalidMaskB,
                                                   "B", outputPath, True)
            
        # set some things up to hold info for our reports
        # this is going to be in the form
        # [var_name] = {"passEpsilonPercent": percent ok with epsilon, "epsilon": epsilon)
        variableComparisons = {}
        
        # go through each of the possible variables in our files
        # and make a report section with images for whichever ones we can
        for name in finalNames:
            
            # pull out the information for this variable analysis run
            varRunInfo = finalNames[name].copy()
            displayName = name
            if (varRunInfo.has_key('display_name')) :
                displayName = varRunInfo['display_name']
            
            # get the data for the variable
            aData = aFile[varRunInfo['variable_name']][:]
            bData = bFile[varRunInfo['variable_name']][:]
            
            # check if this data can be displayed
            if ((aData.shape == bData.shape) and
                (aData.shape == longitude.shape) and
                (bData.shape == longitude.shape)) :
                
                # if we should be making images, then make them for this variable
                if (runInfo['shouldIncludeImages']) :
                    # create the images comparing that variable
                    plot.plot_and_save_figure_comparison(aData, bData, varRunInfo['variable_name'],
                                                         files['file A']['path'], files['file B']['path'],
                                                         latitude, longitude,
                                                         spaciallyInvalidMaskA,
                                                         spaciallyInvalidMaskB,
                                                         spaciallyInvalidMask,
                                                         outputPath, varRunInfo['epsilon'],
                                                         varRunInfo['missing_value'], True,
                                                         displayName)
                
                # generate the report for this variable
                if (runInfo['shouldIncludeReport']) :
                    # get the current time
                    runInfo['time'] = datetime.datetime.ctime(datetime.datetime.now())
                    #get info on the variable
                    variableStats = delta.summarize(aData, bData, varRunInfo['epsilon'],
                                                    (varRunInfo['missing_value'], varRunInfo['missing_value']),
                                                    spaciallyInvalidMask)
                    # hang on to our good % and our epsilon value to describe our comparison
                    variableComparisons[varRunInfo['variable_name']] = \
                            {'pass_epsilon_percent': ((1.0 -
                                                       variableStats['Numerical Comparison Statistics']['diff_outside_epsilon_fraction']) * 100.0),
                             'variable_run_info': varRunInfo
                             }
                    print ('generating report for: ' + varRunInfo['variable_name'])
                    report.generate_and_save_variable_report(files,
                                                             varRunInfo, runInfo,
                                                             variableStats, 
                                                             outputPath, varRunInfo['variable_name'] + ".html")
                                                             
                                                             
                    
            # only log a warning if the user themselves picked the faulty variable
            elif hadUserRequest :
                LOG.warn(varRunInfo['variable_name'] + ' ' +
                         'could not be compared. This may be because the data for this variable does not match in shape ' +
                         'between the two files or the data may not match the shape of the selected longitude and ' +
                         'latitude variables.')
        
        # generate our general report pages once we've looked at all the variables
        if (runInfo['shouldIncludeReport']) :
            print ('generating summary report')
            # get the current time
            runInfo['time'] = datetime.datetime.ctime(datetime.datetime.now())
            # generate the report summary page
            report.generate_and_save_summary_report(files,
                                                    outputPath, 'index.html',
                                                    runInfo,
                                                    variableComparisons, 
                                                    spatialInfo,
                                                    nameStats)
            # make the glossary
            print ('generating glossary')
            report.generate_and_save_doc_page(delta.STATISTICS_DOC, outputPath)
        
        return

    # def build(*args):
    #     """build summary
    #     build extended info
    #     """
    #     LOG.info("building database tables")
    #     
    # def grant(*args):
    #     """grant summary
    #     grant extended info
    #     """
    #     LOG.info("granting permissions for tables")
    #     
    # def index(*args):
    #     """index summary
    #     index extended info
    #     """
    #     LOG.info("creating indices for tables")
        
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
    commands.update(dict(x for x in locals().items() if x[0] not in prior))    
    
    if (not args) or (args[0] not in commands): 
        parser.print_help()
        help()
        return 9
    else:
        locals()[args[0]](*args[1:])

    return 0


if __name__=='__main__':
    sys.exit(main())