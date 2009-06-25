#!/usr/bin/env python
# encoding: utf-8
"""

Top-level routines to compare two files.


Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging, re
from pprint import pprint, pformat
import numpy as np

import glance.io as io
import glance.delta as delta
import glance.plot as plot

LOG = logging.getLogger(__name__)

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
    #plotting related options
    parser.add_option('-p', '--outputpath', dest="outputpath", type='string', default='./',
                    help="set path to output directory")
    parser.add_option('-o', '--longitude', dest="longitudeVar", type='string',
                    help="set name of longitude variable")
    parser.add_option('-a', '--latitude', dest="latitudeVar", type='string',
                    help="set name of latitude variable")
    
                    
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
            lal = list(delta.summarize(aval,bval,epsilon,(amiss,bmiss)).items()) 
            # lal = list(delta.stats(*delta.diff(aval,bval,epsilon,(amiss,bmiss))).items())
            lal.sort()
            for each in lal:
                print '  %s: %s' % each
                if doc_each: print('    ' + delta.STATISTICS_DOC[each[0]])
        if doc_atend:
            print('\n\n' + delta.STATISTICS_DOC_STR)
    
    def plotDiffs(*args) :
        """create comparison images of various variables
        This option creates graphical comparisons between variables in the two given hdf files.
        Images will be created to show the variable data value for each of the two files and to show
        the difference between them.
        Variables to be compared may be specified after the names of the two input files. If no variables
        are specified, all variables that match the shape of the longitude and latitude will be compared.
        Specified variables that do not exist, do not match the correct data shape, or are the longitude/latitude
        variables will be ignored.
        Created images will be stored in the provided path, or if no path is provided, they will be stored
        in the current directory.
        The longitude and latitude variables may be specified with --longitude and --latitude
        If no longitude or latitude are specified the imager_prof_retr_abi_r4_generic1 and
        imager_prof_retr_abi_r4_generic2 variables will be used.
        Examples:
         python -m glance.compare plotDiffs A.hdf B.hdf variable_name_1 variable_name_2 variable_name_3 variable_name_4
         python -m glance.compare --outputpath=/path/where/output/will/be/placed/ plotDiffs A.hdf B.hdf
         python -m glance.compare plotDiffs --longitude=lon_variable_name --latitude=lat_variable_name A.hdf B.hdf variable_name
        """
        # get the file names the user wants to use and open them up
        aFileName, bFileName = args[:2]
        LOG.info("opening %s" % aFileName)
        aFile = io.open(aFileName)
        LOG.info("opening %s" % bFileName)
        bFile = io.open(bFileName)
        
        # get information about the variables stored in the file
        aNames = set(aFile())
        bNames = set(bFile())
        # get the variable names they have in common
        commonNames = aNames.intersection(bNames)
        # pull the ones the user asked for (if they asked for any specifically)
        requestedNames = args[2:] or ['.*']
        finalNames = _parse_varnames(commonNames, requestedNames, options.epsilon, options.missing)
        LOG.debug(str(finalNames))
        hadUserRequest = (len(args) > 2)
        
        # get the output path
        outputPath = options.outputpath
        # TODO, should I validate that the output path is a file path?
        LOG.debug(str(outputPath))
        
        # get information about the longitude/latitude data we will be using
        # to build our plots
        longitudeVariableName = options.longitudeVar or 'imager_prof_retr_abi_r4_generic1'
        latitudeVariableName = options.latitudeVar or'imager_prof_retr_abi_r4_generic2'
        LOG.debug(str("longitude variable: " + longitudeVariableName))
        LOG.debug(str("latitude variable: " + latitudeVariableName))
        # get the actual data
        longitudeA = np.array(aFile[longitudeVariableName][:], dtype=np.float)
        longitudeB = np.array(bFile[longitudeVariableName][:], dtype=np.float)
        latitudeA = np.array(aFile[latitudeVariableName][:], dtype=np.float)
        latitudeB = np.array(bFile[latitudeVariableName][:], dtype=np.float)
        # TODO validate that the two sets are the same (with .shape)
        
        # go through each of the possible variables in our files
        # and make comparison images for whichever ones we can
        for name, epsilon, missing in finalNames:
            
            # get the data for the variable
            aData = aFile[name][:]
            bData = bFile[name][:]
            
            # check if this data can be displayed and is
            # not the lon/lat variable itself
            if ((aData.shape == bData.shape) and
                (aData.shape == longitudeA.shape) and
                (bData.shape == longitudeB.shape) and
                (name != longitudeVariableName) and
                (name != latitudeVariableName)) :
                # since things match, we will try to 
                # create the images comparing that variable
                plot.plot_and_save_figure_comparison(aData, bData, name,
                                                     aFileName, bFileName,
                                                     latitudeA, longitudeA,
                                                     outputPath, epsilon,
                                                     missing)
                # TODO should I pass both sets of lat/lon?
            # only log a warning if the user themselves picked the faulty variables
            elif hadUserRequest :
                LOG.warn(name + ' could not be plotted. This may be because the data for this variable is not the ' +
                         'right shape or because the variable is currently selected as the longitude or latitude ' +
                         'variable for this file.')
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