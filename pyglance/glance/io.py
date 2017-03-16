#!/usr/bin/env python
# encoding: utf-8
"""
I/O routines supporting reading a number of file formats.

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, logging
import numpy as np

LOG = logging.getLogger(__name__)

try:
    import pyhdf
    from pyhdf.SD import SD,SDC, SDS, HDF4Error
except:
    LOG.info('no pyhdf module available for HDF4')
    pyhdf = None
    SD = SDC = SDS = object
    HDF4Error = EnvironmentError
    
try:
    import h5py
    from h5py import h5d
except ImportError:
    LOG.info('no h5py module available for reading HDF5')
    h5py = None

# the newer netCDF library that replaced pycdf
try:
    import netCDF4
except:
    LOG.info("unable to import netcdf4 library")
    netCDF4 = None

""" this is the previous netCDF library, remove this once the new one is fully tested
try:    
    import pycdf
    from pycdf import CDF, NC, strerror
except:
    LOG.info('no pycdf module available')
    pycdf = None
    CDF = NC = object
    def strerror(*args):
        return 'no pycdf module installed'
"""

try:
    import dmv as dmvlib
    LOG.info('loaded dmv module for AERI data file access')
except ImportError:
    LOG.info('no AERI dmv data file format module')
    dmvlib = None

try:
    import adl_blob
    LOG.info('adl_blob module found for JPSS ADL data file access')
except ImportError:
    LOG.info('no adl_blob format handler available')
    adl_blob = None

try :
    from osgeo import gdal
    LOG.info('loading osgeo module for GeoTIFF data file access')
except :
    LOG.info('no osgeo available for reading GeoTIFF data files')
    gdal = None

UNITS_CONSTANT = "units"

fillValConst1 = '_FillValue'
fillValConst2 = 'missing_value'

ADD_OFFSET_STR   = 'add_offset'
SCALE_FACTOR_STR = 'scale_factor'
SCALE_METHOD_STR = 'scaling_method'

UNSIGNED_ATTR_STR = "_unsigned"

SIGNED_TO_UNSIGNED_DTYPES = {
                                np.dtype(np.int8):   np.dtype(np.uint8),
                                np.dtype(np.int16):   np.dtype(np.uint16),
                                np.dtype(np.int32):   np.dtype(np.uint32),
                                np.dtype(np.int64):   np.dtype(np.uint64),
                            }

class IOUnimplimentedError(Exception):
    """
    The exception raised when a requested io operation is not yet available.
    
        msg  -- explanation of the problem
    """
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

class CaseInsensitiveAttributeCache (object) :
    """
    A cache of attributes for a single file and all of it's variables.
    This cache is considered uncased, it will store all attributes it caches
    in lower case and will lower case any strings it is asked to search for
    in the cache.
    When variable or global attribute sets are not yet loaded and something
    from that part of the file is requested the cache will transparently load
    attributes from the file behind the scenes and build the cache for that
    part of the file.
    """
    
    def __init__(self, fileObject) :
        """
        set up the empty cache and hang on to the file object we'll be caching
        """
        
        self.fileToCache             = fileObject
        self.globalAttributesLower   = None
        self.variableAttributesLower = { }
    
    def _load_global_attributes_if_needed (self) :
        """
        load up the global attributes if they need to be cached
        """
        
        # load the attributes from the file if they aren't cached
        if self.globalAttributesLower is None :
            LOG.debug ("Loading file global attributes into case-insensitive cache.")
            tempAttrs                  = self.fileToCache.get_global_attributes(caseInsensitive=False)
            self.globalAttributesLower = dict((k.lower(), v) for k, v in tempAttrs.items())
    
    def _load_variable_attributes_if_needed (self, variableName) :
        """
        load up the variable attributes if they need to be cached
        """
        
        # make a lower cased version of the variable name
        tempVariableName = variableName.lower()
        
        # load the variable's attributes from the file if they aren't cached
        if tempVariableName not in self.variableAttributesLower.keys() :
            LOG.debug ("Loading attributes for variable \"" + variableName + "\" into case-insensitive cache.")
            tempAttrs = self.fileToCache.get_variable_attributes(variableName, caseInsensitive=False)
            # now if there are any attributes, make a case insensitive version
            self.variableAttributesLower[tempVariableName] = dict((k.lower(), v) for k, v in tempAttrs.items())
    
    def get_variable_attribute (self, variableName, attributeName) :
        """
        get the specified attribute for the specified variable,
        if this variable's attributes have not yet been loaded
        they will be loaded and cached
        """
        
        self._load_variable_attributes_if_needed(variableName)
        
        toReturn = None
        tempVariableName  =  variableName.lower()
        tempAttributeName = attributeName.lower()
        if (tempVariableName in self.variableAttributesLower) and (tempAttributeName in self.variableAttributesLower[tempVariableName]) :
            toReturn = self.variableAttributesLower[tempVariableName][tempAttributeName]
        else:
            LOG.debug ("Attribute \"" + attributeName + "\" was not present for variable \"" + variableName + "\".")
        
        return toReturn
    
    def get_variable_attributes (self, variableName) :
        """
        get the variable attributes for the variable name given
        """
        
        self._load_variable_attributes_if_needed(variableName)
        
        toReturn = self.variableAttributesLower[variableName.lower()] if (variableName.lower() in self.variableAttributesLower) else None
        
        return toReturn
    
    def get_global_attribute (self, attributeName) :
        """
        get a global attribute with the given name
        """
        
        self._load_global_attributes_if_needed()
        
        toReturn = self.globalAttributesLower[attributeName.lower()] if (attributeName.lower() in self.globalAttributesLower) else None
        
        return toReturn
    
    def get_global_attributes (self) :
        """
        get the global attributes,
        """
        
        self._load_global_attributes_if_needed()
        
        toReturn = self.globalAttributesLower
        
        return toReturn
    
    def is_loadable_type (self, name) :
        """
        check to see if the indicated variable is a type that can be loaded
        """
        
        # TODO, are there any bad types for these files?
        return True

class hdf (object):
    """wrapper for HDF4 dataset for comparison
    __call__ yields sequence of variable names
    __getitem__ returns individual variables ready for slicing to numpy arrays
    """
    
    _hdf = None
    
    def __init__(self, filename, allowWrite=False):
        
        if pyhdf is None:
            LOG.error('pyhdf is not installed and is needed in order to read hdf4 files')
            assert(pyhdf is not None)
        mode = SDC.READ
        if allowWrite:
            mode = mode | SDC.WRITE
        
        self._hdf = SD(filename, mode)
        self.attributeCache = CaseInsensitiveAttributeCache(self)

    def __call__(self):
        "yield names of variables to be compared"
        return self._hdf.datasets().keys()
    
    # this returns a numpy array with a copy of the full, scaled
    # data for this variable, if the data type must be changed to allow
    # for scaling it will be (so the return type may not reflect the
    # type found in the original file)
    def __getitem__(self, name):
        # defaults
        scale_factor = 1.0
        add_offset = 0.0
        data_type = None 
        scaling_method = None
        
        # get the variable object and use it to
        # get our raw data and scaling info
        variable_object = self.get_variable_object(name)
        raw_data_copy = variable_object[:]
        try :
            # TODO, this currently won't work with geocat data, work around it for now
            scale_factor, scale_factor_error, add_offset, add_offset_error, data_type = SDS.getcal(variable_object)
        except HDF4Error:
            # load just the scale factor and add offset information by hand
            temp = self.attributeCache.get_variable_attributes(name)
            if ADD_OFFSET_STR in temp.keys() :
                add_offset = temp[ADD_OFFSET_STR]
                data_type = np.dtype(type(add_offset))
            if SCALE_FACTOR_STR in temp.keys() :
                scale_factor = temp[SCALE_FACTOR_STR]
                data_type = np.dtype(type(scale_factor))
            if SCALE_METHOD_STR in temp.keys() :
                scaling_method = temp[SCALE_METHOD_STR]
        SDS.endaccess(variable_object)
        
        # don't do lots of work if we don't need to scale things
        if (scale_factor == 1.0) and (add_offset == 0.0) :
            return raw_data_copy
        
        # at the moment geocat has several scaling methods that don't match the normal standards for hdf
        """
        please see constant.f90 for a more up to date version of this information:
            INTEGER(kind=int1) :: NO_SCALE              ! 0
            INTEGER(kind=int1) :: LINEAR_SCALE          ! 1
            INTEGER(kind=int1) :: LOG_SCALE             ! 2
            INTEGER(kind=int1) :: SQRT_SCALE            ! 3 
        """
        if (scaling_method == 0) :
            return raw_data_copy
        if not ((scaling_method is None) or (int(scaling_method) <= 1)) :
            LOG.warn ('Scaling method of \"' + str(scaling_method) + '\" will be ignored in favor of hdf standard method. '
                      + 'This may cause problems with data consistency')
        
        # if we don't have a data type something strange has gone wrong
        assert(not (data_type is None))
        
        # get information about where the data is the missing value
        missing_val = self.missing_value(name)
        missing_mask = np.zeros(raw_data_copy.shape, dtype=np.bool)
        missing_mask[raw_data_copy == missing_val] = True
        
        # create the scaled version of the data
        scaled_data_copy                = np.array(raw_data_copy, dtype=data_type)
        scaled_data_copy[~missing_mask] = (scaled_data_copy[~missing_mask] * scale_factor) + add_offset #TODO, type truncation issues?
        
        return scaled_data_copy 
    
    def get_variable_object(self, name):
        return self._hdf.select(name)
    
    def missing_value(self, name):
        
        return self.get_attribute(name, fillValConst1)
    
    def create_new_variable(self, variablename, missingvalue=None, data=None, variabletocopyattributesfrom=None):
        """
        create a new variable with the given name
        optionally set the missing value (fill value) and data to those given
        
        the created variable will be returned, or None if a variable could not
        be created
        """
        
        raise IOUnimplimentedError('Unable to create variable in hdf file, this functionality is not yet available.')
        
        return None
    
    def add_attribute_data_to_variable(self, variableName, newAttributeName, newAttributeValue) :
        """
        if the attribute exists for the given variable, set it to the new value
        if the attribute does not exist for the given variable, create it and set it to the new value
        """
        
        raise IOUnimplimentedError('Unable add attribute to hdf file, this functionality is not yet available.')
        
        return
    
    def get_variable_attributes (self, variableName, caseInsensitive=True) :
        """
        returns all the attributes associated with a variable name
        """
        
        toReturn = None
        if caseInsensitive :
            toReturn = self.attributeCache.get_variable_attributes(variableName)
        else :
            toReturn = self.get_variable_object(variableName).attributes()
        
        return toReturn
    
    def get_attribute(self, variableName, attributeName, caseInsensitive=True) :
        """
        returns the value of the attribute if it is available for this variable, or None
        """
        toReturn = None
        
        if caseInsensitive :
            toReturn = self.attributeCache.get_variable_attribute(variableName, attributeName)
        else :
            temp_attributes = self.get_variable_attributes(variableName, caseInsensitive=False)
            
            if attributeName in temp_attributes :
                toReturn = temp_attributes[attributeName]
        
        return toReturn
    
    def get_global_attributes(self, caseInsensitive=True) :
        """
        get a list of all the global attributes for this file or None
        """
        
        toReturn = None
        
        if caseInsensitive :
            self.attributeCache.get_global_attributes()
        else :
            toReturn = self._hdf.attributes()
        
        return toReturn
    
    def get_global_attribute(self, attributeName, caseInsensitive=True) :
        """
        returns the value of a global attribute if it is available or None
        """
        
        toReturn = None
        
        if caseInsensitive :
            toReturn = self.attributeCache.get_global_attribute(attributeName)
        else :
            if attributeName in self._hdf.attributes() :
                toReturn = self._hdf.attributes()[attributeName]
        
        return toReturn
    
    def is_loadable_type (self, name) :
        """
        check to see if the indicated variable is a type that can be loaded
        """
        
        # TODO, are there any bad types for these files?
        return True

class nc (object):
    """wrapper for netcdf4-python data access for comparison
    __call__ yields sequence of variable names
    __getitem__ returns individual variables ready for slicing to numpy arrays
    """
    
    _nc = None
    
    def __init__(self, filename, allowWrite=False):
        
        if netCDF4 is None:
            LOG.error('netCDF4 is not installed and is needed in order to read NetCDF files')
            assert(netCDF4 is not None)
        
        mode = 'r'
        if allowWrite :
            mode = 'w'
        
        self._nc = netCDF4.Dataset(filename, mode)
        self.attributeCache = CaseInsensitiveAttributeCache(self)

    def __call__(self):
        "yield names of variables to be compared"
        return self._nc.variables.keys()
    
    # this returns a numpy array with a copy of the full, scaled
    # data for this variable, if the data type must be changed to allow
    # for scaling it will be (so the return type may not reflect the
    # type found in the original file)
    def __getitem__(self, name):
        
        #print ("*** opening variable: " + name)
        
        # defaults
        scale_factor = 1.0
        add_offset = 0.0
        data_type = np.float32 # TODO temporary
        
        # get the variable object and use it to
        # get our raw data and scaling info
        variable_object = self.get_variable_object(name)

        """ # This scaling code is no longer required because the library automatically handles scaling
        raw_data_copy = variable_object[:]
        # load the scale factor and add offset
        temp = self.attributeCache.get_variable_attributes(name)
        if SCALE_FACTOR_STR in temp.keys() :
            scale_factor = temp[SCALE_FACTOR_STR]
        if ADD_OFFSET_STR in temp.keys() :
            add_offset = temp[ADD_OFFSET_STR]
        
        # don't do lots of work if we don't need to scale things
        if (scale_factor == 1.0) and (add_offset == 0.0) :
            return raw_data_copy
        
        # get information about where the data is the missing value
        missing_val = self.missing_value(name)
        missing_mask = np.zeros(raw_data_copy.shape, dtype=np.bool)
        missing_mask[raw_data_copy == missing_val] = True

        # create the scaled version of the data
        scaled_data_copy = np.array(raw_data_copy, dtype=data_type)
        scaled_data_copy[~missing_mask] = (scaled_data_copy[~missing_mask] * scale_factor) + add_offset #TODO, type truncation issues?
        """

        # get our data, save the dtype, and make sure it's a more flexible dtype for now
        scaled_data_copy = np.array(variable_object[:], dtype=data_type)

        temp = self.attributeCache.get_variable_attributes(name)
        if UNSIGNED_ATTR_STR in temp.keys() and str(temp[UNSIGNED_ATTR_STR]).lower() == ( "true" ) :

            LOG.debug("fixing unsigned values in variable " + name)

            # load the scale factor and add offset
            scale_factor = 1.0
            add_offset = 0.0
            temp = self.attributeCache.get_variable_attributes(name)
            if SCALE_FACTOR_STR in temp.keys() :
                scale_factor = temp[SCALE_FACTOR_STR]
            if ADD_OFFSET_STR in temp.keys() :
                add_offset = temp[ADD_OFFSET_STR]

            # get the missing value and figure out the dtype of the original data
            missing_val  = self.missing_value(name)
            orig_dtype   = np.array([missing_val,]).dtype
            needed_dtype = SIGNED_TO_UNSIGNED_DTYPES[orig_dtype] if orig_dtype in SIGNED_TO_UNSIGNED_DTYPES.keys() else None

            if needed_dtype is not None :
                # now figure out where all the corrupted values are, and shift them up to be positive
                needs_fix_mask = (scaled_data_copy < add_offset) & (scaled_data_copy != missing_val)
                # we are adding the 2's complement, but first we're scaling it appropriately
                scaled_data_copy[needs_fix_mask] += ((np.iinfo(np.uint16).max + 1.0) * scale_factor)

        return scaled_data_copy
    
    # TODO, this hasn't been supported in other file types
    def close (self) :
        self._nc.close()
        self._nc = None

    def get_variable_object(self, name):

        return self._nc.variables[name]
    
    def missing_value(self, name):
        
        toReturn = None
        
        temp = self.attributeCache.get_variable_attribute(name, fillValConst1)
        if temp is not None :
            toReturn = temp
        else :
            temp = self.attributeCache.get_variable_attribute(name, fillValConst2)
            if temp is not None :
                toReturn = temp
        
        return toReturn

    def create_new_variable(self, variablename, missingvalue=None, data=None, variabletocopyattributesfrom=None):
        """
        create a new variable with the given name
        optionally set the missing value (fill value) and data to those given
        
        the created variable will be returned, or None if a variable could not
        be created
        """
        
        self._nc.nc_redef()
        
        # if the variable already exists, stop with a warning
        if variablename in self._nc.variables.keys() :
            LOG.warn("New variable name requested (" + variablename + ") is already present in file. " +
                     "Skipping generation of new variable.")
            return None
        # if we have no data we won't be able to determine the data type to create the variable
        if (data is None) or (len(data) <= 0) :
            LOG.warn("Data type for new variable (" + variablename + ") could not be determined. " +
                     "Skipping generation of new variable.")
            return None
        
        dataType = None
        if np.issubdtype(data.dtype, int) :
            dataType = np.int
            #print("Picked INT")
        # TODO, at the moment the fill type is forcing me to use a double, when sometimes I want a float
        #elif np.issubdtype(data.dtype, np.float32) :
        #    dataType = np.float
        #    print("Picked FLOAT")
        elif np.issubdtype(data.dtype, float) :
            dataType = np.float64
            #print("Picked DOUBLE")
        # what do we do if it's some other type?
        
        # create and set all the dimensions
        dimensions = [ ]
        dimensionNum = 0
        for dimSize in data.shape :
            dimensions.append(self._nc.createDimension(variablename + '-index' + str(dimensionNum), dimSize))
            dimensionNum = dimensionNum + 1
        
        # create the new variable
        #print('variable name: ' + variablename)
        #print('data type:     ' + str(dataType))
        #print('dimensions:    ' + str(dimensions))
        newVariable = self._nc.createVariable(variablename, dataType, tuple(dimensions))
        
        # if a missing value was given, use that
        if missingvalue is not None :
            newVariable._FillValue = missingvalue
        
        # if we have a variable to copy attributes from, do so
        if variabletocopyattributesfrom is not None :
            attributes = self.get_variable_attributes(variabletocopyattributesfrom, caseInsensitive=False)

            for attribute in attributes.keys() :
                setattr(newVariable, attribute, attributes[attribute])

        self._nc.nc_enddef()

        # if data was given, use that
        if data is not None :
            newVariable[:](data.tolist())

        return newVariable

    def add_attribute_data_to_variable(self, variableName, newAttributeName, newAttributeValue) :
        """
        if the attribute exists for the given variable, set it to the new value
        if the attribute does not exist for the given variable, create it and set it to the new value
        """
        variableObject = self.get_variable_object(variableName)
        
        self._nc.nc_redef()

        setattr(variableObject, newAttributeName, newAttributeValue)

        self._nc.nc_enddef()

        # TODO, this will cause our attribute cache to be wrong!
        # TODO, for now, brute force clear the cache
        self.attributeCache = CaseInsensitiveAttributeCache(self)
        
        return
    
    def get_variable_attributes (self, variableName, caseInsensitive=True) :
        """
        returns all the attributes associated with a variable name
        """
        
        toReturn = None
        
        if caseInsensitive :
            toReturn = self.attributeCache.get_variable_attributes(variableName)
        else :
            toReturn = { }
            tempVarObj   = self.get_variable_object(variableName)
            tempAttrKeys = tempVarObj.ncattrs()
            for attrKey in tempAttrKeys :
                toReturn[attrKey] = getattr(tempVarObj, attrKey)
        
        return toReturn
    
    def get_attribute(self, variableName, attributeName, caseInsensitive=True) :
        """
        returns the value of the attribute if it is available for this variable, or None
        """
        toReturn = None
        
        if caseInsensitive :
            toReturn = self.attributeCache.get_variable_attribute(variableName, attributeName)
        else :
            temp_attributes = self.get_variable_attributes(variableName, caseInsensitive=False)
            
            if attributeName in temp_attributes :
                toReturn = getattr(self.get_variable_object, attributeName)
        
        return toReturn
    
    def get_global_attributes(self, caseInsensitive=True) :
        """
        get a list of all the global attributes for this file or None
        """
        
        toReturn = None
        
        if caseInsensitive :
            self.attributeCache.get_global_attributes()
        else :
            toReturn = { }
            tempAttrKeys = self._nc.ncattrs()
            for attrKey in tempAttrKeys :
                toReturn[attrKey] = getattr(self._nc, attrKey)
        
        return toReturn
    
    def get_global_attribute(self, attributeName, caseInsensitive=True) :
        """
        returns the value of a global attribute if it is available or None
        """
        
        toReturn = None
        
        if caseInsensitive :
            toReturn = self.attributeCache.get_global_attribute(attributeName)
        else :
            if attributeName in self._nc.ncattrs() :
                toReturn = getattr(self._nc, attributeName)
        
        return toReturn
    
    def is_loadable_type (self, name) :
        """
        check to see if the indicated variable is a type that can be loaded
        """

        return True

nc4 = nc
cdf = nc

# TODO remove
#FIXME_IDPS = [ '/All_Data/CrIS-SDR_All/ES' + ri + band for ri in ['Real','Imaginary'] for band in ['LW','MW','SW'] ] 

class h5(object):
    """wrapper for HDF5 datasets
    """
    _h5 = None
    
    def __init__(self, filename, allowWrite=False):
        self.attributeCache = CaseInsensitiveAttributeCache(self)
        
        mode = 'r'
        if allowWrite :
            mode = 'r+'
        if h5py is None:
            LOG.error('h5py module is not installed and is needed in order to read h5 files')
            assert(h5py is not None)
        self._h5 = h5py.File(filename, mode)
    
    def __call__(self):
        
        variableList = [ ]
        def testFn (name, obj) :
            #print ('checking name: ' + name)
            #print ('object: ' + str(obj))
            
            if isinstance(obj, h5py.Dataset) :
                try :
                    tempType = obj.dtype # this is required to provoke a type error for closed data sets
                    
                    #LOG.debug ('type: ' + str(tempType))
                    variableList.append(name)
                except TypeError :
                    LOG.debug('TypeError prevents the use of variable ' + name
                              + '. This variable will be ignored')
        
        self._h5.visititems(testFn)
        
        LOG.debug('variables from visiting h5 file structure: ' + str(variableList))
        
        return(variableList)
    
    @staticmethod
    def trav(h5,pth): 
        return reduce( lambda x,a: x[a] if a else x, pth.split('/'), h5)
        
    # this returns a numpy array with a copy of the full, scaled
    # data for this variable, if the data type must be changed to allow
    # for scaling it will be (so the return type may not reflect the
    # type found in the original file)
    def __getitem__(self, name):
        
        # defaults
        scale_factor = 1.0
        add_offset = 0.0
        data_type = np.float32 # TODO temporary
        
        # get the variable object and use it to
        # get our raw data and scaling info
        variable_object = self.get_variable_object(name)
        raw_data_copy = variable_object[:]
        
        #print ('*************************')
        #print (dir (variable_object.id)) # TODO, is there a way to get the scale and offset through this?
        #print ('*************************')
        
        # load the scale factor and add offset
        temp = self.attributeCache.get_variable_attributes(name)
        if (SCALE_FACTOR_STR in temp.keys()) :
            scale_factor = temp[SCALE_FACTOR_STR]
        if (ADD_OFFSET_STR in temp.keys()) :
            add_offset = temp[ADD_OFFSET_STR]
        # todo, does cdf have an equivalent of endaccess to close the variable?
        
        # don't do lots of work if we don't need to scale things
        if (scale_factor == 1.0) and (add_offset == 0.0) :
            return raw_data_copy
        
        # get information about where the data is the missing value
        missing_val = self.missing_value(name)
        missing_mask = np.zeros(raw_data_copy.shape, dtype=np.bool)
        missing_mask[raw_data_copy == missing_val] = True
        
        # create the scaled version of the data
        scaled_data_copy = np.array(raw_data_copy, dtype=data_type)
        scaled_data_copy[~missing_mask] = (scaled_data_copy[~missing_mask] * scale_factor) + add_offset #TODO, type truncation issues?
        
        return scaled_data_copy
    
    def get_variable_object(self,name):
        return h5.trav(self._h5, name)
    
    def missing_value(self, name):
        
        toReturn = None
        
        # get the missing value if it has been set
        variableObject = self.get_variable_object(name)
        pListObj = variableObject.id.get_create_plist()
        fillValueStatus = pListObj.fill_value_defined()
        if (h5d.FILL_VALUE_DEFAULT is fillValueStatus) or (h5d.FILL_VALUE_USER_DEFINED is fillValueStatus) :
            temp = np.array((1), dtype=variableObject.dtype)
            pListObj.get_fill_value(temp)
            toReturn = temp
        
        return toReturn
    
    def create_new_variable(self, variablename, missingvalue=None, data=None, variabletocopyattributesfrom=None):
        """
        create a new variable with the given name
        optionally set the missing value (fill value) and data to those given
        
        the created variable will be returned, or None if a variable could not
        be created
        """
        
        raise IOUnimplimentedError('Unable to create variable in hdf 5 file, this functionality is not yet available.')
        
        return None
    
    def add_attribute_data_to_variable(self, variableName, newAttributeName, newAttributeValue) :
        """
        if the attribute exists for the given variable, set it to the new value
        if the attribute does not exist for the given variable, create it and set it to the new value
        """
        
        raise IOUnimplimentedError('Unable to add attribute to hdf 5 file, this functionality is not yet available.')
        
        return
    
    def get_variable_attributes (self, variableName, caseInsensitive=True) :
        """
        returns all the attributes associated with a variable name
        """
        
        toReturn = None
        
        if caseInsensitive :
            toReturn = self.attributeCache.get_variable_attributes(variableName)
        else :
            toReturn = self.get_variable_object(variableName).attrs
        
        return toReturn
    
    def get_attribute(self, variableName, attributeName, caseInsensitive=True) :
        """
        returns the value of the attribute if it is available for this variable, or None
        """
        toReturn = None
        
        if caseInsensitive :
            toReturn = self.attributeCache.get_variable_attribute(variableName, attributeName)
        else :
            temp_attrs = self.get_variable_attributes(variableName, caseInsensitive=False)
            
            if (attributeName in temp_attrs) :
                toReturn = temp_attrs[attributeName]
        
        return toReturn
    
    def get_global_attributes(self, caseInsensitive=True) :
        """
        get a list of all the global attributes for this file or None
        """
        
        toReturn = None
        
        if caseInsensitive :
            self.attributeCache.get_global_attributes()
        else :
            toReturn = self._h5.attrs
        
        return toReturn
    
    def get_global_attribute(self, attributeName, caseInsensitive=True) :
        """
        returns the value of a global attribute if it is available or None
        """
        
        toReturn = None
        
        if caseInsensitive :
            toReturn = self.attributeCache.get_global_attribute(attributeName)
        else :
            if attributeName in self._h5.attrs :
                toReturn = self._h5.attrs[attributeName]
        
        return toReturn
    
    def is_loadable_type (self, name) :
        """
        check to see if the indicated variable is a type that can be loaded
        """
        
        # TODO, are there any bad types for these files?
        return True


class aeri(object):
    """wrapper for AERI RNC/SUM/CXS/etc datasets
    """
    _dmv = None
    _vectors = { }
    _scalars = { }
    
    @staticmethod
    def _meta_mapping(fp):
        ids = fp.metaIDs()
        names = [fp.queryMetaDescString(1, id_, fp.SHORTNAME) for id_ in ids]
        assert len(ids) == len(names)
        return (dict((n, i) for n, i in zip(names, ids)))
    
    def _inventory(self):
        fp = self._dmv
        assert(fp is not None)
        # get list of vectors and scalars
        self._vectors = dict( (fp.queryVectorDescString(n,fp.SHORTNAME), n) for n in fp.vectorIDs() )
        self._scalars = self._meta_mapping(fp)

    def __init__(self, filename, allowWrite=False):
        assert(allowWrite==False)
        if dmvlib is None:
            LOG.error('cannot open AERI files without dmv module being available')
            return
        self._dmv = dmvlib.dmv()
        rc = self._dmv.openFile(filename)
        if rc!=0:
            LOG.error("unable to open file, rc=%d" % rc)
            self._dmv = None        
        else:
            self._inventory()
    
    def __call__(self):
        return list(self._vectors.keys()) + list(self._scalars.keys())
        
    def __getitem__(self, name):
        fp = self._dmv
        assert(fp is not None)
        if 'DMV_RECORDS' in os.environ:
            nrecs = int(os.environ['DMV_RECORDS'])
            LOG.warning('overriding dmv record count to %d' % nrecs)
        else:
            nrecs = self._dmv.recordCount()
        recrange = range(1, nrecs+1)
        if name in self._vectors:
            vid = self._vectors[name]
            vdata = [ fp.vectorDepValues(rec, vid) for rec in recrange ]
            return np.array(vdata)
        elif name in self._scalars:
            vdata = fp.metaValueMatrix(recrange, [self._scalars[name]])
            return np.array(vdata)
        else:
            raise LookupError('cannot find variable %s' % name)
       
    def get_variable_object(self,name):
        return None
    
    def missing_value(self, name):
        return float('nan')
    
    def create_new_variable(self, variablename, missingvalue=None, data=None, variabletocopyattributesfrom=None):
        """
        create a new variable with the given name
        optionally set the missing value (fill value) and data to those given
        
        the created variable will be returned, or None if a variable could not
        be created
        """
        
        raise IOUnimplimentedError('Unable to create variable in aeri file, this functionality is not yet available.')
        
        return None
    
    def add_attribute_data_to_variable(self, variableName, newAttributeName, newAttributeValue) :
        """
        if the attribute exists for the given variable, set it to the new value
        if the attribute does not exist for the given variable, create it and set it to the new value
        """
        
        raise IOUnimplimentedError('Unable to add attribute to aeri file, this functionality is not yet available.')
        
        return
    
    def get_variable_attributes (self, variableName, caseInsensitive=True) :
        """
        returns all the attributes associated with a variable name
        """
        toReturn = { }
        
        # TODO
        LOG.warn('Glance does not yet support attribute retrieval in AERI files. None will be used.')
        
        return toReturn
    
    def get_attribute(self, variableName, attributeName, caseInsensitive=True) :
        """
        returns the value of the attribute if it is available for this variable, or None
        """
        toReturn = None
        
        # TODO
        LOG.warn('Glance does not yet support attribute retrieval in AERI files. None will be used.')
        
        return toReturn
    
    def get_global_attributes(self, caseInsensitive=True) :
        """
        get a list of all the global attributes for this file or None
        """
        
        toReturn = None
        
        # TODO
        LOG.warn('Glance does not yet support attribute retrieval in AERI files. None will be used.')
        
        return toReturn
    
    def get_global_attribute(self, attributeName, caseInsensitive=True) :
        """
        returns the value of a global attribute if it is available or None
        """
        
        toReturn = None
        
        # TODO
        LOG.warn('Glance does not yet support attribute retrieval in AERI files. None will be used.')
        
        return toReturn
    
    def is_loadable_type (self, name) :
        """
        check to see if the indicated variable is a type that can be loaded
        """
        
        # TODO, are there any bad types for these files?
        return True

# handle the variety of file suffixes by building aliases to aeri class
cxs = rnc = cxv = csv = spc = sum = uvs = aeri

class tiff (object):
    """wrapper for to open GeoTIFF data sets for comparison
    __call__ yields sequence of variable names
    __getitem__ returns individual variables ready for slicing to numpy arrays
    """
    
    _tiff = None
    
    GRAY_NAME  = "grayscale value"
    RED_NAME   = "red"
    GREEN_NAME = "green"
    BLUE_NAME  = "blue"
    IR_NAME    = "infrared"
    ALPHA_NAME = "alpha"
    
    
    # if we are using meaningful names, we will translate between
    # the band index numbers and these names (otherwise bands use generic names)
    EXPECTED_BAND_NAME_KEY = {
                                1: [GRAY_NAME],
                                2: [GRAY_NAME, ALPHA_NAME],
                                3: [RED_NAME, GREEN_NAME, BLUE_NAME],
                                4: [RED_NAME, GREEN_NAME, BLUE_NAME, ALPHA_NAME],
                                5: [RED_NAME, GREEN_NAME, BLUE_NAME, IR_NAME, ALPHA_NAME],
                             }
    
    # a reverse look up to help disambigurate what meaningful name goes with
    # which number (one of these dictionaries will be selected based on the
    # number of bands in the geotiff)
    REV_INFO               = {
                                1: {
                                    GRAY_NAME:  1,
                                    },
                                2: {
                                    GRAY_NAME:  1,
                                    ALPHA_NAME: 2,
                                    },
                                3: {
                                    RED_NAME:   1,
                                    GREEN_NAME: 2,
                                    BLUE_NAME:  3,
                                    },
                                4: {
                                    RED_NAME:   1,
                                    GREEN_NAME: 2,
                                    BLUE_NAME:  3,
                                    ALPHA_NAME: 4,
                                    },
                                5: {
                                    RED_NAME:   1,
                                    GREEN_NAME: 2,
                                    BLUE_NAME:  3,
                                    IR_NAME:    4,
                                    ALPHA_NAME: 5,
                                    },
                                
                              }
    
    def _get_generic_band_name (self, number) :
        """get a generic band name for this number"""
        
        return ("band at index " + str(number))
    
    def _get_band_index_from_name (self, name) :
        """get an index for the band from a name
        
        name may be either a meaningful name from the list that shows
        up in the reverse index keys or a generic name that was
        generated by _get_generic_band_name
        """
        
        to_return = None
        
        if name in self.revIndex.keys() :
            to_return = self.revIndex[name]
        else :
            to_return = int(name.split(' ')[-1])
        
        return to_return
    
    def __init__(self, filename, allowWrite=False, useMeaningfulNames=True):
        
        if gdal is None:
            LOG.error('gdal is not installed and is needed in order to read GeoTIFF files')
            assert(gdal is not None)
        
        if allowWrite:
            LOG.warn("Write access requested, but is not currently supported for GeoTIFF files. File will be opened read-only.")
        
        self._tiff     = gdal.Open(filename)
        self.niceNames = useMeaningfulNames
        self.revIndex  = self.REV_INFO[self._tiff.RasterCount] if self._tiff.RasterCount in self.REV_INFO else { }

    def __call__(self):
        "yield names of variables to be compared"
        
        # GeoTIFF files don't actually have named variables, so get something appropriate based on the numbering of bands
        num_bands = self._tiff.RasterCount
        
        to_return = [ ]
        if self.niceNames and (num_bands in self.EXPECTED_BAND_NAME_KEY.keys()) :
            to_return = self.EXPECTED_BAND_NAME_KEY[num_bands][:]
        else :
            for bandNumber in range(1, num_bands + 1) :
                to_return.append(self._get_generic_band_name(bandNumber))
        
        return to_return
    
    # this returns a numpy array with a copy of the full, scaled
    # data for this variable, if the data type must be changed to allow
    # for scaling it will be (so the return type may not reflect the
    # type found in the original file)
    def __getitem__(self, name):
        
        LOG.debug("opening variable: " + name)
        
        # first figure out the index for this variable
        var_index = self._get_band_index_from_name(name)
        
        # get the data out of the file
        temp_band = self._tiff.GetRasterBand(var_index)
        temp_data = temp_band.ReadAsArray()
        
        # there is no standard scaling procedure for GeoTIFFs, so skip that!
        
        return temp_data
    
    # TODO, this hasn't been supported in other file types
    def close (self) :
        # Dave couldn't find any explicit way to close it other
        # than let the garbage collector take care of it
        self._tiff = None
    
    def get_variable_object(self, name):
        return None
    
    def missing_value(self, name):
        return None
    
    def create_new_variable(self, variablename, missingvalue=None, data=None, variabletocopyattributesfrom=None):
        """
        TODO, this is not yet supported for GeoTIFFs
        
        create a new variable with the given name
        optionally set the missing value (fill value) and data to those given
        
        the created variable will be returned, or None if a variable could not
        be created
        """
        
        LOG.warn ("GeoTIFF io class does not yet support writing. Unable to create new variable in GeoTIFF file.")
        
        return None
    
    def add_attribute_data_to_variable(self, variableName, newAttributeName, newAttributeValue) :
        """
        TODO this is not yet supported for GeoTIFFs
        
        if the attribute exists for the given variable, set it to the new value
        if the attribute does not exist for the given variable, create it and set it to the new value
        """
        
        LOG.warn ("GeoTIFF io class does not yet support writing. Unable to add attribute information to GeoTIFF file.")
        
        return
    
    def get_variable_attributes (self, variableName, caseInsensitive=True) :
        """
        returns all the attributes associated with a variable name
        """
        
        # FUTURE, GeoTIFF files do have attributes, but this isn't hooked up yet
        
        return { }
    
    def get_attribute(self, variableName, attributeName, caseInsensitive=True) :
        """
        returns the value of the attribute if it is available for this variable, or None
        """
        
        # FUTURE, GeoTIFF files do have attributes, but this isn't hooked up yet
        
        return None
    
    def get_global_attributes(self, caseInsensitive=True) :
        """
        get a list of all the global attributes for this file or None
        """
        
        # FUTURE, GeoTIFF files do have attributes, but this isn't hooked up yet
        
        return { }
    
    def get_global_attribute(self, attributeName, caseInsensitive=True) :
        """
        returns the value of a global attribute if it is available or None
        """
        
        # FUTURE, GeoTIFF files do have attributes, but this isn't hooked up yet
        
        return None
    
    def is_loadable_type (self, name) :
        """
        check to see if the indicated variable is a type that can be loaded
        """
        
        # TODO, are there any bad types for these files?
        return True

# people also name tiff files with one f...
tif  = tiff
# Nick has special tiff files with alpha...
tifa = tiff

def _search_xml(pathname):
    xs = '.xml'
    yield pathname + xs
    yield os.path.splitext(pathname)[0] + xs
    yield pathname.replace('-', '_') + xs
    yield os.path.splitext(pathname)[0].replace('-', '_') + xs

class jpss_adl(object):
    """wrapper for JPSS ADL BLOBs 
    This is a somewhat unique case in that the BLOB loader requires both an XML path and a BLOB path.
    In this case, it is assumed that a softlinked pathname.xml exists for a given pathname.
    FORMAT=jpss_adl glance stats truth/ATMS-FSDR.BE ATMS-FSDR
    """
    _blob = None

    def __init__(self, filename, allowWrite=False):
        assert(allowWrite==False)
        for xmlname in _search_xml(filename):
            if not os.path.exists(xmlname): 
                continue
            LOG.info('using %s for %s' % (xmlname, filename))
            break
        if not os.path.exists(xmlname):
            LOG.error(xmlname + ' needs to provide layout for ' + filename)
            return            
        if adl_blob is None:
            LOG.error('cannot open JPSS ADL files without adl_blob module in $PYTHONPATH')
            return
        if filename.lower().endswith('.be'):
            endian = adl_blob.BIG_ENDIAN
        elif filename.lower().endswith('.be'):
            endian = adl_blob.LITTLE_ENDIAN
        else:
            endian = adl_blob.NATIVE_ENDIAN
        LOG.debug('endianness of %s is %s' % (filename, endian))
        self._blob = adl_blob.map(xmlname, filename, writable=False, endian=endian)        
    
    def __call__(self):
        fieldnames = [name for name,field in self._blob._fields_]
        return fieldnames
        
    def __getitem__(self, name):
        field = getattr(self._blob, name)
        if not hasattr(field,'_length_'): # FUTURE: is this rigorous? 
            LOG.info('creating numpy array out of singleton value for %s' % name)
            return np.array([field])
        return np.array(field)
       
    def get_variable_object(self,name):
        return None
    
    def missing_value(self, name):
        return float('nan')
    
    def create_new_variable(self, variablename, missingvalue=None, data=None, variabletocopyattributesfrom=None):
        """
        create a new variable with the given name
        optionally set the missing value (fill value) and data to those given
        
        the created variable will be returned, or None if a variable could not
        be created
        """
        
        raise IOUnimplimentedError('Unable to create variable in JPSS ADL file, this functionality is not yet available.')
        
        return None
    
    def add_attribute_data_to_variable(self, variableName, newAttributeName, newAttributeValue) :
        """
        if the attribute exists for the given variable, set it to the new value
        if the attribute does not exist for the given variable, create it and set it to the new value
        """
        
        raise IOUnimplimentedError('Unable to add attribute to JPSS ADL file, this functionality is not yet available.')
        
        return
    
    def get_variable_attributes (self, variableName, caseInsensitive=True) :
        """
        returns all the attributes associated with a variable name
        """
        toReturn = { }
        
        # TODO
        LOG.warn('Glance does not yet support attribute retrieval in JPSS ADL files. None will be used.')
        
        return toReturn
    
    def get_attribute(self, variableName, attributeName, caseInsensitive=True) :
        """
        returns the value of the attribute if it is available for this variable, or None
        """
        toReturn = None
        
        # TODO
        LOG.warn('Glance does not yet support attribute retrieval in JPSS ADL files. None will be used.')
        
        return toReturn
    
    def get_global_attribute(self, attributeName, caseInsensitive=True) :
        """
        returns the value of a global attribute if it is available or None
        """
        
        toReturn = None
        
        # TODO
        LOG.warn('Glance does not yet support attribute retrieval in JPSS ADL files. None will be used.')
        
        return toReturn
    
    def get_global_attributes(self, caseInsensitive=True) :
        """
        get a list of all the global attributes for this file or None
        """
        
        toReturn = None
        
        # TODO
        LOG.warn('Glance does not yet support attribute retrieval in JPSS ADL files. None will be used.')
        
        return toReturn
    
    def is_loadable_type (self, name) :
        """
        check to see if the indicated variable is a type that can be loaded
        """
        
        # TODO, are there any bad types for these files?
        return True

def open(pathname, allowWrite=False):
    suffix = os.path.splitext(pathname)[1][1:].lower()

    # Just test we can open the file so we automatically raise a suitable
    # error if we can't access it.
    from __builtin__ import open
    with open(pathname):
        pass

    if (not suffix) or (suffix not in globals()):
        # this ican be used to specify a format on the command line by setting the
        # environment variable FORMAT, for example:
        #           export FORMAT=nc
        suffix = os.environ.get('FORMAT', None)
        LOG.info('overriding unknown load format to "%s"' % suffix)
    cls = globals()[suffix]
    return cls(pathname, allowWrite=allowWrite)

if __name__=='__main__':
    import doctest
    doctest.testmod()
