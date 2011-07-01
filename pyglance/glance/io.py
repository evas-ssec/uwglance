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

try:    
    import pycdf
    from pycdf import CDF, NC, strerror
except:
    LOG.info('no pycdf module available')
    pycdf = None
    CDF = NC = object
    def strerror(*args):
        return 'no pycdf module installed'

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


fillValConst1 = '_FillValue'
fillValConst2 = 'missing_value'

class IOUnimplimentedError(Exception):
    """
    The exception raised when a requested io operation is not yet available.
    
        msg  -- explanation of the problem
    """
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

class hdf(SD):
    """wrapper for HDF4 dataset for comparison
    __call__ yields sequence of variable names
    __getitem__ returns individual variables ready for slicing to numpy arrays
    """
    
    def __init__(self, filename, allowWrite=False):
        if pyhdf is None:
            LOG.error('pyhdf is not installed and is needed in order to read hdf4 files')
            assert(pyhdf is not None)
        mode = SDC.READ
        if allowWrite:
            mode = mode | SDC.WRITE
        super(self.__class__,self).__init__(filename, mode)

    def __call__(self):
        "yield names of variables to be compared"
        return self.datasets().keys()
    
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
        
        print ("***** getting " + name + " from file")
        
        # get the variable object and use it to
        # get our raw data and scaling info
        variable_object = self.get_variable_object(name)
        raw_data_copy = variable_object[:]
        print ("****** raw data loaded")
        try :
            # TODO, this currently won't work with geocat data, work around it for now
            scale_factor, scale_factor_error, add_offset, add_offset_error, data_type = SDS.getcal(variable_object)
        except HDF4Error:
            # load just the scale factor and add offset
            temp_attributes = variable_object.attributes()
            if ('add_offset' in temp_attributes) :
                add_offset = temp_attributes['add_offset']
                data_type = np.dtype(type(add_offset))
            if ('scale_factor' in temp_attributes) :
                scale_factor = temp_attributes['scale_factor']
                data_type = np.dtype(type(scale_factor))
            if ('scaling_method' in temp_attributes) :
                scaling_method = temp_attributes['scaling_method']
        SDS.endaccess(variable_object)
        
        print ("***** scaling information loaded")
        
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
        scaled_data_copy = np.array(raw_data_copy, dtype=data_type)
        #scaled_data_copy[~missing_mask] = (scaled_data_copy[~missing_mask] - add_offset) * scale_factor #TODO, type truncation issues?
        scaled_data_copy[~missing_mask] = (scaled_data_copy[~missing_mask] * scale_factor) + add_offset #TODO, type truncation issues?
        
        return scaled_data_copy 
    
    def get_variable_object(self, name):
        return self.select(name)
    
    def missing_value(self, name):
        variable_object = self.select(name)
        
        to_return = None
        if hasattr(variable_object, fillValConst1) :
            to_return = getattr(variable_object, fillValConst1, None)
        SDS.endaccess(variable_object)
        
        return to_return
    
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
    
    def get_attribute(self, variableName, attributeName) :
        """
        returns the value of the attribute if it is available for this variable, or None
        """
        toReturn = None
        
        variable_object = self.get_variable_object(variableName)
        temp_attributes = variable_object.attributes()
        
        if attributeName in temp_attributes :
            toReturn = temp_attributes[attributeName]
        
        return toReturn

class nc(CDF):
    """wrapper for NetCDF3/4/opendap dataset for comparison
    __call__ yields sequence of variable names
    __getitem__ returns individual variables ready for slicing to numpy arrays
    """
    
    def __init__(self, filename, allowWrite=False):
        
        if pycdf is None:
            LOG.error('pycdf is not installed and is needed in order to read NetCDF files')
            assert(pycdf is not None)

        
        mode = NC.NOWRITE
        if allowWrite :
            mode = NC.WRITE
        
        super(self.__class__,self).__init__(filename, mode)

    def __call__(self):
        "yield names of variables to be compared"
        return self.variables().keys()
    
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
        # load the scale factor and add offset
        temp_attributes = variable_object.attributes()
        if ('scale_factor' in temp_attributes) :
            scale_factor = temp_attributes['scale_factor']
        if ('add_offset' in temp_attributes) :
            add_offset = temp_attributes['add_offset']
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
        #scaled_data_copy[~missing_mask] = (scaled_data_copy[~missing_mask] - add_offset) * scale_factor #TODO, type truncation issues?
        scaled_data_copy[~missing_mask] = (scaled_data_copy[~missing_mask] * scale_factor) + add_offset #TODO, type truncation issues?
        
        return scaled_data_copy 
    
    def get_variable_object(self, name):
        return self.var(name)
    
    def missing_value(self, name):
        
        variable_object = self.var(name)
        
        to_return = None
        if hasattr(variable_object, fillValConst1) \
           or \
           hasattr(variable_object, fillValConst2) :
            to_return = getattr(variable_object, fillValConst1,
                                getattr(variable_object, fillValConst2, None))
        
        return to_return
    
    def create_new_variable(self, variablename, missingvalue=None, data=None, variabletocopyattributesfrom=None):
        """
        create a new variable with the given name
        optionally set the missing value (fill value) and data to those given
        
        the created variable will be returned, or None if a variable could not
        be created
        """
        
        self.redef()
        
        # if the variable already exists, stop with a warning
        if variablename in self.variables().keys() :
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
            dataType = NC.INT
            #print("Picked INT")
        # TODO, at the moment the fill type is forcing me to use a double, when sometimes I want a float
        #elif np.issubdtype(data.dtype, np.float32) :
        #    dataType = NC.FLOAT
        #    print("Picked FLOAT")
        elif np.issubdtype(data.dtype, float) :
            dataType = NC.DOUBLE
            #print("Picked DOUBLE")
        # what do we do if it's some other type?
        
        # create and set all the dimensions
        dimensions = [ ]
        dimensionNum = 0
        for dimSize in data.shape :
            dimensions.append(self.def_dim(variablename + '-index' + str(dimensionNum), dimSize))
            dimensionNum = dimensionNum + 1
        
        # create the new variable
        #print('variable name: ' + variablename)
        #print('data type:     ' + str(dataType))
        #print('dimensions:    ' + str(dimensions))
        newVariable = self.def_var(variablename, dataType, tuple(dimensions))
        
        # if a missing value was given, use that
        if missingvalue is not None :
            newVariable._FillValue = missingvalue
        
        # if we have a variable to copy attributes from, do so
        if variabletocopyattributesfrom is not None :
            tocopyfrom = self.get_variable_object(variabletocopyattributesfrom)
            attributes = tocopyfrom.attributes()
            for attribute in attributes.keys() :
                newVariable.__setattr__(attribute, attributes[attribute])
        
        self.enddef()
        
        # if data was given, use that
        if data is not None :
            newVariable.put(data.tolist()) 
        
        return newVariable
    
    def add_attribute_data_to_variable(self, variableName, newAttributeName, newAttributeValue) :
        """
        if the attribute exists for the given variable, set it to the new value
        if the attribute does not exist for the given variable, create it and set it to the new value
        """
        variableObject = self.get_variable_object(variableName)
        
        self.redef()
        
        variableObject.__setattr__(newAttributeName, newAttributeValue)
        
        self.enddef()
        
        return
    
    def get_attribute(self, variableName, attributeName) :
        """
        returns the value of the attribute if it is available for this variable, or None
        """
        toReturn = None
        
        variable_object = self.get_variable_object(variableName)
        temp_attributes = variable_object.attributes()
        
        if attributeName in temp_attributes :
            toReturn = temp_attributes[attributeName]
        
        return toReturn


nc4 = nc
cdf = nc

# TODO remove
#FIXME_IDPS = [ '/All_Data/CrIS-SDR_All/ES' + ri + band for ri in ['Real','Imaginary'] for band in ['LW','MW','SW'] ] 

class h5(object):
    """wrapper for HDF5 datasets
    """
    _h5 = None
    
    def __init__(self, filename, allowWrite=False):
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
        if ('scale_factor' in variable_object.attrs) :
            scale_factor = variable_object.attrs['scale_factor']
        if ('add_offset' in variable_object.attrs) :
            add_offset = variable_object.attrs['add_offset']
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
        #scaled_data_copy[~missing_mask] = (scaled_data_copy[~missing_mask] - add_offset) * scale_factor #TODO, type truncation issues?
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
    
    def get_attribute(self, variableName, attributeName) :
        """
        returns the value of the attribute if it is available for this variable, or None
        """
        toReturn = None
        
        variable_object = self.get_variable_object(variableName)
        if (attributeName in variable_object.attrs) :
            toReturn = variable_object.attrs[attributeName]
        
        return toReturn




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
    
    def get_attribute(self, variableName, attributeName) :
        """
        returns the value of the attribute if it is available for this variable, or None
        """
        toReturn = None
        
        # TODO
        LOG.warn('Glance does not yet support attribute retrieval in AERI files. None will be used.')
        
        return toReturn

# handle the variety of file suffixes by building aliases to aeri class
cxs = rnc = cxv = csv = spc = sum = uvs = aeri


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
    
    def get_attribute(self, variableName, attributeName) :
        """
        returns the value of the attribute if it is available for this variable, or None
        """
        toReturn = None
        
        # TODO
        LOG.warn('Glance does not yet support attribute retrieval in JPSS ADL files. None will be used.')
        
        return toReturn




def open(pathname, allowWrite=False):
    suffix = os.path.splitext(pathname)[1][1:].lower()
    if (not suffix) or (suffix not in globals()):
        suffix = os.environ.get('FORMAT', None)
        LOG.info('overriding unknown load format to "%s"' % suffix)
    cls = globals()[suffix]
    return cls(pathname, allowWrite=allowWrite)



if __name__=='__main__':
    import doctest
    doctest.testmod()
