#!/usr/bin/env python
# encoding: utf-8
"""
I/O routines supporting reading a number of file formats.

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging

from pyhdf.SD import SD,SDC, SDS, HDF4Error
try:
    import h5py
except ImportError:
    pass
from pycdf import CDF, NC

import numpy as np

LOG = logging.getLogger(__name__)

class hdf(SD):
    """wrapper for HDF4 dataset for comparison
    __call__ yields sequence of variable names
    __getitem__ returns individual variables ready for slicing to numpy arrays
    """
    
    def __init__(self,filename):
        super(self.__class__,self).__init__(filename, SDC.READ)

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
        
        # get the variable object and use it to
        # get our raw data and scaling info
        variable_object = self.get_variable_object(name)
        raw_data_copy = variable_object[:]
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
        missing_value_attr_name = '_FillValue'
        variable_object = self.select(name)
        
        to_return = None
        if hasattr(variable_object, missing_value_attr_name) :
            to_return = getattr(variable_object, missing_value_attr_name, None)
        SDS.endaccess(variable_object)
        
        return to_return
        

class nc(CDF):
    """wrapper for NetCDF3/4/opendap dataset for comparison
    __call__ yields sequence of variable names
    __getitem__ returns individual variables ready for slicing to numpy arrays
    """
    
    def __init__(self,filename):
        super(self.__class__,self).__init__(filename, NC.NOWRITE)

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
        
        missing_value_attr_name_1 = '_FillValue'
        missing_value_attr_name_2 = 'missing_value'
        variable_object = self.var(name)
        
        to_return = None
        if hasattr(variable_object, missing_value_attr_name_1) \
           or \
           hasattr(variable_object, missing_value_attr_name_2) :
            to_return = getattr(variable_object, missing_value_attr_name_1,
                                getattr(variable_object, missing_value_attr_name_2, None))
        
        return to_return

nc4 = nc
cdf = nc

# TODO remove
#FIXME_IDPS = [ '/All_Data/CrIS-SDR_All/ES' + ri + band for ri in ['Real','Imaginary'] for band in ['LW','MW','SW'] ] 

class h5(object):
    """wrapper for HDF5 datasets
    """
    _h5 = None
    
    def __init__(self,filename):
        self._h5 = h5py.File(filename,'r')
    
    def __call__(self):
        
        variableList = [ ]
        def testFn (name, obj) :
            #print ('checking name: ' + name)
            print ('object: ' + str(obj))
            
            if isinstance(obj, h5py.Dataset) :
                try :
                    tempType = obj.dtype # this is required to provoke a type error for closed data sets
                    #TODO, why are there closed data sets?!?
                    LOG.debug ('type: ' + str(tempType))
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
        return None
        

def open(pathname):
    cls = globals()[os.path.splitext(pathname)[1][1:]]
    return cls(pathname)



if __name__=='__main__':
    import doctest
    doctest.testmod()
