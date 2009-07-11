#!/usr/bin/env python
# encoding: utf-8
"""
I/O routines supporting reading a number of file formats.

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging

from pyhdf.SD import SD,SDC
try:
    import h5py
except ImportError:
    pass
from pycdf import CDF, NC

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
    
    def __getitem__(self, name):
        return self.select(name)
        
    def missing_value(self, name):
        return getattr(self.select(name),'_FillValue',None)
        

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
    
    def __getitem__(self, name):
        return self.var(name)
        
    def missing_value(self, name):
        return getattr(self.var(name),'_FillValue',getattr(self.var(name),'missing_value',None))
nc4 = nc
cdf = nc


FIXME_IDPS = [ '/All_Data/CrIS-SDR_All/ES' + ri + band for ri in ['Real','Imaginary'] for band in ['LW','MW','SW'] ] 

class h5(object):
    """wrapper for HDF5 datasets
    """
    _h5 = None
    
    def __init__(self,filename):
        self._h5 = h5py.File(filename,'r')
    
    def __call__(self):
        "FIXME: this should return the real list of variables, which will include slashes"
        return set(FIXME_IDPS)
    
    @staticmethod
    def trav(h5,pth): 
        return reduce( lambda x,a: x[a] if a else x, pth.split('/'), h5)
        
    def __getitem__(self,name):
        return h5.trav(self._h5, name)
    
    def missing_value(self, name):
        return None
        

def open(pathname):
    cls = globals()[os.path.splitext(pathname)[1][1:]]
    return cls(pathname)



if __name__=='__main__':
    import doctest
    doctest.testmod()
