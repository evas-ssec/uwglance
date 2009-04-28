#!/usr/bin/env python
# encoding: utf-8
"""
I/O routines supporting reading a number of file formats.

Created by rayg Apr 2009.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging

from pyhdf.SD import SD,SDC

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
        

class h5(object):
    pass

def open(pathname):
    cls = globals()[os.path.splitext(pathname)[1][1:]]
    return cls(pathname)



if __name__=='__main__':
    import doctest
    doctest.testmod()
