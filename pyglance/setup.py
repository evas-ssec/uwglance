#!/usr/bin/env python
"""
$Id: setup.py 66 2009-04-24 18:47:55Z rayg $
see http://peak.telecommunity.com/DevCenter/setuptools

distribution:
python setup.py develop --install-dir=$HOME/Library/Python
python setup.py sdist
python setup.py bdist_egg
(cd dist; rsync -Cuav * larch.ssec.wisc.edu:/home/httpd/html/eggs/repos/glance/)

use: 
python setup.py install --install-dir=$HOME/Library/Python
easy_install -d $HOME/Library/Python -vi http://larch.ssec.wisc.edu/eggs/repos glance
"""

# changed to support egg distribution
from setuptools import setup, find_packages

setup( name="glance", 
       version="0.2.6.1", 
       zip_safe = True,
       entry_points = { 'console_scripts': [ 'glance = glance.compare:main' ] },
       packages = find_packages('.'),
       install_requires=[ 'numpy' ],
       package_data = {'': ['*.txt'], }
       )

