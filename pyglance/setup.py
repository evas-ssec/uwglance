#!/usr/bin/env python
"""
$Id: setup.py 66 2009-04-24 18:47:55Z rayg $
see http://peak.telecommunity.com/DevCenter/setuptools

note: if you want to develop this code and run from code on the command line,
please run the following line when you update to a new version of the code.
python setup.py develop --install-dir=$HOME/Library/Python

distribution:
python setup.py develop --install-dir=$HOME/Library/Python
python setup.py sdist
python setup.py bdist_egg
(cd dist; rsync -Cuav * larch.ssec.wisc.edu:/var/apache/larch/htdocs/eggs/repos/glance)

use: 
python setup.py install --install-dir=$HOME/Library/Python
easy_install -d $HOME/Library/Python -vi http://larch.ssec.wisc.edu/eggs/repos glance
"""

# changed to support egg distribution
from setuptools import setup, find_packages

setup( name="glance", 
       version="0.2.6.19", 
       zip_safe = True,
       entry_points = { 'console_scripts': [ 'glance = glance.compare:main' ] },
       packages = find_packages('.'),
       install_requires=[ 'numpy' ],
       package_data = {'': ['*.txt', '*.gif'], }
       )

