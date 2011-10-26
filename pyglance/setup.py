#!/usr/bin/env python
"""
$Id$
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
       version="0.2.6.29", 
       zip_safe = True,
       entry_points = { 'console_scripts': [ 'glance = glance.compare:main' ] },
       packages = find_packages('.'),
       install_requires=[ 'numpy', 'matplotlib' ],
       package_data = {'': ['*.txt', '*.gif'], }
       )

