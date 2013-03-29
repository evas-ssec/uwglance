#!/usr/bin/env python
# encoding: utf-8
"""
An auxillary module for the glance GUI that holds information on constants. 

Created by evas June 2012.
Copyright (c) 2012 University of Wisconsin SSEC. All rights reserved.
"""

import logging
#import numpy as np

LOG = logging.getLogger(__name__)

"""
The constants module handles generally useful constants in the Glance GUI.
It not only stores values that are intended to be unchangable and globally available.
"""

# constants for the possible image types
ORIGINAL_A = "Original A Data"
ORIGINAL_B = "Original B Data"
ABS_DIFF   = "Abs. Difference"
RAW_DIFF   = "Raw Difference"
HISTOGRAM  = "Histogram"
MISMATCH   = "Mismatch Areas"
SCATTER    = "Scatter Plot"
HEX_PLOT   = "Hex Plot"

# a list of all the image types, for convenience
IMAGE_TYPES = [ORIGINAL_A,
               ORIGINAL_B,
               ABS_DIFF,
               RAW_DIFF,
               HISTOGRAM,
               MISMATCH,
               SCATTER,
               HEX_PLOT
              ]

# a list of image types that require both the A and B data
COMPARISON_IMAGES = [ABS_DIFF,
                     RAW_DIFF,
                     HISTOGRAM,
                     MISMATCH,
                     SCATTER,
                     HEX_PLOT
                    ]

# constants for possible types of data handling
SIMPLE_2D = "Simple Two Dimensional"
MAPPED_2D = "Mapped Two Dimensional"
ONLY_1D   = "One Dimensional"

# a list of data handling types, for conveinence
DATA_FORMS = [SIMPLE_2D,
#              MAPPED_2D, TODO add this back in when the mapped stuff is ready
              ONLY_1D]

# the default names that the model will try to select for the latitude and longitude
DEFAULT_LONGITUDE = 'pixel_longitude'
DEFAULT_LATITUDE  = 'pixel_latitude'

# the number of bins to use for histograms
DEFAULT_NUM_BINS = 50

# colormaps that are available in the GUI
# TODO, this needs to be upkept when the list of colormaps in the figure manager changes
CM_RAINBOW     = "Rainbow"
CM_GRAY        = "Grayscale"
COLORMAP_NAMES = [CM_RAINBOW, CM_GRAY]

NO_DATA_MESSAGE = "Requested data was not available or did not exist."
UNKNOWN_DATA_FORM = "An ununsupported plot format was requested. Aborting attempt to plot data."


if __name__=='__main__':
    import doctest
    doctest.testmod()