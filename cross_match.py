#! /usr/bin/env python
#
# Program: query_vizier
#
# Author: DARKphot
#
# Usage: ./query_vizier
#
# Description: Code to query the photometric table of Vizier
# Import Libraries
# import numpy as np
# import pandas as pd
# import os
# import warnings
# import requests
# from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

# Set up loggers
import logging
# Create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# Create console handler
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.WARNING)
# Create text log handler
# file_handler = logging.FileHandler('logs/cross_match.log', mode='w')
# file_handler.setLevel(logging.DEBUG)
# Create formatter for both handlers
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(funcName)s - %(message)s')
console_handler.setFormatter(formatter)
# file_handler.setFormatter(formatter)
# Add both handlers to logger
logger.addHandler(console_handler)
# logger.addHandler(file_handler)

def cross_match(pos, full_table, name_ra_col = 'RA', name_dec_col='Dec'):
    """ """
    ra_input, dec_input = zip(*pos)
    coord_input = SkyCoord(ra=ra_input*u.degree, dec=dec_input*u.degree)

    ra_catalog = np.array(full_table[name_ra_col])
    dec_catalog = np.array(full_table[name_dec_col])e
    coord_catalog = SkyCoord(ra=ra_catalog*u.degree, dec=dec_catalog*u.degree)

    idx, d2d, d3d = coord_input.match_to_catalog_sky(coord_catalog)

    return idx

