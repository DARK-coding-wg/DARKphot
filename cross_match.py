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
import darkphot_local_catalog as dp_local

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


def cross_match(desired_pos, full_table, name_ra_col = 'RA', name_dec_col='Dec', match_radius=1.5):
    """
    Keywords:
        match_radius - matching radius in arcsec
    """
    ra_input, dec_input = zip(*desired_pos)
    coord_input = SkyCoord(ra=ra_input*u.degree, dec=dec_input*u.degree)

    ra_catalog = np.array(full_table[name_ra_col])
    dec_catalog = np.array(full_table[name_dec_col])
    coord_catalog = SkyCoord(ra=ra_catalog*u.degree, dec=dec_catalog*u.degree)

    idx, d2d, d3d = coord_input.match_to_catalog_sky(coord_catalog)

    return idx



#####################
###### Run Code #####
#####################
if __name__ == '__main__':

    local_cat = dp_local.LocalCatalog()
    ## A couple sources, including one that doesn't have any counterparts
    fname_catalog = 'example_data/test_cat.csv'
    fname_param = 'example_data/input_local_cat.param'
    fname_object_names_sel = 'example_data/input_names.txt'
    local_cat.add_local_catalog(fname_catalog, fname_param, fname_object_names_sel, name_id_col='NAME')
    full_catalog = local_cat._read_text_catalog('csv')
    target_pos = local_cat.object_sel_coordinates



