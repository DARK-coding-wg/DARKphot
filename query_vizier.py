#! /usr/bin/env python
# 
# Program: query_vizier
#
# Author: Nick Lee
#
# Usage: ./query_vizier
#
# Description: Code to query the photometric table of Vizier
#
# To Do:
#    Tests
#    Additional functionality for different formats of RA/Dec
#    Querying by name
#

# Import Libraries
import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import math
import pdb
import os
import sys
import argparse
import requests
from astropy.table import Table

## Debugging module and intital setup
import logging
logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(levelname)s- %(message)s')

###############
### Classes ###
###############


#################
### Functions ###
#################
def get_all_dataframes(ra_list, dec_list):
    """ Sloppy way to run through a list of sources """
    list_of_frames = []
    for ra, dec in zip(ra_list, dec_list):
        url = create_url_from_pos(ra, dec, radius=1.5)
        download_from_vizier(url, 'test_vo_table.vot')
        try:
            data = read_vo_table('test_vo_table.vot')
        except ValueError:
            print 'no object at coordinates {:} {:+f}'.format(ra, dec)
            list_of_frames.append(-99)
        else:
            list_of_frames.append(data)
        finally:
            os.remove('test_vo_table.vot')
    return list_of_frames

def read_vo_table(filename):
    """
    Reads in a VOTable and returns a pandas dataframe containing the data.

    Args:
        filename - Full path filename of the VOTable. Should have a file extension of .vot
    Returns:
        pandas_table - A pandas dataframe containing the data from a VOTable
    Raises:
        ValueError: The VOTable exists but is empty. 
        IOError: Improper filename provided. Most likely the file does not exist.
    """
    try:
        tab = Table.read(filename)
    except ValueError as exc:
        logging.critical('Input VOTable at {:} is empty'.format(filename))
        raise 
    except IOError as exc:
        logging.critical('Invalid filename passed to read_vo_table')
        raise 
    return tab.to_pandas()

def download_from_vizier(url, filename):
    """ Download the VOTable from Vizier Photometric Table and save to designated filename.

    Args:
        url - Full url to Vizier Photometric Table query
        filename - local filename to save the VOTable
    Returns:
        None
    Side Effects:
        Downloads and saves the VOTable found at url to filename.
    """
    response = requests.get(url, stream=True)
    response.raise_for_status()       # Make sure url loaded properly
    with open(filename, 'wb') as vo_table:
        for ii, chunk in enumerate(response.iter_content(chunk_size=1024)):
            if chunk:                 # filter out keep-alive new chunks
                vo_table.write(chunk)
                #f.flush()

def create_url_from_pos(ra, dec, radius=1.5):
    """ Create the URL used to query Vizier's photometric Viewer based on astrometry

    Uses the format found at http://vizier.u-strasbg.fr/vizier/sed/doc/

    Args:
        ra - Right Ascension in Decimal Degrees
        dec - Declination in Decimal Degrees
    Keywords:
        radius - Search radius in arcsec. Default is 1.5 arcsec. 
    Returns:
        url - URL where VOTable from Vizier Photometric Table can be found
    """
    url = 'http://vizier.u-strasbg.fr/viz-bin/sed?-c={:f}{:+f}&-c.rs={:f}'.format(ra,dec,radius)
    return url

####################
####################
def parse_args():
    '''
    Read command line arguments
    '''
    ## Define a parser
    parser = argparse.ArgumentParser(description='Modules to query the Vizier Photometric Catalog')

    ## Add arguments
    parser.add_argument('-d','--debug',action='store_true',help='Turn on debugging messages')

    ## Use parser to interpret the command line arguments
    inputs = vars(parser.parse_args())

    ## Set debugging options
    if inputs['debug']==False:
        logging.disable(logging.INFO)  # Turns off debug messages
    else:
        logging.disable(logging.NOTSET)
        logging.debug('Start of program')

#####################
###### Run Code #####
#####################
if __name__ == '__main__':
    ## Process command line arguments
    parse_args()

    ## A couple sources, including one that doesn't have any counterparts
    ra_list = [187.27832916, 150.231314, 149.426254]
    dec_list = [2.05199, -4.1234,  2.073906]

    all_photometry = get_all_dataframes(ra_list, dec_list)
  

        
