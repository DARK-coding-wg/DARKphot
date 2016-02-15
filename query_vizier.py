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
from astropy.coordinates import SkyCoord

## Debugging module and intital setup
import logging
logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(levelname)s- %(message)s')

###############
### Classes ###
###############
class VizierCatalog(object):
    """ Class for creating a photometric catalog from Vizier """
        
    def get_all_dataframes(self, source_list, source_id=None, radius=1.5):
        """ Create a catalog containing all the photometric data from Vizier for a list of sources

        Returns full output from Vizier, with a single added column of source_id.
        Args:
            source_list - list of souces. Can be a list of (ra, dec) tuples or a list of names as strings. RA & Dec must be in decimal degrees.
        Keywords:
            radius - Defines region to search for counterparts in arcsec. Default is 1.5"
            source_id - list of source_ids to use. If not provided, will use the index of the source_list
        Returns:
            all_frames - a pandas dataframe that contains all observations from Vizier and an additional source_id column for identifying sources
        """
        if source_id is None:
            source_id = np.arange(len(source_list))
        list_of_frames = []
        for ii, source in enumerate(source_list):
            url = self._create_url(source, radius=radius)
            self._download_from_vizier(url, 'test_vo_table.vot')
            try:
                photometry = self._read_vo_table('test_vo_table.vot')
            except ValueError:
                logging.warning('No observations found via Vizier for source_id = {:}'.format(source_id[ii]))
            else:
                n_rows = photometry.shape[0]
                photometry['source_id'] = np.repeat(source_id[ii], n_rows)
                list_of_frames.append(photometry)
            finally:
                os.remove('test_vo_table.vot')
        all_frames = self._replace_frequency_with_wavelength(pd.concat(list_of_frames))
        return all_frames  

    def _read_vo_table(self, filename):
        """
        Reads in a VOTable and returns a pandas dataframe containing the data (private).
    
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
            raise 
        except IOError as exc:
            logging.critical('Invalid filename passed to read_vo_table')
            raise 
        return tab.to_pandas()
    
    def _download_from_vizier(self, url, filename):
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

    def _create_url(self, source, radius=1.5):
        """ Make the correct URL for the Vizier query

        Uses the format found at http://vizier.u-strasbg.fr/vizier/sed/doc/
        Args:
            source - Either a tuple of (ra, dec) in decimal degrees or a string that represents the name of the source
        Keywords:
            radius - Search radius in arcsec. Default is 1.5 arcsec. 
        Returns:
            url - URL where VOTable from Vizier Photometric Table can be found
        """
        if isinstance(source, tuple):
            url = 'http://vizier.u-strasbg.fr/viz-bin/sed?-c={:f}{:+f}&-c.rs={:f}'.format(source[0], source[1], radius)
        elif isinstance(source, str):
            url = 'http://vizier.u-strasbg.fr/viz-bin/sed?-c={:}&-c.rs={:f}'.format(source, radius)
        else:
            raise IOError('Provided Source is not of correct type. Must be either a tuple of (ra, dec) or a string of name')
        return url
        
    def _replace_frequency_with_wavelength(self, catalog):
        """ Converts the column sed_freq to sed_wave with the appropriate units """
        df = catalog.copy()
        df['sed_wave'] = 2.998e5/df['sed_freq']  # Converts from GHz to microns
        df.drop('sed_freq', axis=1, inplace=True)
        return df
    
    def _convert_coords_to_degrees(self, RA, dec):
        """
        Convert RA n h:m:s and dec in d:m:s to decimal degrees.
        Input can be list-like (e.g. RA = [0,42,30]), or strings (e.g. dec = '+41d12m00s' or
        '+41 12 00' or similar).
        """
        if not isinstance(RA, str):
            RA =  '{:d}h{:d}m{:f}s'.format(RA[0], RA[1], RA[2])
        if not isinstance(dec, str):
            dec = '{:+d}d{:d}m{:f}s'.format(dec[0], dec[1], dec[2])
        coords = SkyCoord(RA, dec)
        return coords.ra.deg, coords.dec.deg

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
    name_list = ['vega', 'ic348', 'not_a_source']
    
    vizier = VizierCatalog()
    pos_phot = vizier.get_all_dataframes(zip(ra_list, dec_list))
    named_phot = vizier.get_all_dataframes(name_list)

        
