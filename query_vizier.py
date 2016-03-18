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
import warnings
import requests
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import pdb

## Debugging module and intital setup
import logging
logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(levelname)s- %(message)s')

###############
### Classes ###
###############
class VizierCatalog(object):
    """ Class for creating a photometric catalog from Vizier """

    def query_vizier(self, source_list, radius=1.5):
        """ Creates photometric table and source information table

        Args:
            source_list - list of souces. Can be a list of (ra, dec) tuples or a list of names as strings. RA & Dec must be in decimal degrees.
        Keywords:
            radius - Defines region to search for counterparts in arcsec. Default is 1.5
        Returns:
            vizier_data - pandas dataframe containing all observations from Vizier and an additional source_id column that identifies which source each obsersvation is for
            source_info - pandas dataframe containing the names, positions, and source_id of each input source. 
        """
        ## Figure out if inputs are names or positions
        names = []
        pos = []
        source_id = []
        
        for ii, source in enumerate(source_list):
            ## Single string = Source name
            if isinstance(source, str):
                names.append(source)
                pos.append((np.nan, np.nan))
                source_id.append(ii)
            else:
                names.append(np.nan)
                pos.append(self._check_coords(source))
                source_id.append(ii)
        
        ## Create pandas dataframe containing source_id & source name and/or position
        self.source_info = pd.DataFrame({'names':names, 'positions':pos, 'source_id':source_id, 'input':source_list})
        ## Go through Vizier queries 
        vizier_data = self.get_all_dataframes(source_list, source_id=source_id, radius=radius)
        return vizier_data, self.source_info

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
        vot_name = 'test_vo_table.vot'
        for ii, source in enumerate(source_list):
            url = self._create_url(source, radius=radius)
            self._download_from_vizier(url, vot_name)
            try:
                photometry = self._read_vo_table(vot_name)
            except ValueError:
                logging.warning('No observations found via Vizier for source_id = {:}'.format(source_id[ii]))
            else:
                n_rows = photometry.shape[0]
                photometry['source_id'] = np.repeat(source_id[ii], n_rows)
                list_of_frames.append(photometry)
            finally:
                os.remove(vot_name)
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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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
        if isinstance(source, str):
            url = 'http://vizier.u-strasbg.fr/viz-bin/sed?-c={:}&-c.rs={:4.2f}'.format(source, float(radius))
        else:
            try:
                coords = self._check_coords(source)
            except:
                raise IOError('Provided Source is not of correct type. Must be either a tuple of (ra, dec) or a string of name')
            else:
                url = 'http://vizier.u-strasbg.fr/viz-bin/sed?-c={:f}{:+f}&-c.rs={:4.2f}'.format(coords[0], coords[1], float(radius))
        return url

    def _replace_frequency_with_wavelength(self, catalog):
        """ Converts the column sed_freq to sed_wave with the appropriate units.

        Treats any negative frequencies as non-detections and replaces them with NaNs """
        df = catalog.copy()
        speed_of_light = 299792.458  # [um * GHz]
        df['sed_wave'] = speed_of_light/df['sed_freq']  # Converts from GHz to microns
        df.loc[df['sed_wave'] < 0, 'sed_wave'] = np.nan
        df.drop('sed_freq', axis=1, inplace=True)
        return df

    def _check_coords(self, position):
        """ Check that a source position is in the correct format. If not change to default coords.

        Args:
            position - source position. Must be 2-element list-like object with position[0] = ra, position[1] = dec.
                       Acceptable formats for ra & dec are decimal degrees as single floats or ints, or
                       sexagisimal as 3-element tuple/lists of (h, m, s), or strings with h/m/s
        Returns:
            source_pos - Tuple of (ra,dec) in degrees.
        """
        in_degrees = isinstance(position[0], (float, int)) & isinstance(position[1], (float, int))
        
        if isinstance(position[0],(np.ndarray, list)):
             position = (tuple(position[0]), tuple(position[1]))
        if not in_degrees:
            try:
                coords = SkyCoord(position[0], position[1], unit=(u.hourangle, u.deg))
                new_pos = (coords.ra.deg, coords.dec.deg)
            except:
                logging.critical('WARNING: ra and dec format cannot be read')
                raise
        else:
            new_pos = (float(position[0]), float(position[1]))
        return new_pos
        

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
    ra_list = [187.27832916, 150.231314, 149.426254, '15 34 57.224'] 
    dec_list = [2.05199, -4.1234,  2.073906, '+23 30 11.610'] 
    #ra_list = ['14:23:45.45'] #
    #dec_list = ['02:10:45.45'] #
    name_list = ['vega', 'ic348', 'not_a_source']
    
    vizier = VizierCatalog()
    pos_phot, pos_sources = vizier.query_vizier(zip(ra_list, dec_list))
    named_phot, named_sources = vizier.query_vizier(name_list)



        
