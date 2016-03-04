#! /usr/bin/env python
# 
# Program: query_vizier_test
#
# Author: Nick Lee
#
# Usage: python -m unittest query_vizier_test
#
# Description: Unit Tests for the Vizier Query modules
#
# To Do:
#    
#

import unittest
import pandas as pd
import numpy as np
import pandas.util.testing as pdt
from query_vizier import VizierCatalog
import pdb

##################
### Unit Tests ###
##################
class TestQuery(unittest.TestCase):
    """ Tests for the querying of Vizier data """
    def setUp(self):
        self.vizier = VizierCatalog()

    ### Test for _replace_frequency_with_wavelength ###
    def test_freq_to_wave(self):
        """ Check conversion from frequency in GHz to wavelength in microns """
        input_frame = pd.DataFrame({'source_id':np.arange(5), 'sed_freq':[100, 0.01, 0, np.nan, -99]})
        new_frame = self.vizier._replace_frequency_with_wavelength(input_frame)
        wavelengths = pd.DataFrame({'source_id':np.arange(5), 'sed_wave':[2.997925e3, 2.997925e7, np.inf, np.nan, np.nan]})
        pdt.assert_series_equal(wavelengths['sed_wave'], new_frame['sed_wave'])
    
    ### Test for _create_url ###
    def test_create_url_name_input(self):
        """ Check _create_url with string name input"""
        url = self.vizier._create_url('vega')
        self.assertTrue(url=="http://vizier.u-strasbg.fr/viz-bin/sed?-c=vega&-c.rs=1.50")
    
    def test_create_url_tuple_input(self):
        """ Check _create_url with tuple of floats"""
        url = self.vizier._create_url((279.234733,38.783))
        self.assertTrue(url=="http://vizier.u-strasbg.fr/viz-bin/sed?-c=279.234733+38.783000&-c.rs=1.50")
    
    def test_create_url_zero_input(self):
        """ Check _create_url with tuple of 0 ints"""
        url = self.vizier._create_url((0,0))
        print url
        self.assertTrue(url=="http://vizier.u-strasbg.fr/viz-bin/sed?-c=0.000000+0.000000&-c.rs=1.50")
    
    def test_create_url_one_element_list(self):
        """ Check _create_url with tuple of floats"""
        self.assertRaises(IOError,self.vizier._create_url,(279.234733))
    
    def test_create_url_negative_dec(self):
        """ Check that negative declination is handled correctly """
        url = self.vizier._create_url((1.0,-38.783),radius=10.)
        self.assertTrue(url=="http://vizier.u-strasbg.fr/viz-bin/sed?-c=1.000000-38.783000&-c.rs=10.00")
    
    def test_create_url_negative_ra(self):
        """ Check that negative declination is handled correctly """
        url = self.vizier._create_url((-100.0,-38.783),radius=10.)
        self.assertTrue(url=="http://vizier.u-strasbg.fr/viz-bin/sed?-c=-100.000000-38.783000&-c.rs=10.00")
        

    ### Test for _check_coords ###
    def test_check_coords_conversion(self):
        """Check that various representations of Vega coords return the same result"""
        formats = ['tuples',
                    'lists',
                    'arrays',
                    'colon strings',
                    'signed coloncstrings',
                    'hms dms',
                    'hms dms with spaces',
                   'space strings']
        RA      = pd.Series([(18,36,56.3364),
                             [18,36,56.3364],
                    np.array([18,36,56.3364]),
                             '18:36:56.3364',
                            '+18:36:56.3364',
                             '18h36m56.3364s',
                             '18h 36m 56.3364s',
                             '18 36 56.3364'], index=formats, name='RA')
        dec     = pd.Series([(38,47,1.291),
                             [38,47,1.291],
                    np.array([38,47,1.291]),
                             '38:47:1.291',
                            '+38:47:1.291',
                             '38d47m1.291s',
                             '38d 47m 1.291s',
                             '38 47 1.291'],   index=formats, name='dec')
        tests   = pd.concat([RA, dec], axis=1)  #Concat Series to a DataFrame
        result  = (279.234735, 38.783691944)    #This should be the result of all calls
        for r,d in zip(tests['RA'], tests['dec']):
            current_result = self.vizier._check_coords((r,d))
            pdt.assert_almost_equal(current_result, result)
