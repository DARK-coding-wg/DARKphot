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

##################
### Unit Tests ###
##################
class TestQuery(unittest.TestCase):
    """ Tests for the querying of Vizier data
    """
    def setUp(self):
        self.vizier = VizierCatalog()
        
    def test_freq_to_wave(self):
        """ Check conversion from frequency in GHz to wavelength in microns """
        input_frame = pd.DataFrame({'source_id':np.arange(5), 'sed_freq':[100, 0.01, 0, np.nan, -99]})
        new_frame = self.vizier._replace_frequency_with_wavelength(input_frame)
        wavelengths = pd.DataFrame({'source_id':np.arange(5), 'sed_wave':[2.997925e3, 2.997925e7, np.inf, np.nan, np.nan]})
        pdt.assert_series_equal(wavelengths['sed_wave'], new_frame['sed_wave'])
    
    
    """ Test for _create_url """
    def test_create_url_name_input(self):
        """ Check _create_url with string name input"""
        vizier = VizierCatalog()
        url_return = vizier._create_url('vega')
        self.assertTrue(url_return=="http://vizier.u-strasbg.fr/viz-bin/sed?-c=vega&-c.rs=1.50")
    
    def test_create_url_tuple_input(self):
        """ Check _create_url with tuple of floats"""
        vizier = VizierCatalog()
        url_return = vizier._create_url((279.234733,38.783))
        self.assertTrue(url_return=="http://vizier.u-strasbg.fr/viz-bin/sed?-c=279.234733+38.783000&-c.rs=1.50")
    
    def test_create_url_one_element_list(self):
        """ Check _create_url with tuple of floats"""
        vizier = VizierCatalog()
        self.assertRaises(IOError,vizier._create_url,(279.234733))