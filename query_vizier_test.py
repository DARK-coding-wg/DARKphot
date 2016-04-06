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
import os
import pandas as pd
import numpy as np
import pandas.util.testing as pdt
from query_vizier import VizierCatalog
import pdb

import logging
logger = logging.getLogger()
logger.setLevel(logging.WARNING)
logging.disable(logging.CRITICAL)  # Uncomment to turn off error messages

##################
### Unit Tests ###
##################
class TestQuery(unittest.TestCase):
    """ Tests for the querying of Vizier data """
    def setUp(self):
        self.vizier = VizierCatalog()

    def assert_equal_pos(self, pos1, pos2, delta=None):
        """ Check two positions and make sure they're the same within a certain distance """
        if delta is None:
            delta = 0.01/3600.  # 1/10 arcsec
        pos1 = self.vizier._check_coords(pos1)
        pos2 = self.vizier._check_coords(pos2)
        self.assertAlmostEqual(pos1[0], pos2[0], delta=delta)
        self.assertAlmostEqual(pos1[1], pos2[1], delta=delta)

    ### Test for _replace_frequency_with_wavelength ###
    def test_freq_to_wave(self):
        """ Check conversion from frequency in GHz to wavelength in microns """
        input_frame = pd.DataFrame({'source_id': np.arange(5), 'sed_freq': [100, 0.01, 0, np.nan, -99]})
        new_frame = self.vizier._replace_frequency_with_wavelength(input_frame)
        wavelengths = pd.DataFrame({'source_id': np.arange(5), 'sed_wave': [2.997925e3, 2.997925e7, np.inf, np.nan, np.nan]})
        pdt.assert_series_equal(wavelengths['sed_wave'], new_frame['sed_wave'])

    ## Tests for _read_vo_table
    def test_read_vega_table(self):
        """ Make sure the vega table is being read correctly """
        vega_path = 'test_data/vega_vizier_votable.vot'
        if not os.path.exists(vega_path):
            vega_url = 'http://vizier.u-strasbg.fr/viz-bin/sed?-c=Vega&-c.rs=1.5'
            self.vizier._download_from_vizier(vega_url, vega_path)
        vega_table = self.vizier._read_vo_table(vega_path)
        vega_pos = (279.23473333333334, 38.78368888888889)
        self.assert_equal_pos((vega_table['_RAJ2000'].mean(), vega_table['_DEJ2000'].mean()), vega_pos, delta=1.5/3600.)
    def test_read_non_existant_file(self):
        """ Check to make sure the right Exceptions are raised if trying to read a non-existant filename"""
        bad_path = 'test_data/non_existant_file.vot'
        if os.path.exists(bad_path):
            os.remove(bad_path)
        with self.assertRaises(IOError):
            self.vizier._read_vo_table(bad_path)
    def test_empty_file(self):
        """ Check to make sure the right Exceptions are raised if trying to read a votable with no data """
        empty_path = 'test_data/empty_votable.vot'
        if not os.path.exists(empty_path):
            empty_url = 'http://vizier.u-strasbg.fr/viz-bin/sed?-c=not_a_source&-c.rs=1.5'
            self.vizier._download_from_vizier(empty_url, empty_path)
        with self.assertRaises(ValueError):
            empty_table = self.vizier._read_vo_table(empty_path)

    ### Test for _create_url ###
    def test_create_url_name_input(self):
        """ Check _create_url with string name input"""
        url = self.vizier._create_url('vega')
        self.assertTrue(url == "http://vizier.u-strasbg.fr/viz-bin/sed?-c=vega&-c.rs=1.50")
    def test_create_url_tuple_input(self):
        """ Check _create_url with tuple of floats"""
        url = self.vizier._create_url((279.234733, 38.783))
        self.assertTrue(url == "http://vizier.u-strasbg.fr/viz-bin/sed?-c=279.234733+38.783000&-c.rs=1.50")
    def test_create_url_zero_input(self):
        """ Check _create_url with tuple of 0 ints"""
        url = self.vizier._create_url((0, 0))
        self.assertTrue(url == "http://vizier.u-strasbg.fr/viz-bin/sed?-c=0.000000+0.000000&-c.rs=1.50")
    def test_create_url_negative_dec(self):
        """ Check that negative declination is handled correctly """
        url = self.vizier._create_url((1.0, -38.783), radius=10.)
        self.assertTrue(url == "http://vizier.u-strasbg.fr/viz-bin/sed?-c=1.000000-38.783000&-c.rs=10.00")

    ### Test for _check_coords ###
    def test_check_coords_formats(self):
        """Check that various representations of Vega coords return the same result"""
        formats = ['tuples',
                    'lists',
                    'arrays',
                    'colon strings',
                    'signed coloncstrings',
                    'hms dms',
                    'hms dms with spaces',
                    'space strings']
        ra      = pd.Series([(18,36,56.3364),
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
                             '38 47 1.291'], index=formats, name='dec')
        tests = pd.concat([ra, dec], axis=1)  #Concat Series to a DataFrame
        result = (279.234735, 38.783691944)  #This should be the result of all calls
        for rr, dd in zip(tests['RA'], tests['dec']):
            current_result = self.vizier._check_coords((rr, dd))
            pdt.assert_almost_equal(current_result, result)

    def test_check_coords_ValueError(self):
        """Test that wacky numbers and numbers outside sensible range raise an error"""
        RAfail = [np.nan, -1, 361, np.inf, 1e30, (25, 30, 0), 'wallah'] #, (-12,30,0)] <-- Okay so RA = -12 is allowed. It shouldn't be...
        decfail = [np.nan, -91, 91, np.inf, 1e30, (91, 30, 0), (-91, 30, 0), 'billah']
        for ra in RAfail:
            with self.assertRaises(ValueError):
                self.vizier._check_coords((ra, 30.))
        for dec in decfail:
            with self.assertRaises(ValueError):
                self.vizier._check_coords((30., dec))

    def test_check_coords_TypeError(self):
        """Test that non-sensible strings in tuples and Nones raise an error"""
        RAfail = [('wallah', 30, 0), None, (None, 30, 0)]
        decfail = [('billah', 30, 0), None, (None, 30, 0)]
        for ra in RAfail:
            with self.assertRaises(TypeError):
                self.vizier._check_coords((ra, 30.))
        for dec in decfail:
            with self.assertRaises(TypeError):
                self.vizier._check_coords((30., dec))

    def test_check_coords_NameError(self):
        """Test that undefined variables raise an error"""
        with self.assertRaises(NameError):
            self.vizier._check_coords(((wallah, 30, 0), 30.))
        with self.assertRaises(NameError):
            self.vizier._check_coords((wallah, 30.))
        with self.assertRaises(NameError):
            self.vizier._check_coords((30., (billah, 30, 0)))
        with self.assertRaises(NameError):
            self.vizier._check_coords((30., billah))

    def test_check_coords_ValidInput(self):
        """Test that edge cases work (is 'edge case' the wrong term here?)"""
        ra_okay = [0, 1e-30, (1e-30, 30, 0)]
        dec_okay = [0, 1e-30, (1e-30, 30, 0)]
        for ra in ra_okay:
            self.vizier._check_coords((ra, 30.))
        for dec in dec_okay:
            self.vizier._check_coords((30., dec))

    def test_download_from_vizier(self):
        """ Check the download of a url """
        vega_url = 'http://vizier.u-strasbg.fr/viz-bin/sed?-c=Vega&-c.rs=1.5'
        vega_path = 'test_data/test_download_from_vizier.vot'
        self.vizier._download_from_vizier(vega_url, vega_path)
        self.assertTrue(os.path.exists(vega_path))
        try:
            os.remove(vega_path)
        except:
            print vega_path + " didn't download properly."

    def test_download_vega_table(self):
        """ Make sure the vega table is being downloaded correctly """
        vega_path = 'test_data/vega_downloaded_table.vot'
        vega_url = 'http://vizier.u-strasbg.fr/viz-bin/sed?-c=Vega&-c.rs=1.5'
        self.vizier._download_from_vizier(vega_url, vega_path)
        vega_table = self.vizier._read_vo_table(vega_path)
        vega_pos = (279.23473333333334, 38.78368888888889)
        self.assert_equal_pos((vega_table['_RAJ2000'].mean(), vega_table['_DEJ2000'].mean()), vega_pos, delta=1.5/3600.)
        try:
            os.remove(vega_path)
        except:
            print vega_path + " didn't download properly"

    ### Test get_all_dataframes ###
    def test_source_id_not_provided(self):
        source_list = [(187.27832916, 2.05199), 'not_a_source', 'Vega']
        phot = self.vizier._get_all_dataframes(source_list, source_id=range(len(source_list)), radius=1.5)
        self.assertTrue(all(phot['source_id'].unique() == [0, 2]))

    def test_source_id_provided(self):
        source_list = [(187.27832916, 2.05199), 'not_a_source', 'Vega']
        source_id = [-6, '23', 4.5]
        phot = self.vizier._get_all_dataframes(source_list, source_id=source_id, radius=1.5)
        for ii, src_id in enumerate([-6, 4.5]):
            self.assertAlmostEqual(src_id, phot['source_id'].unique()[ii])

    ### Test query_vizier ###
    def test_query_vizier_source_id(self):
        source_list = [(187.27832916, 2.05199), 'not_a_source', ('18h 36m 56.3364s', '+38:47:1.291'), (1, 2, 3), 'Vega']
        data, source_info = self.vizier.query_vizier(source_list)
        self.assertTrue(all(source_info['source_id'] == [0, 1, 2, 3, 4]))

    def test_query_vizier_sources(self):
        source_list = [(187.27832916, 2.05199), 'not_a_source', ('18h 36m 56.3364s', '+38:47:1.291'), (1, 2, 3), 'Vega']
        data, source_info = self.vizier.query_vizier(source_list)
        for ii, src_id in enumerate(source_list):
            self.assertAlmostEqual(src_id, source_info.loc[ii, 'source'])

    def test_query_vizier_data(self):
        source_list = [(187.27832916, 2.05199), 'not_a_source', ('18h 36m 56.3364s' ,'+38:47:1.291'), (1, 2, 3), 'Vega']
        data, source_info = self.vizier.query_vizier(source_list)
        self.assertTrue(all(data['source_id'].unique() == [0, 2, 4]))







