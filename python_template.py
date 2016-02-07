#! /usr/bin/env python
# 
# Program: python_template
#
# Author: Nick Lee
#
# Usage: ./python_template
#
# Description: Examples of header, docstrings, etc. to use in python code
#
# To Do:
#    Actually write some code. 
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

## Debugging module and intital setup
import logging
logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(levelname)s- %(message)s')

###############
### Classes ###
###############
class Sample_Class(object):
    """Summary of class here.

    Longer class information....
    Longer class information....

    Attributes:
        likes_spam: A boolean indicating if we like SPAM or not.
        eggs: An integer count of the eggs we have laid.
    """

    def __init__(self, likes_spam=False):
        """Inits SampleClass with blah."""
        self.likes_spam = likes_spam
        self.eggs = 0

    def public_method(self):
        """Performs operation blah."""
        pass

#################
### Functions ###
#################
def fetch_bigtable_rows(big_table, keys, other_silly_variable=None):
    """ Single sentence summary. e.g. Fetches rows from a Bigtable.

    Longer summary here. E.g.
    Retrieves rows pertaining to the given keys from the Table instance
    represented by big_table.  Silly things may happen if
    other_silly_variable is not None.

    Args:
        big_table: An open Bigtable Table instance.
        keys: A sequence of strings representing the key of each table row
            to fetch.
        other_silly_variable: Another optional variable, that has a much
            longer name than the other args, and which does nothing.

    Returns:
        A dict mapping keys to the corresponding table row data
        fetched. Each row is represented as a tuple of strings. For
        example:

        {'Serak': ('Rigel VII', 'Preparer'),
         'Zim': ('Irk', 'Invader'),
         'Lrrr': ('Omicron Persei 8', 'Emperor')}

        If a key from the keys argument is missing from the dictionary,
        then that row was not found in the table.

    Raises:
        IOError: An error occurred accessing the bigtable.Table object.
    """
    pass

####################
####################
def parse_args():
    '''
    Read command line arguments
    '''
    ## Define a parser
    parser = argparse.ArgumentParser(description='Example python program')

    ## Add arguments
    parser.add_argument('-d','--debug',action='store_true',help='Turn on debugging messages')

    ## Use parser to interpret the command line arguments
    inputs = vars(parser.parse_args())

    ## Set debugging options
    if inputs['debug']==False:
        logging.disable(logging.CRITICAL)  # Turns off debug messages
    else:
        logging.disable(logging.NOTSET)
        logging.debug('Start of program')

#####################
###### Run Code #####
#####################
if __name__ == '__main__':
    ## Process command line arguments
    parse_args()

    
