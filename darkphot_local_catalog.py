import numpy as np
from types import *
############################
# Main class: LocalCatalog #
############################

# TODO input checks

class LocalCatalog(object):
    """
    Class for creating a local catalog in DARKphot
    """

    def add_local_catalog(self, fname_input_catalog, fname_parameter_file,
                          fname_object_sel_names=None,
                          name_ra_col = 'RA', name_dec_col='Dec',
                          name_id_col=None):
        """

        Args:
            fname_input_catalog (str): Filename of the input catalog
            fname_parameter_file (str): Filename of the parameter file,
                specifying the passband which need to be imported
            fname_object_sel_names (str): Filename including object names which shall
                be taken from the catalog
            name_ra_col (str; default: RA): Name of the RA column
                in the input catalog
            name_dec_col (str; default: Dec): Name of the Dec column
                in the input catalog
            name_id_col (str; default: None): Name of the source identifer
                column in the input catalog
        """

        import pandas as pd

        self.fname_parameter_file = fname_parameter_file
        self.fname_input_catalog = fname_input_catalog
        self.fname_object_sel_names = fname_object_sel_names

        col_flux_names, col_fluxerror_names, central_wavelengths = \
            self._read_parameter_file()

        # Checking for the file_format (currently based on the file ending)
        file_format = fname_input_catalog.split('.')[-1]

        if file_format in ['csv']:
            inter_catalog = self._read_text_catalog(file_format)
        else:
            raise(ValueError('fname_input_catalog has not an accpeted file ending. (Allowed: csv)'))

        if self.fname_object_sel_names is not None:

            if self._read_sel_object_names(): # Returns False, if there is no object
                mask = get_match_mask_names(self.object_sel_names,
                                            inter_catalog[name_id_col])
                inter_catalog = inter_catalog[mask]

        for col_name_flux, col_name_fluxerror, central_wavelength in \
                zip(col_flux_names, col_fluxerror_names, central_wavelengths):
            fluxes = inter_catalog[col_name_flux]
            flux_errors = inter_catalog[col_name_fluxerror]
            ra_values, dec_values = inter_catalog[name_ra_col], inter_catalog[name_dec_col]
            id_names = inter_catalog[name_id_col]
            #freq = 2.998E18 / float(central_wavelength) # Central wavelength is in AA
            cen_wavelenght = central_wavelength * 1.e-4 # Wavelength from AA to micron

            #Creating dataframe for unique names of objects
            names_frame = pd.DataFrame({'name':id_names})
            names_frame.drop_duplicates

            # Extract IDs from names frame
            IDs = names_frame.name[names_frame.name == id_names].index.tolist()

            inter_data_frame = pd.DataFrame({'_ID':IDs,
                                             '_RAJ2000':ra_values,
                                             '_DEJ2000':dec_values,
                                             #'sed_freq':len(ra_values)*[freq],
                                             'sed_wave':cen_wavelenght,
                                             'sed_flux':fluxes,
                                             'sed_eflux':flux_errors})

            if hasattr(self, 'data_frame'):
                self.data_frame = pd.concat([self.data_frame, inter_data_frame])
            else:
                self.data_frame = inter_data_frame

        # TODO Need to think more about the index creation and about the cross reference table between index and object
        # identifier
        self.data_frame.index = range(len(self.data_frame))

    def _read_parameter_file(self):
        """ Read in the parameter file (private) """

        # TODO: Needs checking of potential problems with the input file
        extension_name = self.fname_parameter_file.split(".")[-1]
        assert extension_name == "param", \
            "Input file extension needs to be .param. Currently it is ."+str(extension_name)

        with open(self.fname_parameter_file, 'r') as f:

            col_flux_names = []
            col_fluxerror_names = []
            central_wavelengths = []
            ii = 0
            for line in f:
                if line[0] == '#':
                    ii += 1
                    continue
                else:
                    ii += 1
                    inter = line.strip(' \n\t\n')
                    inter = inter.split()

                    try:
                        float(inter[0])
                    except ValueError:
                        raise ValueError("Input type of 1st column at line "+str(ii)+" contains only digits. It should be a string")
                        
                    col_flux_names.append(inter[0])
                    col_fluxerror_names.append(inter[1])
                    central_wavelengths.append(float(inter[2]))


        return col_flux_names, col_fluxerror_names, central_wavelengths

    def _read_sel_object_names(self):
        """
        Reading in the list of object names which shall be put into the catalog
        Returns:

        """
        # TODO: Needs checking of potential problems with the input file
        with open(self.fname_object_sel_names, 'r') as f:
            self.object_sel_names = []
            name_sel = False
            for line in f:
                if line[0] == '#':
                    continue
                else:
                    inter = line.strip(' \n\t\n')
                    self.object_sel_names.append(inter)
                    name_sel = True # Only if there is a line, in the file specifying an object the selection will be activated

        return name_sel

    def _read_text_catalog(self, file_format):
        """ Read in the catalog (private) """

        import astropy.io.ascii as asciitable

        # Converting potential file endings to asciitable input identifiers
        file_format_str_list = ['csv']
        asciitable_format_str_list = ['csv']
        asciitable_format = \
            asciitable_format_str_list[file_format_str_list.index(file_format)]

        inter_catalog = asciitable.read(self.fname_input_catalog,
                                        format=asciitable_format)

        return inter_catalog

#####################
### Help functions ##
#####################

def get_match_mask_names(sel_list, master_array):
    """
    Create a mask for those elements in a master_array which are included in a
    a selection list
    Args:
        sel_list (list):
        master_array (array):

    Returns:

    """
    # TODO Think about how this function can be made more efficient

    master_mask = np.zeros(len(master_array), dtype=bool)
    for element in sel_list:
        master_mask = master_mask | (master_array == element)

    return master_mask





#####################
###### Run Code #####
#####################
if __name__ == '__main__':
    ## Process command line arguments

    ## A couple sources, including one that doesn't have any counterparts
    fname_catalog = 'example_data/test_cat.csv'
    fname_param = 'example_data/input_local_cat.param'
    fname_object_names_sel = 'example_data/input_names.txt'

    lc1 = LocalCatalog()
    lc1.add_local_catalog(fname_catalog, fname_param, fname_object_names_sel, name_id_col='NAME')
    print lc1.data_frame

