import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

import logging
logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(levelname)s- %(message)s')

def check_coords(position):
    """
    Check that a source position is in the correct format.
    If not, change to default coords.

    Args:
        position - source position. Must be 2-element list-like object with
                   position[0] = RA, position[1] = dec.
                   Acceptable formats for RA & dec are decimal degrees as
                   single floats or ints, or sexagisimal as 3-element
                   list-like object of (h,m,s), or strings with h/m/s.

    Returns:
        source_pos - Tuple of (RA,dec) in degrees.
    """
    in_degrees = isinstance(position[0], (float, int)) & isinstance(position[1], (float, int))
    if isinstance(position[0], (np.ndarray, list)):
        position = (tuple(position[0]), tuple(position[1]))
    if not in_degrees:
        try:
            coords = SkyCoord(position[0], position[1], unit=(u.hourangle, u.deg))
            new_pos = (coords.ra.deg, coords.dec.deg)
        except ValueError as exc:
            print exc
            raise SystemExit
        except:
            logging.critical('WARNING: RA and dec format cannot be read')
            raise
    else:
        new_pos = (float(position[0]), float(position[1]))

    if not (0 <= new_pos[0] <= 360):
        raise ValueError('RA should be in the range [0,360]')
    if not (-90 <= new_pos[1] <= 90):
        raise ValueError('dec should be in the range [-90,90]')

    return new_pos