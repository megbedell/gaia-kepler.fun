import numpy as np
from astropy import units as u
from astropy.io import ascii
from astroquery.xmatch import XMatch

def xmatch(file, ra_name='RA', dec_name='Dec'):
    # do a CDS XMatch between targets in file and Gaia; return astropy table
    gaia_cat = 'vizier:I/337/tgas' # DR1 TGAS
    table = XMatch.query(cat1=open(file), cat2=gaia_cat, max_distance=5 * u.arcsec, colRA1=ra_name, colDec1=dec_name)
    return table
    
def save_output(table, filename):
    # take astropy table and write out
    ascii.write(table, filename)
    
if __name__ == "__main__":
    filename = '/Users/mbedell/Documents/Research/Stars/Kepler/kepler_sc_targets.csv'
    ra, dec = np.genfromtxt(filename, usecols=(5,6), delimiter=',', unpack=True, dtype=None, skip_header=2)
    # TODO: convert sexagesimal to decimal & rewrite file
    table = xmatch('/Users/mbedell/Documents/Research/Stars/Kepler/kepler_sc_targets.csv', 
                    ra_name='RA (J2000)', dec_name='Dec (J2000)')