import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits
from astroquery.xmatch import XMatch
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column, MaskedColumn
from tqdm import tqdm

def xmatch(file, ra_name='RA', dec_name='Dec', dist=5, cat='src'):
    # do a CDS XMatch between targets in file and Gaia; return astropy table
    if cat=='tgas':
        gaia_cat = 'vizier:I/337/tgas' # DR1 TGAS
    elif cat=='src':
        gaia_cat = 'vizier:I/337/gaia'
    table = XMatch.query(cat1=open(file), cat2=gaia_cat, max_distance=dist * u.arcsec, colRA1=ra_name, colDec1=dec_name)
    return table
    
def save_output(table, filename):
    # take astropy table and write out
    table.write(filename, format='fits')
    
if __name__ == "__main__":
    # OG Kepler long-cadence targets:
    lc_file = 'kic_lc_coords.csv'  
     
    try: # file exists
        lc = np.genfromtxt(lc_file, delimiter=',', unpack=True, dtype=None, skip_header=1)
        kic, ra, dec = lc['f0'], lc['f1'], lc['f2'] 
    except: # make the file
        hdu = fits.open('/Users/mbedell/Documents/Research/Stars/Kepler/kepler_lc_targets.fit')
        data = hdu[1].data
        ra, dec = data['_RA'], data['_DE']
        kic = data['KIC']
        f = open(lc_file, 'w')
        f.write('KIC, RA, Dec\n')
        for k,(r,d) in zip(kic, zip(ra, dec)):
            f.write('{0}, {1:.8f}, {2:.8f}\n'.format(k, r, d))
        f.close()
        print('file made')
    
    dist = 3
    cat = 'src'
    lc_table = xmatch(lc_file, ra_name='RA', dec_name='Dec', dist=dist, cat=cat)
    save_output(lc_table, 'lc_{c}_{d}arcsec.fits'.format(d=dist, c=cat))
    
    lc_matches = [np.sum(lc_table['KIC'] == k) for k in tqdm(kic)]
    plt.hist(lc_matches)
    plt.ylabel('Sources')
    plt.xlabel('Matches')
    plt.savefig('lc_{c}_matches_{d}arcsec.png'.format(d=dist, c=cat))
    plt.clf()
    
    plt.hist(lc_table['angDist'], bins=np.arange(0.,dist,0.1))
    plt.ylabel('Sources')
    plt.xlabel('Angular Dist (arcsec)')
    plt.savefig('lc_{c}_angdist_{d}arcsec.png'.format(d=dist, c=cat))
    plt.clf()
    
    
    # OG Kepler short-cadence targets:
    sc_file = '/Users/mbedell/Documents/Research/Stars/Kepler/kepler_sc_kicnames.txt'
    sc_names = np.genfromtxt(sc_file, delimiter=',', unpack=True, dtype='int64')
    sc_table = lc_table.copy()
    delete = []
    for i in tqdm(range(len(sc_table))):
        if not sc_table['KIC'][i] in sc_names:
            delete.append(i)
    sc_table.remove_rows(delete)
    save_output(lc_table, 'sc_{c}_{d}arcsec.fits'.format(d=dist, c=cat))
    
    sc_matches = [np.sum(sc_table['KIC'] == k) for k in tqdm(sc_names)]
    plt.hist(sc_matches)
    plt.ylabel('Sources')
    plt.xlabel('Matches')
    plt.savefig('sc_{c}_matches_{d}arcsec.png'.format(d=dist, c=cat))
    plt.clf()
    
    
    plt.hist(sc_table['angDist'], bins=np.arange(0.,dist,0.1))
    plt.ylabel('Sources')
    plt.xlabel('Angular Dist (arcsec)')
    plt.savefig('sc_{c}_angdist_{d}arcsec.png'.format(d=dist, c=cat))
    plt.clf()
    
    
        
