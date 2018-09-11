import numpy as np
from nasa_tables import *
from astropy.io import ascii

if __name__ == "__main__":
    # OG Kepler long-cadence targets:
    print("making Kepler coords file...")
    lc_file = 'kic_lc_coords.csv'  
    select = 'kepid,ra,dec'
    lc_nasa_table = get_keplerstellar_table(select=select, cache=False)
    lc_nasa_table.rename_column('ra', 'ra_kic')
    lc_nasa_table.rename_column('dec', 'dec_kic')
    ascii.write(lc_nasa_table, lc_file, overwrite=True, format='csv')

    # K2 targets:
    print("making K2 coords file...")
    k2_file = 'epic_coords.csv'  
    select = 'epic_number,ra,dec'
    k2_nasa_table = get_k2targets_table(select=select, cache=False)
    k2_nasa_table.rename_column('ra', 'ra_epic')
    k2_nasa_table.rename_column('dec', 'dec_epic')
    ascii.write(k2_nasa_table, k2_file, overwrite=True, format='csv')
    
    # confirmed exoplanets:
    print("making confirmed exoplanets file...")
    exo_file = 'confirmed_coords.csv'  
    select = 'pl_hostname,pl_letter,ra,dec'
    exo_nasa_table = get_confirmed_planets_table(select=select, cache=False)
    exo_nasa_table.rename_column('ra', 'ra_nasa')
    exo_nasa_table.rename_column('dec', 'dec_nasa')
    exo_nasa_table.remove_columns(['pl_hostname', 'pl_letter'])
    ascii.write(exo_nasa_table, exo_file, overwrite=True, format='csv')