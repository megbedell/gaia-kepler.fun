## gaia-kepler.fun

### how I made these matches

I've gone through a few iterations of cross-match techniques.

The original code for generating matches through CDS XMatch is contained in `query.py`. This code is now deprecated. I left off using CDS XMatch because it does not return the full Gaia DR2 data model. Around the same time, I extended the matching to include distances from [Bailer-Jones](https://arxiv.org/abs/1804.10121).

Here's how I do it now:

#### 1. generating coordinate lists from [Nasa Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/)

I use a custom code built from `astroquery` to make CSV files of coordinates for the following catalogs:

- Kepler long-cadence targets
- K2 targets
- confirmed planets

A script to get these CSV files is in `make_coord_files.py` with supplementary functions in `nasa_tables.py`, although note that subsequent changes to the NASA Exoplanet Archive interface have made these a little buggy.

#### 2. Gaia cross-matching with the [online archive](http://gea.esac.esa.int/archive/)

I upload the CSV files to the archive, then execute the following ADQL query:
```
SELECT * , distance(
  POINT('ICRS', kepler.ra_kic, kepler.dec_kic),
  POINT('ICRS', gaia.ra, gaia.dec)) AS angDist
FROM user_mbedell.kepler AS kepler
JOIN gaiadr3.gaia_source_lite AS gaia
  ON 1=CONTAINS(
    POINT('ICRS', kepler.ra_kic, kepler.dec_kic),
    CIRCLE('ICRS', gaia.ra, gaia.dec, 30./3600))
```
(the result of this query is saved as `data/kepler_30arcsec_gaia.fits`)

#### 3. cross-matching Gaia sources with Bailer-Jones distances -- NOT YET DONE FOR DR3

With each of the tables generated in the previous step, I open them in TOPCAT and do the following:

- trim the table down to contain only the Gaia `source_ids` (this will make the query go faster)

- connect to ARI-Gaia with TAP and execute the following query:
```
SELECT *
FROM TAP_UPLOAD.t2
JOIN gaiadr2_complements.geometric_distance USING (source_id)
```

- take the resulting table and save it as (e.g.) `data/kepler_30arcsec_dist.fits`

#### 4. combining tables into final data products

I then run the `build_tables_w_dist.py` script to trim unnecessary data columns, combine the three data sources (Gaia, Bailer-Jones, and NASA Exoplanet Archive), calculate angular separations between Gaia and NASA sources including propagation of proper motions between reference epochs (for Gaia sources with net proper motion at or above the level of 50 mas/yr), and save the appropriate match subsets.