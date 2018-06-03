### gaia-kepler.fun

This repo contains source code and demos for the Gaia-Kepler crossmatches at [gaia-kepler.fun](http://gaia-kepler.fun).

The website is built on [Bootstrap](https://getbootstrap.com/) using a slightly modified version of the theme [Creative](https://github.com/BlackrockDigital/startbootstrap-creative).

#### contents

Various utility functions for cross-matching and data munging written in python are located under `scripts`. A full explanation of how the cross-matches were done is in `scripts/HOWTO.md`. 

A demo Jupyter notebook which reads in the cross-match and makes cool plots (including color-magnitude diagrams displayed on the website, interactive bokeh versions, and a fun animation) can be found at `notebooks/demo.ipynb`.