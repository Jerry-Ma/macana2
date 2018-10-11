'''
@package aztec
Python tools for aztec_c++ pipeline (codename: macana). Created on 04/09/2014

@author: David Sanchez-Arguelles
@license: GPL
@contact: domars@inaoep.mx
@requires: numpy scipy astropy 
'''

__version__="$Revision: 1.1 $"

import numpy as npy
import astropy as astro
import astropy.units as units
import astropy.constants as constants
import astropy.io.fits as pyfits
import astropy.coordinates as coordinates
import os

#Some constants

sigma2beam = 2.0*npy.sqrt(2.0*npy.log(2.0))
beam2sigma = 1.0/sigma2beam

LMT_position = coordinates.EarthLocation.from_geodetic (-97.31472222222222224,18.985, height=4640.)
ASTE_position = coordinates.EarthLocation.from_geodetic (-67.703304,-22.971657,  height=4861.9)


def plotmap(filename, wcut=70):
	from aztec.map import AztecMap
	from aztec.plots import MapPlotter
	
	m = AztecMap(filename)
	m.wcut(wcut)
	pl = MapPlotter(m)
	pl.plot()
	return pl
