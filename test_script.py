import sys
import os
import glob

import numpy as np

from astropy.io import fits
from astropy.io import ascii
from astropy.table import QTable

from astropy.cosmology import FlatLambdaCDM

import astropy.units as u

import pyspeckit

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style as style

import pylab as pl
for ii in pl.get_fignums():
    pl.close(ii)

datafile = fits.open('PKS1017-325_heliocorrect_galextinct.fits')
n = len(datafile[0].data)
crpix1 = datafile[0].header['CRPIX1']
crval1 = datafile[0].header['CRVAL1']
cd1_1 = datafile[0].header['CD1_1'] # The data have been resampled to 0.85 angstroms / pixel
wavelength = ((np.arange(n) + 1.0) - crpix1) * cd1_1 + crval1
flux = datafile[0].data
errorfile = fits.open('PKS1017-325_error_heliocorrect_galextinct.fits')
error = errorfile[0].data

sp = pyspeckit.Spectrum(data=flux, xarr=wavelength, error=error, unit='erg/s/cm^2/AA', xarrkwargs={'unit':'AA'},header=datafile[0].header)
#
sp.plotter(xmin=6050, xmax=6750, errstyle='fill')
sp.plotter.axis.set_xlabel(r'Wavelength $(\AA)$')
sp.plotter.axis.set_ylabel(r'Flux $(\mathrm{erg/s/cm^2/\AA})$')
sp.plotter.refresh()
# Clear line detections in this window, from blue to red: 

o3region_linewavs = {'HeII' : 4685.7,  
                     'Hbeta' : 4861.34, 
                     'OIII4959' :4958.9,
                     'OIII5007' : 5006.8
                    }


#o3region_linewavs = {  
#                     'Hbeta' : 4861.34, 
#                     'OIII4959' :4958.9,
#                     'OIII5007' : 5006.8
#                    }


z_guess = 0.31868
z_guess_factor = (1 + z_guess)

fwhm_guess = 15.

# I just quickly estimated this in IRAF SPLOT.
# I could also estimate this in python, with 
# width_guess = data.sum() / amplitude_guess / np.sqrt(2*np.pi)

o3region_amplitudes = {}
for line in o3region_linewavs:
    sp.specfit.selectregion(xmin = o3region_linewavs[line] *z_guess_factor - 30, xmax = o3region_linewavs[line]*z_guess_factor + 30)
    amp = np.max(sp.data[sp.specfit.xmin:sp.specfit.xmax])
    d = {line : amp}
    o3region_amplitudes.update(d)
    
o3region_redshifted_lines = {}
for line in o3region_linewavs:
    redline = o3region_linewavs[line] * z_guess_factor
    d = {line : redline}
    o3region_redshifted_lines.update(d)    
    

o3region_redshifted_lines

guesses = []
# Remember, dictonaries aren't ordered. Force it to be with sorted() so that we go blue --> red
for line in sorted(o3region_linewavs):
    guess = [o3region_amplitudes[line], o3region_redshifted_lines[line], fwhm_guess]
    print(line)
    print(guess)
    guesses += guess
    
guesses


linemask = []
for line in o3region_linewavs:
    xmin = o3region_linewavs[line] - 30.0
    xmax = o3region_linewavs[line] + 30.0
    mask = [xmin, xmax]
    linemask += mask
    
linemask

# Then fit baseline with
sp.baseline(xmin=6050, xmax=6750, exclude=linemask, subtract=False, reset_selection=False, highlight_fitregion=True, order=1)



sp.plotter()
sp.specfit(guesses=list(guesses), fittype='gaussian', xmin=6100, xmax=6700)
print(sp.specfit.parinfo)
#sp.specfit.parinfo['SHIFT1'].value = 6178
#sp.specfit.parinfo['SHIFT1'].limited = (True,True)
#sp.specfit.parinfo['SHIFT1'].limits = (6150, 6200)
#sp.specfit.parinfo['WIDTH1'].value = 5
#print(sp.specfit.parinfo)
#sp.specfit(parinfo=sp.specfit.parinfo, fittype='gaussian', xmin=6100, xmax=6700, guesses=None)
#print(sp.specfit.parinfo)
