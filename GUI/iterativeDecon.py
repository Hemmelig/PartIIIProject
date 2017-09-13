#!/usr/bin/env python3

import obspy
from obspy import read
from obspy.core import Stream
#imports
from obspy.core import event
from obspy.taup.taup import getTravelTimes
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import os.path
import time
import glob
import shutil
import numpy as np
import scipy
import cmath
from pylab import semilogy
import warnings
import time

def iterative_deconvolution(comp2, comp1, maxbumps=200, dt=0.04, filt='cosine', fmax=0.5, timeshift=25.):
    '''
Ligorria & Ammon 1999, Kikuchi & Kanamori 1982
After Chuck Ammons code 1998- Jenny Andrews 2006
    '''
    print('Performing iterative deconvolution, with a maximum of ' + str(maxbumps) + ' bumps and a ' + filt + ' filter with an fmax/gausswidth of ' + str(fmax) + '. The components  should be cut ' + str(timeshift) + ' seconds before the main arrival.')
 
    # Frequency domain
    NFFT = len(comp2) 
    frq  = np.fft.fftfreq(NFFT, dt)
    
    # Filters
    if (filt == 'gauss'):
        filterf1 = np.exp(-frq * frq / (4. * fmax * fmax))   # Gauss filter, fmax is the width of the filter and not the maximum frequency present...
        filterf = np.exp(-np.pi * np.pi * frq * frq / (fmax * fmax)) # Factor of 4 as written in Jenny's script, fortran fft uses angular frequency? 
    else:
        if (filt == 'cosine'):
            filterf = np.square(np.cos(np.pi * frq / (2. * fmax))) # fmax is the cut-off frequency
            filterf[np.abs(frq) > fmax] = 0.
        else:
            exit('Choose Gauss or Cosine filter for your iterative deconvolution')  

    # Phaseshift 
    phaseshift = [cmath.exp(1j * (frq[x] * 2. * np.pi * -(timeshift))) for x in range(len(frq))]

    # Filter components
    comp2f = np.real(np.fft.ifft(filterf * np.fft.fft(comp2, NFFT)))
    comp1f = np.real(np.fft.ifft(filterf * np.fft.fft(comp1, NFFT)))

    # Computer power in numerator
    power = sum(comp1f * comp1f)
    # Initialize values
    maxind = np.empty([maxbumps, 1], dtype=int)
    maxval = np.empty([maxbumps, 1], dtype=float)
    residual = comp2f # Initial residual is horizontal component
    fit_old  = 0 # Initial fit = 0%

    for peak in range(maxbumps):
    # correlate signals and find peak
        corr = (np.correlate(residual, comp1f, 'full'))
        corr = corr[len(comp2) - int(timeshift / dt):-int(timeshift / dt)] # Set maximum Xcorr to be at the predicted main arrival
        maxind[peak] = np.argmax(np.absolute(corr)) # Peak index
        maxval[peak] = corr[maxind[peak]] / (power) # Peak amplitude

        # build deconvolution result
        decon = np.zeros_like(comp2)
        for ind in range(peak + 1):
            decon[maxind[ind]] += maxval[ind]

        decon = np.real(np.fft.ifft(filterf * np.fft.fft(decon)))    # Filter deconvolution result

        # Reconvolve with vertical component
        conv = np.real(np.convolve(decon, comp1, 'full'))
        conv = conv[int(timeshift / dt):int(timeshift / dt) + NFFT]

        #plt.subplot(211)
        #plt.subplot(2,1,1)
        #plt.plot(np.arange(len(conv)), conv)
        #plt.subplot(2,1,2)
        #plt.plot(np.arange(len(decon)), decon)

        #plt.show()

        #print("Plot closed")

        # calculate residual and fit
        residual = comp2f - conv
        fit = 100.* (1. - sum(residual * residual) / sum(comp2f * comp2f))
        if ((fit - fit_old) < 0.01):
            break
        fit_old = fit

    print('Stopped at' + str(peak) + ' when fitting ' + str(fit) + ' %')

    return decon, fit
