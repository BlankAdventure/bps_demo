# -*- coding: utf-8 -*-
"""
Demo code of creating a real-valued whtie noise signal with a specified
power spectral density
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from scipy.fft import fftshift

# Generate a real-value vector of AWGN samples with specified PSD
# fs -> sample rate [Hz]
# npsd -> desired noise power spectral density [dBW/Hz]
# dur -> signal duration
def get_noise(fs, npsd, dur=10):
    npwr = np.sqrt(10**(npsd/10)*fs)
    nsamps = dur*fs
    return np.random.randn(nsamps)*npwr

# Generate a spectrogram 
# sig -> signal
# fs -> sample rate [Hz]
# nfft -> # points to use for fft
# noverlap -> # samples to overlap (sliding window is used)
def get_psd(sig, fs, nfft=512, noverlap=384):
    _,_,sxx = spectrogram(sig,fs=fs,nperseg=nfft,noverlap=noverlap,nfft=nfft,mode='psd',return_onesided=False,window='rect',detrend=None) 
    ff = np.linspace(-fs/2, fs/2, nfft)
    P = fftshift(np.mean(sxx,axis=1))
    return (ff,P)


# Example values
npsd = -50 # dBW/Hz
fs = 15000 # Hz
ns = get_noise(10000,-50)
ff, P = get_psd(ns, 10000)


# Observe in the plot our PSD if flat and falls at -50 [dbW/Hz] as requested
plt.plot(ff,10*np.log10(P))
print ( 10*np.log10(np.mean(P)) )