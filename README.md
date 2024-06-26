# bps_demo
A small app to visually demonstrate the concept of bandpass sampling and alias regions. 

[Click here for live version](http://www.blankadventure.com/bps)


![image](https://github.com/BlankAdventure/bps_demo/assets/24900496/f6b9cf0b-fdc7-4383-8484-b2dc91844202)



The default settings (shown above) depict the two-sided power spectrum for a real signal with an "RF" carrier of 3500 Hz and a lower and upper sideband occupying a bandwidth of 1000 Hz. In this case the maximum frequency is 4000 Hz and we start with a sampling rate of 12000 Hz, yielding a Nyquist frequency of 6000 Hz. Thus we are initially oversampling the signal.

The sampling rate slider can be adjusted to see the impact of different sampling rates on the reconstructed signal. The black line depics the original (oversampled) signal and the green line depicts the modified signal spectrum. The two dashed vertical lines indicate the location of the Nyquist frequency (i.e., +- Fs/2). The upper left legend lists the alias-free sampling rates, which are also indicated graphically by the green strips below the slider (the red strips indicate the aliased regions).

In this example the maximum frequency is 4000 Hz therefore anything above 8000 Hz would be considered oversampling. In the usual understanding of Nyquist sampling, one would consider that we must sample at 8000 Hz+ if we want to properly reconstruct the signal. Indeed, by positioning the slider below 8000 Hz, such aliasing artifacts appear as expected:
![image](https://github.com/BlankAdventure/bps_demo/assets/24900496/f1bc1388-a64c-4e3b-ac9f-c6c63ccda69e)



However, for the case of bandpass signals, it turns out there are certain under-sampling rates that preserve the signal. Examples from the 3 valid zones are depicted below:


![image](https://github.com/BlankAdventure/bps_demo/assets/24900496/3c7d5447-1ff1-44a3-b5ab-44b61997412c)
![image](https://github.com/BlankAdventure/bps_demo/assets/24900496/119aed1e-9c5e-499d-9927-9e2be361213d)
![image](https://github.com/BlankAdventure/bps_demo/assets/24900496/d1624639-bde1-469b-89a8-22f67d2c4ddb)


Whether or not these Nyquest zones are useful is another question. Notice the spectra may become inverted and they don't necessarily land at DC/0 Hz. 

Two additional slides are provided to see how the carrier frequency and bandwidth affect the Nyquist zones.
