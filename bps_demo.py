# -*- coding: utf-8 -*-
"""

Bandpass sampling demo tool
"""

# The below two lines were required to deal with a weird issue described here:
#https://stackoverflow.com/questions/15457786/ctrl-c-crashes-python-after-importing-scipy-stats

import os
os.environ['FOR_DISABLE_CONSOLE_CTRL_HANDLER'] = '1'


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import chirp, spectrogram
from scipy.fft import fftshift
from nicegui import ui

font = {'family': 'consolas',
        'color':  'black',
        'weight': 'normal',
        'size': 9,
        }

def noise_power(npsd, fs):    
    return np.sqrt(10**(npsd/10)*fs)

def calc_valid_ranges (Fs, Fc, Bw):
    M = 1    
    f_max = 2*(Fc + Bw/2) 
    ranges = [(f_max, Fs)]    
    while f_max/(M+1) >= 2*Bw:
       FsMin = (2*Fc + Bw)/(M+1);
       FsMax = (2*Fc - Bw)/M;
       ranges.append ( (FsMin,FsMax) )
       M = M + 1;
    return ranges

def build_str(ranges):
    mystr = '--ALIAS-FREE REGIONS (FS)--\n'
    label_dict = {}
    for idx, entry in enumerate(ranges):
        U,L = entry
        if idx == 0:
            mystr = mystr + f' Oversampled: {U:.0f}+ [Hz]\n'
            label_dict[int(U)] = 'OS'
        else:
            mystr = mystr + f' Zone {idx}: {U:.0f} to {L:.0f} [Hz]\n'
            label_dict[int(U)] = f'Z{idx}'
            label_dict[int(L)] = f'Z{idx}'
    return mystr, label_dict


check_in_range = lambda x,r: any(lower <= x <= upper for (lower, upper) in r)
clr_dict = {True: 'limegreen', False: 'r'}

class BandpassApp():
    def __init__(self, fc=3500, bw=1000, dur=10, npsd=-60):
        self.fc = fc # RF carrier freq, Hz
        self.bw = bw # Signal bandwidth, Hz
        self.dur = dur # seconds
        self.npsd = npsd 
        self.fu = fc + bw/2
        self.fl = fc - bw/2
        self.line1 = None
        self.line2 = None
        self.axvline1 = None
        self.axvline2 = None
        self.min_fs = bw*2
        self.base_fs = self.fu*3
        self.base_ff, self.base_psd  = self.get_psd(self.base_fs)       
        self.ranges = calc_valid_ranges(self.base_fs, fc, bw)        
        print(self.ranges)

        self.setup()
    def get_psd(self, fs, nfft=512, noverlap=384):
        t = np.arange(0, self.dur, 1/fs)  
        y1 = chirp(t, f0=self.fl, f1=self.fc, t1=t[-1], method='hyperbolic')
        y2 = 0.50*chirp(t, f0=self.fc, f1=self.fu, t1=t[-1], method='linear')
        y = y1 + y2
        y = y + np.random.randn(np.size(y))*noise_power(self.npsd, fs)
        _,_,sxx = spectrogram(y,fs=fs,nperseg=nfft,noverlap=noverlap,nfft=nfft,mode='psd',return_onesided=False,window='hann',detrend=None) 
        ff = np.linspace(-fs/2, fs/2, nfft)
        Pdb = fftshift(10*np.log10(np.mean(sxx,axis=1)))
        return (ff, Pdb)
    
    def setup(self):
        zone_str, zone_labels = build_str(self.ranges)
        ff, Pdb  = self.get_psd(self.base_fs)
        with ui.card().classes('bg-yellow-50'):
            with ui.column().classes('gap-0'): #items-center
                ui.label('Bandpass Sampling Demo').style('font-size: 120%; font-weight: bold;').classes('w-full text-center')
                ui.label(f'(fc={self.fc:.0f} Hz | BW={self.bw:.0f} Hz | fmax={self.fu:.0f} Hz)').style('font-size: 90%; font-weight: normal;').classes('w-full text-center')
                self.main_plot = ui.pyplot(figsize=(9, 5))
                with self.main_plot:
                    self.main_plot.fig.patch.set_alpha(0)                    
                    self.line1, = plt.plot(self.base_ff,self.base_psd, color='k')
                    self.line2, = plt.plot(ff,Pdb,color=clr_dict[check_in_range(self.fu*2,self.ranges)])          
                    self.axvline1 = plt.axvline( self.base_fs/2, linestyle='--', color='grey', alpha=0.60)
                    self.axvline2 = plt.axvline(-self.base_fs/2, linestyle='--', color='grey', alpha=0.60)
                    ax = plt.gca()
                    ax.patch.set_alpha(0)
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    ax.spines['left'].set_visible(True)
                    ax.spines['bottom'].set_visible(True)
                    plt.xlim(-self.base_fs/2,self.base_fs/2)
                    plt.ylim(-65, -25)                    
                    plt.xlabel('Freq [Hz]')
                    plt.ylabel('PSD [dBW/Hz]')
                    plt.margins(x=0,y=0,tight=True)                    
                    self.main_plot.fig.tight_layout()
                
                # Plot annotation of aliasinsg zones
                ui.label(zone_str).classes('bg-gray-50').style("position: absolute; top: 13%; left: 12%; white-space: pre; font-size: 85%; font-weight: 500; border: 1px solid; padding: 3px;")

                # Setup the slider
                ui.label('Sampling Rate [Hz]:').classes('text-left italic')
                ui.slider(min=self.min_fs, max=self.base_fs, step=5, value=self.base_fs).props('label-always') \
                    .on('update:model-value', lambda e: self.update_plot(e.args),throttle=0.4).classes('w-full').props()
                
                # Setup the indicator bar
                with ui.row().classes('w-full gap-0 bg-red-300').style('position: relative; top: -10px;'):
                    ofs = 0
                    w = 0
                    for L,U in reversed(self.ranges):
                        ofs = ofs + w
                        w = 100*(U-L)/(self.base_fs - self.min_fs)
                        p = 100*(L-self.min_fs)/(self.base_fs - self.min_fs)
                        if w == 0:
                            w = 0.1
                        ui.element('div').classes('bg-green-400').style(f'position: relative; left: {p-ofs}%; height: 15px; width: {w}%;')


    def update_plot(self, fs):
        ff, Pdb  = self.get_psd(fs)        
        with self.main_plot:
            self.line2.set_xdata(ff)
            self.line2.set_ydata(Pdb)
            self.line2.set( color=clr_dict[check_in_range(fs,self.ranges)] )
            self.axvline1.set_data([fs/2, fs/2], [0, 1])
            self.axvline2.set_data([-fs/2, -fs/2], [0, 1])

    def run(self,port=5000):
        ui.run(port=port)



app = BandpassApp()
app.run()            

