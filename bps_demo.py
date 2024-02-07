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


labels = { 1000: 'AA', 4000: 'BB' }

font = {'family': 'consolas',
        'color':  'black',
        'weight': 'normal',
        'size': 9,
        }

def calc_noise(npsd, fs):    
    return 10**(npsd/10)*fs

def calc_valid_ranges (Fc, Bw):
    M = 1    
    f_max = 2*(Fc + Bw/2) 
    ranges = [(f_max, np.inf)]    
    while f_max/(M+1) >= 2*Bw:
       FsMin = (2*Fc + Bw)/(M+1);
       FsMax = (2*Fc - Bw)/M;
       ranges.append ( (FsMin,FsMax) )
       M = M + 1;
    return ranges

def build_str(ranges):
    mystr = 'ALIAS-FREE REGIONS (FS)\n'
    label_dict = {}
    for idx, entry in enumerate(ranges):
        U,L = entry
        if idx == 0:
            mystr = mystr + f'Oversampled: {U:.0f}+ [Hz]\n'
            label_dict[int(U)] = 'OS'
        else:
            mystr = mystr + f'Zone {idx}: {U:.0f} to {L:.0f} [Hz]\n'
            label_dict[int(U)] = f'Z{idx}'
            label_dict[int(L)] = f'Z{idx}'
    return mystr, label_dict


check_in_range = lambda x,r: any(lower <= x <= upper for (lower, upper) in r)
clr_dict = {True: 'limegreen', False: 'r'}

class BandpassApp():
    def __init__(self, fc=3500, bw=1000, dur=10, ns=0.1):
        self.fc = fc # RF carrier freq, Hz
        self.bw = bw # Signal bandwidth, Hz
        self.dur = dur # seconds
        self.ns = ns 
        self.fu = fc + bw/2
        self.fl = fc - bw/2
        self.line1 = None
        self.line2 = None
        self.axvline1 = None
        self.axvline2 = None
        self.base_fs = self.fu*4
        self.base_ff, self.base_psd  = self.get_psd(self.base_fs)
        
        self.ranges = calc_valid_ranges(fc, bw)
        
        self.setup()

    def get_psd(self, fs, nfft=512, noverlap=384):
        t = np.arange(0, self.dur, 1/fs)  
        y1 = chirp(t, f0=self.fl, f1=self.fc, t1=t[-1], method='hyperbolic')
        y2 = 0.50*chirp(t, f0=self.fc, f1=self.fu, t1=t[-1], method='linear')
        y = y1 + y2
        y = y + np.random.randn(np.size(y))*self.ns
        _,_,sxx = spectrogram(y,fs=fs,nperseg=nfft,noverlap=noverlap,nfft=nfft,mode='psd',return_onesided=False,window='hann',detrend=None) 
        ff = np.linspace(-fs/2, fs/2, nfft)
        Pdb = fftshift(10*np.log10(np.mean(sxx,axis=1)))
        return (ff,Pdb)
    
    def setup(self):
        zone_str, zone_labels = build_str(self.ranges)
        ui.add_body_html('<script>' + f"labels = {zone_labels}" + '</script>')
        ff, Pdb  = self.get_psd(self.base_fs)
        with ui.card().classes('bg-yellow-50'):
            with ui.column().classes('items-center'):
                ui.label('Bandpass Sampling Demo').style('font-size: 120%; font-weight: bold;')
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
                    plt.margins(x=0,y=0,tight=True)                    
                    plt.text(-7650, -33, zone_str, fontdict=font)
                    self.main_plot.fig.tight_layout()
                ui.label('Sampling Rate [Hz]:')
                ui.slider(min=1000, max=self.base_fs, step=10, value=self.base_fs).props('label-always') \
                    .on('update:model-value', lambda e: self.update_plot(e.args),throttle=0.4).classes('w-11/12').props('markers :marker-labels="labels"')
                    

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
