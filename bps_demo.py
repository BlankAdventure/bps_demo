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
from nicegui import ui, run
import asyncio


async def worker(name, queue):
    while True:
        res = await queue.get()
        res()
        #await make_async(res)()
        queue.task_done()

def make_async(function_to_decorate):
    async def async_wrap(*args,**kwargs):
        return await asyncio.to_thread(function_to_decorate, *args)
    return async_wrap

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


def noise_power(npsd, fs):    
    return np.sqrt(10**(npsd/10)*fs)

# Convert list of valid ranges into a readable string
def build_str(ranges):
    mystr = '--ALIAS-FREE REGIONS (FS)--\n'
    for idx, entry in enumerate(ranges):
        U,L = entry
        if idx == 0:
            mystr = mystr + f' Oversampled: {U:.0f}+ [Hz]\n'
        else:
            mystr = mystr + f' Zone {idx}: {U:.0f} to {L:.0f} [Hz]\n'
    return mystr

def get_psd(fs, fc, fl, fu, dur, npsd, nfft=512, noverlap=384):
    t = np.arange(0, dur, 1/fs)  
    y1 = chirp(t, f0=fl, f1=fc, t1=t[-1], method='hyperbolic')
    y2 = 0.50*chirp(t, f0=fc, f1=fu, t1=t[-1], method='linear')
    y = y1 + y2
    y = y + np.random.randn(np.size(y))*noise_power(npsd, fs)
    _,_,sxx = spectrogram(y,fs=fs,nperseg=nfft,noverlap=noverlap,nfft=nfft,mode='psd',return_onesided=False,window='hann',detrend=None) 
    ff = np.linspace(-fs/2, fs/2, nfft)
    Pdb = fftshift(10*np.log10(np.mean(sxx,axis=1)))
    return (ff, Pdb)


class psdClass():
    def __init__(self, fc=3500, bw=1000, dur=5, npsd=-60):
        self.dur = dur 
        self.npsd = npsd         
        self.update_ref_psd(fc, bw)
        self.update_test_psd(self.base_fs)

    def calc_base_vals(self, new_fc, new_bw):
        self.fc = new_fc
        self.bw = new_bw
        self.fu = self.fc + self.bw/2
        self.fl = self.fc - self.bw/2
        self.min_fs = self.bw*2
        self.base_fs = self.fu*3
        self.ranges = calc_valid_ranges(self.base_fs, self.fc, self.bw)        

    def update_ref_psd(self, new_fc, new_bw):
        self.calc_base_vals(new_fc, new_bw)
        ff, pdb = get_psd(self.base_fs, self.fc, self.fl, self.fu, self.dur, self.npsd)
        self.ref_ff = ff
        self.ref_psd = pdb
        return ff, pdb
        
    def update_test_psd(self, new_fs):
        ff, pdb = get_psd( new_fs, self.fc, self.fl, self.fu, self.dur, self.npsd)
        self.test_ff = ff
        self.test_psd = pdb
        return ff, pdb


check_in_range = lambda x,r: any(lower <= x <= upper for (lower, upper) in r)
clr_dict = {True: 'limegreen', False: 'r'}


def test_func(x):
    print(f'this was called with {x}')
    


@ui.page('/main')
def main():
    
    psd = psdClass()
    queue = asyncio.LifoQueue(maxsize=1)
    
    async def ur(fc, bw):
        try: 
            queue.put_nowait( lambda: update_ref_psd(fc, bw) )
        except:
            print('q full')

    async def ut(fs):
        try: 
            queue.put_nowait( lambda: update_test_psd(fs) )
        except:
            print('q full')


    def draw_zonebar():
        with zone_bar:
            ofs = 0
            w = 0
            for L,U in reversed(psd.ranges):
                ofs = ofs + w
                w = 100*(U-L)/(psd.base_fs - psd.min_fs)
                p = 100*(L-psd.min_fs)/(psd.base_fs - psd.min_fs)
                if w == 0:
                    w = 0.1
                ui.element('div').classes('bg-green-400').style(f'position: relative; left: {p-ofs}%; height: 15px; width: {w}%;')    
    
    def update_ref_psd(new_fc, new_bw):
        psd.update_ref_psd(new_fc, new_bw)
        title.set_text(f'(fc={psd.fc:.0f} Hz | BW={psd.bw:.0f} Hz | fmax={psd.fu:.0f} Hz)')
        regions.set_text(build_str(psd.ranges))
        samp_slider.props(f'markers :min={psd.min_fs}')
        samp_slider.props(f'markers :max={psd.base_fs}')        
        if samp_slider.value > psd.base_fs:
            samp_slider.value = psd.base_fs        
        zone_bar.clear()
        draw_zonebar()        

        with main_plot:
            line1.set_xdata(psd.ref_ff)
            line1.set_ydata(psd.ref_psd)
            plt.xlim(-psd.base_fs/2,psd.base_fs/2)        
        
        update_test_psd(samp_slider.value)

    def update_test_psd(new_fs):
        psd.update_test_psd(new_fs)
        with main_plot:
            line2.set_xdata(psd.test_ff)
            line2.set_ydata(psd.test_psd)
            line2.set( color=clr_dict[check_in_range(new_fs,psd.ranges)] )
            axvline1.set_data([new_fs/2, new_fs/2], [0, 1])
            axvline2.set_data([-new_fs/2, -new_fs/2], [0, 1])

    
    
    with ui.card().classes('bg-yellow-50'):
        with ui.column().classes('gap-0'): 
            ui.label('Bandpass Sampling Demo').style('font-size: 120%; font-weight: bold;').classes('w-full text-center')
            title = ui.label(f'(fc={psd.fc:.0f} Hz | BW={psd.bw:.0f} Hz | fmax={psd.fu:.0f} Hz)').style('font-size: 90%; font-weight: normal;').classes('w-full text-center')
            with ui.pyplot(figsize=(9, 5)) as main_plot:
                main_plot.fig.patch.set_alpha(0)                    
                line1, = plt.plot(psd.ref_ff,psd.ref_psd, color='k')
                line2, = plt.plot(psd.ref_ff,psd.ref_psd, color=clr_dict[check_in_range(psd.fu*2,psd.ranges)])          
                axvline1 = plt.axvline( psd.base_fs/2, linestyle='--', color='grey', alpha=0.60)
                axvline2 = plt.axvline(-psd.base_fs/2, linestyle='--', color='grey', alpha=0.60)
                ax = plt.gca()
                ax.patch.set_alpha(0)
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.spines['left'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
                plt.xlim(-psd.base_fs/2,psd.base_fs/2)
                plt.ylim(-65, -25)                    
                plt.xlabel('Freq [Hz]')
                plt.ylabel('PSD [dBW/Hz]')
                plt.margins(x=0,y=0,tight=True)                    
                main_plot.fig.tight_layout()
            
            # Plot annotation of aliasinsg zones
            regions = ui.label(build_str((psd.ranges))).classes().style("position: absolute; top: 13%; left: 12%; white-space: pre; font-size: 85%; font-weight: 500; border: 1px solid; padding: 3px; background-color: rgba(249, 250, 251, 0.88);")
            
            # Setup the slider
            ui.label('Sampling Rate [Hz]:').classes('text-left italic')
            samp_slider = ui.slider(min=psd.min_fs, max=psd.base_fs, step=5, value=psd.base_fs).props('label-always') \
                .on('update:model-value', lambda e: ut(e.args),throttle=1,leading_events=False).classes('w-full').props()
            
            # Setup the indicator bar
            with ui.row().classes('w-full gap-0 bg-red-300').style('position: relative; top: -10px;') as zone_bar:
                draw_zonebar()
        
            with ui.row().classes('w-full items-center justify-left'):
                ui.label('Carrier Freq [Hz]:').classes('italic')
                ui.slider(min=2500,max=4500,step=50,value=psd.fc).style('width: 35%;').props('label-always switch-label-side').on('update:model-value', lambda e: ur(e.args,psd.bw),throttle=1,leading_events=False)
                ui.label('Bandwidth [Hz]:').classes('italic')
                ui.slider(min=500,max=1500,step=50,value=psd.bw).style('width: 35%;').props('label-always switch-label-side').on('update:model-value', lambda e: ur(psd.fc,e.args),throttle=1,leading_events=False)


    task = asyncio.create_task(worker('name', queue))        
    
    


ui.run(port=5000, on_air=False,title='Bandpass Sampling Demo',host='0.0.0.0')

