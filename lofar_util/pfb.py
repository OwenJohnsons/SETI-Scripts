#! usr/bin/python3
import numpy as np
import scipy
from scipy.signal import firwin, freqz, lfilter
import matplotlib.pyplot as plt
import os
import sys

"""
Polyphase filterbank for use on .rawudp files as processed at LOFAR
Based on tutorial by Danny Price
https://github.com/telegraphic/pfb_introduction
Author: Charlie Giese
"""


def load_data(path, prefix, suffix, subbands, chunks)
    #file crab_giant_pulses1_2020-03-17T17:06:12,_17901776.rawudp

    file0=open(path+prefix'0'+suffix, 'rb')
    file1=open(path+prefix'1'+suffix, 'rb')
    file2=open(path+prefix'2'+suffix, 'rb')
    file3=open(path+prefix'3'+suffix, 'rb')

    file0=np.fromfile(file0, dtype=np.int8)
    file1=np.fromfile(file1, dtype=np.int8)
    file2=np.fromfile(file2, dtype=np.int8)
    file3=np.fromfile(file3, dtype=np.int8)

    x_real=file0.reshape(-1, subbands)
    x_imag=file1.reshape(-1, subbands)
    y_real=file2.reshape(-1, subbands)
    y_imag=file3.reshape(-1, subbands)

    chunksize=int(x_real.shape[0]/chunks)
    x_complex=np.empty_like(x_real).astype(np.complex64)
    y_complex=np.empty_like(x_real).astype(np.complex64)
    for i in range(chunks):
        x_complex[i*chunksize: (i+1)*chunksize, :]=x_real[i*chunksize: (i+1)*chunksize, :] + (x_imag[i*chunksize: (i+1)*chunksize, :]*1j)
        y_complex[i*chunksize: (i+1)*chunksize, :]=y_real[i*chunksize: (i+1)*chunksize, :] + (y_imag[i*chunksize: (i+1)*chunksize, :]*1j)
    return x_complex, y_complex

def x_input(t_start, t_stop, sample_rate, x):
    #this file is (8203104 , 488) in size, want to  loop through frequency channels
    #for the time series around a pulse, 1-2 seconds in this case

    samp_start=int(t_start/sample_rate)
    samp_stop= int(t_stop/sample_rate)
    x=x[samp_start : samp_stop, :]
    return x

def y_input(t_start, t_stop, sample_rate, y):
    #this file is (8203104 , 488) in size, want to  loop through frequency channels
    #for the time series around a pulse, 1-2 seconds in this case

    samp_start=int(t_start/sample_rate)
    samp_stop= int(t_stop/sample_rate)
    y=y[samp_start : samp_stop, :]
    return y


def db(x):
    """ Convert linear value to dB value """
    return 10*np.log10(x)

def generate_win_coeffs(M, P, window_fn="hamming"):
    win_coeffs = scipy.signal.get_window(window_fn, M*P)
    sinc       = scipy.signal.firwin(M * P, cutoff=1.0/P, window="rectangular")
    win_coeffs *= sinc
    return win_coeffs

def pfb_fir_frontend(x, win_coeffs, M, P):
    W = x.shape[0] // M // P
    x_p = x.reshape((W*M, P)).T
    h_p = win_coeffs.reshape((M, P)).T
    #print(M, W, M*W -M +1)
    x_summed = np.zeros((P, M*W-M + 1)).astype(np.complex64)
    for t in range(0, M*W-M + 1):
        x_weighted = x_p[:, t:t+M] * h_p
        x_summed[:, t] = x_weighted.sum(axis=1).astype(np.complex64)
    return x_summed.T

def fft(x_p, P, axis=1):
    return np.fft.fft(x_p, P, axis=axis)

def pfb_filterbank(x, win_coeffs, M, P):
    x = x[:int(len(x)//(M*P))*M*P] # Ensure it's an integer multiple of win_coeffs
    x_fir = pfb_fir_frontend(x, win_coeffs, M, P)
    x_pfb = fft(x_fir, P)
    return x_pfb

def StokesI_forming(x, y):
    x_complex_channelised=abs(x)
    y_complex_channelised=abs(y)
    stokesI=np.empty_like(x_complex_channelised, dtype=np.int32)
    chunksize=int(stokesI.shape[0] / 256)
    for i in range(256):
        stokesI[i*chunksize:(i+1)*chunksize] =(
        np.square(y_complex_channelised[i * chunksize: (i+1) * chunksize, :].astype(np.int32)) +
        np.square(x_complex_channelised[i * chunksize: (i+1) * chunksize, :].astype(np.int32)))
    return stokesI

def pfb_spectrometer(x, n_taps, n_chan, n_int, window_fn="hamming"):

    M = n_taps
    P = n_chan
    # Generate window coefficients
    win_coeffs = generate_win_coeffs(M, P, window_fn)
    pg = np.sum(np.abs(win_coeffs)**2)
    win_coeffs /= pg**.5 # Normalize for processing gain
    # Apply frontend, take FFT, then take power (i.e. square)
    x_pfb = pfb_filterbank(x, win_coeffs, M, P)
    x_psd = np.real(x_pfb * np.conj(x_pfb))

    # Trim array so we can do time integration
    x_psd = x_psd[:np.round(x_psd.shape[0]//n_int)*n_int]

    # Integrate over time, by reshaping and summing over axis (efficient)
    x_psd = x_psd.reshape(x_psd.shape[0]//n_int, n_int, x_psd.shape[1])
    x_psd = x_psd.mean(axis=1)

    return x_psd


def spec_plot(t_start, t_stop, sample_time, x, y, n_taps, n_chan, n_int, path, prefix, suffix
             subbands, chunks):
    x, y=load_data(path, prefix, suffix, subbands, chunks)
    x_data=x_input(1., 2., 5.12*(10**(-6)), x)
    y_data=y_input(1., 2., 5.12*(10**(-6)), y)

    x_psd_trial=pfb_spectrometer(x_data[:, 0], n_taps, n_chan, n_int, window_fn="hamming")
    y_psd_trial=pfb_spectrometer(y_data[:, 0], n_taps, n_chan, n_int, window_fn="hamming")

    x_array=np.zeros(shape=(x_psd_trial.shape[0], x_data.shape[1]*n_chan))
    y_array=np.zeros(shape=(y_psd_trial.shape[0], y_data.shape[1]*n_chan))
    print(np.shape(x_array))
    for i in range(x_data.shape[1]):
        x=x_data[:, i]
        y=y_data[:, i]
        x_psd_i=pfb_spectrometer(x, n_taps, n_chan, n_int, window_fn="hamming")
        y_psd_i=pfb_spectrometer(y, n_taps, n_chan, n_int, window_fn="hamming")
        x_array[:, i*n_chan: (i+1)*n_chan]=x_psd_i
        y_array[:, i*n_chan: (i+1)*n_chan]=y_psd_i

    stokesI=StokesI_forming(x_array, y_array)

    mean_list=[]
    chan_list=[]
    for i in range(stokesI.shape[1]):
        chan=stokesI[:, i]
        temp=np.mean(chan, axis=0)
        mean_list.append(temp)
        chan_list.append(i)
    means=np.array(mean_list)
    indices=np.argwhere(means >= 6000)
    stokesI[:, indices]=1

    indices=[2825, 3904, 4922, 5503, 7419, 7404, 7408, 7409, 7410, 7418, 7403, 7402, 7447, 7443,
    7380, 7381, 7382, 7383, 7384, 7385, 7393, 7394, 7395, 7556, 7538, 7680, 7695]
    stokesI[:, indices]=1
    stokesI[:, 7350: ]=1

    plt.imshow(db(stokesI), cmap='viridis', aspect='auto')
    plt.colorbar()
    plt.xlabel("Channel")
    plt.ylabel("Time")
    plt.savefig(n_chan+"_pfb_StokesI.png")
    stokesI.astype(np.float32).tofile(n_chan+"_pfb_stokesI.fil")
    np.save(n_chan+'_pfb_stokesI.npy', stokesI.astype(np.int32))

output=spec_plot(t_start, t_stop, sample_time, n_taps, n_chan, n_int, path, prefix,
                     suffix, subbands, chunks):
"""
t_start=time in seconds of where to start pfb input
t_stop=time in seconds of where to stop pfb input
sample_time=sample time of data in seconds
n_taps=number of taps of pfb
n_chan=number of channels for channelisation
n_int=duration of time integration in seconds
path=path to directory of data
prefix=file name before port
suffix=file name after port
subbands=number of subbands in raw data
chunks=number of chunks to split data into to avoid memory issues
"""
