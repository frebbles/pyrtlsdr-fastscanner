# Python RTLSDR Scanner and fm listen.
#
# Author: Farran Rebbeck
# Date  : 20161214

# imports
# NOTE: Requires pyrtlsdr, numpy

from rtlsdr import RtlSdr
import numpy as np
sdr = RtlSdr()
import scipy.signal as signal
import sounddevice as sd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Config variables (bracketing for bands etc)
samplecount = 1024
freqsteps = 12500
sdr_scan_sample_rate =   2000000
sdr_scan_center_freq = 145000000
tune_offset = 20000
tune_sample_rate = 300000
# Specific tuning offset...... possibly something weird
rtl_tune_offset = 4000 

# Can be 'AUTO' or in dB
sdr.gain = 4

scanning = 0

# Main program loop
while True:

  if (scanning == 0):
    # Set the center freq for scanning
    if (sdr.center_freq != sdr_scan_center_freq):
      sdr.center_freq = sdr_scan_center_freq

    # Set the sampling rate for scanning
    if (sdr.sample_rate != sdr_scan_sample_rate):
      sdr.sample_rate = sdr_scan_sample_rate
    
    scanning = 1

  # Take some samples from the RTLSDR
  samples = sdr.read_samples(samplecount)
  
  # FFT for the sample set
  fft = np.fft.fft(samples)
  freqs = np.absolute(np.fft.fft(samples))
  freqs = freqs[1:-1]

  # Shift freqs to align them with the next bunch of calculations
  freqs = np.fft.fftshift(freqs)

  # Convert the FFT to dB
  freqs = 20.0*np.log10(freqs)

  # Figure out start/finish frequencies for numerical correspondence
  step = sdr.sample_rate / samplecount
  start = sdr.center_freq - (sdr.sample_rate/2)
  curfreq=start
  fcount = len(freqs)

  # Variables for the high count algorithm
  highfreq = 0;
  highfreqv = -100;
  hfind = 0;

  # Loop through frequency buckets find the highest signal.
  for idx, f in enumerate(freqs):
    curfreq += step
    if f > highfreqv:
      highfreq = curfreq
      highfreqv = f
      hfind = idx

  # For the highest frequency sensed, we hone in on a multiple of {freqsteps} by using a simple rounding calculation
  highfreq = int(highfreq / freqsteps) +1
  highfreq = highfreq * freqsteps

  print 'High freq = %d, %d dB (%d/%d)' % (highfreq, highfreqv, hfind, fcount)

  # Determine if the peak is worth investigating and tuning to
  if (highfreqv > 40):
  
    print "Retuning"
    # Tune to the frequency
    sdr.sample_rate = tune_sample_rate
    sdr.center_freq = highfreq + tune_offset

    scanning = 0

    print "...Retuned"
    print "Reading Samples"

    F_offset = -tune_offset - rtl_tune_offset
    Fs = tune_sample_rate

    # Capture samples
    psamples = sdr.read_samples(512*1024)
    x1 = np.array(psamples).astype("complex64")
    
    print "...Read " + str(len(x1)) + " samples"

    print "Mixdown"

    # Mix data down
    fc1 = np.exp(-1.0j * 2.0 * np.pi * F_offset/Fs * np.arange(len(x1)))
    x2 = x1 * fc1
    
    # NFM Signal processing
    f_bw = 12500
    n_taps = 2

    print "Filtering"

    lpf = signal.remez(n_taps, [0, f_bw, f_bw+(Fs/2-f_bw)/4, Fs/2], [1,0], Hz=Fs)
    x3 = signal.lfilter(lpf, 1.0, x2)

    dec_rate = int(Fs / f_bw)

    print "Decimating"
    
    x4 = x3[0::dec_rate]
    
    Fs_y = Fs/dec_rate

    print "FM Demod"
 
    y5 = x4[1:] * np.conj(x4[:-1])
    x5 = np.angle(y5)

    print "FM Deemp"    

    # De emphasis
    d = Fs_y * 75e-6
    x = np.exp(-1/d)
    b = [1-x]
    a = [1,-x]
    x6 = signal.lfilter(b,a,x5)

    print "Audio Dec"

    audio_freq = 12500
    dec_audio = int(Fs_y/audio_freq)
    Fs_audio = Fs_y / dec_audio

    x7 = signal.decimate(x6, dec_audio)

    x7 *= 100000 / np.max(np.abs(x7))
    
    print "Audio out"

    sd.play(x7.astype("int16"), audio_freq)



