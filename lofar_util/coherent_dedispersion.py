
import os
import numpy as np
import mmap

"""
Slightly Modified version of Griffin Foster's dedispersion.py script for coherent dedispersion.
https://github.com/griffinfoster/frb180301-analysis
Modified by Charlie Giese
"""



subbands=488
f_top=197.55859375e6
f_off=-0.1953125e6
sample_rate=1/5.12e-6
dm=56.77
try:
	import os
	import pickle
	import pyfftw

	from lsl.common.paths import data as dataPath

	# Enable the PyFFTW cache
	if not pyfftw.interfaces.cache.is_enabled():
		pyfftw.interfaces.cache.enable()
		pyfftw.interfaces.cache.set_keepalive_time(60)

	# Read in the wisdom (if it exists)
	wisdomFilename = os.path.join(dataPath, 'pyfftw-wisdom.pkl')
	if os.path.exists(wisdomFilename):
		fh = open(wisdomFilename, 'r')
		wisdom = pickle.load(fh)
		fh.close()

		pyfftw.import_wisdom(wisdom)
		useWisdom = True
	else:
		useWisdom = False

	usePyFFTW = True

except ImportError:
	usePyFFTW = False
	useWisdom = False


__version__ = '0.5'
__revision__ = '$Rev$'
__all__ = ['delay', 'incoherent', 'getCoherentSampleSize', 'coherent', '__version__', '__revision__', '__all__']


# Dispersion constant in MHz^2 s / pc cm^-3
_D = 4.148808e3


# Coherent dedispersion N and chirp cache
_coherentCache = {}

def delay(freq, dm):
	"""
	Calculate the relative delay due to dispersion over a given frequency
	range in Hz for a particular dispersion measure in pc cm^-3.  Return
	the dispersive delay in seconds.

	.. versionchanged:: 1.1.1
		If only a single frequency is provided, the returned delay is
		relative to infinite frequency.
	"""

	# Validate in input frequencies
	## Right Type?
	try:
		freq.size
	except AttributeError:
		freq = np.array(freq, ndmin=1)
	## Right size?
	singleFreq = False
	if freq.size == 1:
		singleFreq = True
		freq = np.append(freq, np.inf)

	# Delay in s
	tDelay = dm*_D*((1e6/freq)**2 - (1e6/freq.max())**2)

	# Cleanup
	if singleFreq:
		tDelay = tDelay[0]

	return tDelay


def getCoherentSampleSize(centralFreq, sampleRate, dm):
	"""
	Estimate the number of samples needed to successfully apply coherent
	dedispersion to a data stream.
	"""

	# Roughly estimate the number of points we need to look at to do the dedispersion
	# correctly.  Based on the the relative dispersion delay between the high and low
	# ends of an observational band.
	F0 = centralFreq
	BW = sampleRate

	delayBand = dm*_D *((1e6/(F0-BW/2.0))**2 - (1e6/(F0+BW/2.0))**2)	# Dispersion delay across the band
	samples = delayBand*BW									# Conversion to samples
	samples = 2**(np.ceil(np.log(samples)/np.log(2)))		# Conversion to next largest power of 2
	samples *= 2											# Correction for the 'wings' of the convolution
	return int(samples)


def __taperFunction(freq):
	"""
	Taper function based Equation (1) of "Pulsar Coherent De-dispersion
	Experiment at Urumqi Observatory" CJA&A, 2006, S2, 53.
	"""

	freqMHz = freq / 1e6
	fMHz0 = freqMHz.mean()
	fMHz1 = freqMHz - fMHz0
	BW = fMHz1.max() - fMHz1.min()

	taper = 1.0 / np.sqrt( 1.0 + ( np.abs(fMHz1) / (0.47*BW) )**80 )

	return taper


def __chirpFunction(freq, dm, taper=False):
	"""
	Chip function for coherent dedispersion for a given set of frequencies (in Hz).
	Based on Equation (6) of "Pulsar Observations II -- Coherent Dedispersion,
	Polarimetry, and Timing" By Stairs, I. H.
	"""

	freqMHz = freq / 1e6
	fMHz0 = freqMHz.mean()
	fMHz1 = freqMHz - fMHz0
	BW = fMHz1.max() - fMHz1.min()

	chirp = np.exp(-2j*np.pi*_D*1e6 / (fMHz0**2*(fMHz0 + fMHz1)) * dm*fMHz1**2)
	if taper:
		chirp *= __taperFunction(freq)

	return chirp


def coherent(t, timeseries, centralFreq, sampleRate, dm, taper=False, previousTime=None, previousData=None, nextTime=None, nextData=None, enableCaching=True):
	"""
	Simple coherent dedispersion of complex-valued time-series data at a given central
	frequency and sample rate.  A tapering function can also be applied to the chirp of
	the form:

		:math:`\\sqrt{1 + \\left(\\frac{\\Delta f_{MHz}}{0.47 \\times \\mbox{BW}}\\right)^{80}}`,

	where :math:`\\Delta f_{MHz}` is the frequency difference in MHz from the band
	center and BW is the bandwidth in MHz.

	.. note::
		At the large fractional bandwidths of LWA, the window size needed for coherent
		dedispersion can be prohibitive.  For example, at 74 MHz with 19.6 MS/s and a
		DM or 10 pc / cm^3 this function uses a window size of about 268 million points.

	.. versionchanged:: 1.0.1
		Added support for using PyFFTW instead of NumPy for the FFTs and iFFTs.
		Added a cache for storing the chrip function between subsequent calls

	.. versionchanged:: 0.6.4
		Added support for keeping track of time through the dedispersion process.
	"""

	# Caching for N and the chirp function
	try:
		pair = (centralFreq*1.0, sampleRate*1.0, dm*1.0, taper, timeseries.dtype)

		# Get an idea of how many samples we need to do the dedispersion correctly
		# Compute the chirp function
		N, chirp = _coherentCache[pair]
	except KeyError:
		# Get an idea of how many samples we need to do the dedispersion correctly
		N = getCoherentSampleSize(centralFreq, sampleRate, dm)

		# Compute the chirp function
		freq = np.fft.fftfreq(N, d=1.0/sampleRate) + centralFreq
		chirp = __chirpFunction(freq, dm, taper=taper)
		chirp = chirp.astype(timeseries.dtype)

		# Update the cache
		if enableCaching:
			_coherentCache[pair] = (N, chirp)

	# Figure out the output array size
	nSets = len(timeseries) // N
	outT = np.zeros(timeseries.size, dtype=t.dtype)
	outD = np.zeros(timeseries.size, dtype=timeseries.dtype)

	if nSets == 0:
		RuntimeWarning("Too few data samples for proper dedispersion")

	if usePyFFTW:
		if timeseries.dtype not in (np.complex64, np.complex128):
			raise RuntimeError("Unsupported data type for timeseries: %s" % str(timeseries.dtype))

		dd = timeseries.dtype
		di1 = np.empty(N, dtype=dd)
		do1 = np.empty(N, dtype=dd)
		do2 = np.empty(N, dtype=dd)

		forwardPlan = pyfftw.FFTW(di1, do1, direction='FFTW_FORWARD', flags=('FFTW_ESTIMATE', 'FFTW_UNALIGNED'))
		backwardPlan = pyfftw.FFTW(do1, do2, direction='FFTW_BACKWARD', flags=('FFTW_ESTIMATE', 'FFTW_UNALIGNED'))

	# Go!
	for i in range(2*nSets+1):
		start = int(i*N/2 - N/4)
		stop = start + N

		if start < 0:
			timeIn = np.zeros(N, dtype=t.dtype)
			dataIn = np.zeros(N, dtype=timeseries.dtype)

			if previousData is not None:
				try:
					timeIn[:-start] = previousTime[start:]
					dataIn[:-start] = previousData[start:]
				except ValueError:
					raise RuntimeError("Too few data samples for proper start buffering")

			timeIn[-start:N] = t[0:N+start]
			dataIn[-start:N] = timeseries[0:N+start]

		elif stop > timeseries.size:
			timeIn = np.zeros(N, dtype=t.dtype)
			dataIn = np.zeros(N, dtype=timeseries.dtype)

			if nextData is not None:
				ns = nextData.size
				df = dataIn.size - (timeseries.size-start)

				try:
					timeIn[timeseries.size-start:] = nextTime[:(dataIn.size-(timeseries.size-start))]
					dataIn[timeseries.size-start:] = nextData[:(dataIn.size-(timeseries.size-start))]
				except ValueError:
					raise RuntimeError("Too few data samples for proper end buffering")

			if start < timeseries.size:
				timeIn[0:(timeseries.size-start)] = t[start:]
				dataIn[0:(timeseries.size-start)] = timeseries[start:]

		else:
			timeIn = t[start:stop]
			dataIn = timeseries[start:stop]

		timeOut = timeIn
		if usePyFFTW:
			forwardPlan(dataIn, do1)
			do1 *= chirp
			backwardPlan(do1, do2)
			dataOut = do2

		else:
			dataOut = np.fft.fft( dataIn )
			dataOut *= chirp
			dataOut = np.fft.ifft( dataOut )

		# Get the output data ranges
		outStart  = i*N//2
		outStop   = outStart + N//2
		dataStart = N//4
		dataStop  = dataStart + N//2

		# Make sure we don't fall off the end of the array
		if outStop >= outD.size:
			diff = outStop - outD.size
			outStop  -= diff
			dataStop -= diff

		if i == 0 and previousData is None:
			continue
		if i == (2*nSets) and nextData is None:
			continue

		outT[outStart:outStop] = timeOut[dataStart:dataStop]
		outD[outStart:outStop] = dataOut[dataStart:dataStop]

	return outT, outD

def times(data, sample_rate):
    times=np.arange(0, data.shape[0], 1)/sample_rate
    return times

def freq_array(nbands, f_top, f_off): #in MHz
    f_min=f_top+(nbands*f_off)
    freqs=np.arange(f_top, f_min, f_off)
    return freqs

def stokesI(frequency, sample_rate, dm, nbands):
    x=x_complex.astype(np.complex64)
    y=y_complex.astype(np.complex64)
    xdd=np.zeros_like(x)
    ydd=np.zeros_like(y)
    t=times(x, sample_rate)

    temp= coherent(t, x, frequency, sample_rate, dm, taper=False, previousTime=None,
    previousData=None, nextTime=None, nextData=None, enableCaching=True)
    xdd=temp[0]

    temp= coherent(t, y, frequency, sample_rate, dm, taper=False, previousTime=None,
    previousData=None, nextTime=None, nextData=None, enableCaching=True)
    ydd=temp[0]


    x_coherently_dd=abs(xdd)
    y_coherently_dd=abs(ydd)
    stokesI=np.empty_like(x_coherently_dd, dtype=np.int32)
    stokesI=(np.square(x_coherently_dd.astype(np.int32)) + np.square(y_coherently_dd.astype(np.int32)))
    return stokesI

N=getCoherentSampleSize(np.min(freq_array(subbands, f_top, f_off)), sample_rate, dm)
print(N)
path="/mnt/ucc4_data2/data/David/crab_extracted_pulses/"

with open(path+"/crab_giant_pulses0_2020-03-17T19:25:13,_17901936.rawudp", "rb") as f:
    map0 = np.memmap(f, mode="r", dtype=np.int8)
Nsets=int(map0.shape[0]/N)
del map0
f.close()
band=0
for i in range(Nsets):
    with open(path+"/crab_giant_pulses0_2020-03-17T19:25:13,_17901936.rawudp", "rb") as f:
        map0 = np.memmap(f, mode="r", dtype=np.int8)
        x_real=np.array(map0[i*N:(i+1)*N]).astype(np.int8)
        del map0
    with open(path+"/crab_giant_pulses1_2020-03-17T19:25:13,_17901936.rawudp", "rb") as f:
        map1 = np.memmap(f, mode="r", dtype=np.int8)
        x_imag=np.array(map1[i*N:(i+1)*N]).astype(np.int8)
        del map1
    with open(path+"/crab_giant_pulses2_2020-03-17T19:25:13,_17901936.rawudp", "rb") as f:
        map2 = np.memmap(f, mode="r", dtype=np.int8)
        y_real=np.array(map2[i*N:(i+1)*N]).astype(np.int8)
        del map2
    with open(path+"/crab_giant_pulses3_2020-03-17T19:25:13,_17901936.rawudp", "rb") as f:
        map3 = np.memmap(f, mode="r", dtype=np.int8)
        y_imag=np.array(map3[i*N:(i+1)*N]).astype(np.int8)
        del map3
    f.close()
    x_complex=np.empty_like(x_real, dtype=np.complex64)
    y_complex=np.empty_like(x_real, dtype=np.complex64)
    frequency=freq_array(subbands, f_top, f_off)
    x_complex=(x_real + (x_imag*1j)).astype(np.complex64)
    y_complex=(y_real + (y_imag*1j)).astype(np.complex64)
    temp=stokesI(frequency[band], sample_rate, dm, subbands)
    output=open('/mnt/ucc2_data1/data/giesec/crab/coherent_dedisp_1936I.fil', 'ab')
    output.write(temp)
    output.close()
    if (i % 8203104*N==0):
        band += 1
