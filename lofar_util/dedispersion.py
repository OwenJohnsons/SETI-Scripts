
import os
import numpy as np

"""
NOT WORKING
Based on Griffin Foster's notebook for incoherent dedispersion
"""



path='/mnt/ucc2_data1/data/giesec/crab/'
output_file=path+'incoherent_dd/1936I_dd.fil'
input_file=path+'raw_filterbank_files/1936_I1.fil'
dm=56.77
f_top=197.55859375
f_off=-0.1953125
nchans=488
tInt=1



# Dispersion constant in MHz^2 s / pc cm^-3
_D = 4.148808e3


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
		freq = np.append(freq, numpy.inf)

	# Delay in s
	tDelay = dm*_D*((1e6/freq)**2 - (1e6/freq.max())**2)

	# Cleanup
	if singleFreq:
		tDelay = tDelay[0]
	return tDelay

def freq_array(nchans, f_top, f_off): #in MHz
    f_min=f_top+(nchans*f_off)
    freqs=np.arange(f_top, f_min, f_off)
    return freqs

def incoherent(freq, waterfall, tInt, dm, boundary='wrap', fill_value=np.nan):
	"""
	Given a list of frequencies in Hz, a 2-D array of spectra as a function of
	time (time by frequency), and an integration time in seconds, perform 
	incoherent dedispersion on the data.

	The 'boundary' keyword is used to control how the boundary is treated.  The
	two options are:
	  * wrap - Wrap the time boundary for each frequency (default)
	  * fill - Fill the data not within the dedispersed time range with
	    the value specified by the 'fill_value' keyword

	.. versionchanged:: 1.0.3
		Added in options for handling the boundary conditions
	"""

	# Validate the boundary mode
	if boundary not in ('wrap', 'fill'):
		raise ValueError("Unknown boundary handling type '%s'" % boundary)

	# Compute the dispersive delay for the given frequency range
	tDelay = delay(freq, dm)

	# Convert the delays to integration periods
	tDelay = np.round(tDelay / tInt)
	tDelay = tDelay.astype(numpy.int32)

	# Roll the various frequency bins by the right amount
	ddWaterfall = waterfall*0.0
	for i,d in enumerate(tDelay):
		ddWaterfall[:,i] = np.roll(waterfall[:,i], -d)
		if boundary == 'fill' and d > 0:
			ddWaterfall[-d:,i] = fill_value

	# Return
	return ddWaterfall


with open(input_file, "rb") as f:
    map0 = np.memmap(f, mode="r", dtype=np.int8)
    file=np.array(map-1).astype(np.int32)
    del map0
f.close()
frequencies=freq_array(nchans, f_top, f_off)
stokesI_dd=incoherent(frequencies, file, tInt, dm)
output=open(output_file, 'ab')
output.write(stokesI_dd)
output.close()


