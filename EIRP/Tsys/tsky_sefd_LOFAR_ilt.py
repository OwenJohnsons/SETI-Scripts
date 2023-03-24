#%%
from astropy.coordinates import SkyCoord
from scipy import optimize as opt
import argparse
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pygdsm

# Number of active tiles during observations
N_TILES = 94

# Patch in helpers from tsky_orig.py
stationDiameter = 56.5 # metres
scaleFactor = 1.02 # van Haarlem et al. Tab. B1
c = 299792458 # m/s
rad2Deg = 57.2958 # deg/rad

hwhm = lambda freq: rad2Deg * scaleFactor * (c / (freq * 1e6)) / stationDiameter / 2
gauss2d = lambda sigma, x, y: np.exp(-(np.square(x)/(2 * sigma * sigma) + np.square(y) / (2 * sigma * sigma)))

def powerl(x, a ,b):
	return a * np.power(x, b)

# Kondratiev et al.
def lofar_tinst_range(band = 'HBA', freqs = None, dv = 0.):
	"""
	calculates the LOFAR HBA/LBA average Tinst using polynomial expressions 
	for Tinst from fit to Wijnholds (2011) between frequencies f1 and f2 (in MHz).
	Return value is Tinst in Kelvins.
	If frequency array 'freqs' is given, then average Tinst will be calculated for each
	frequency range f0-f1, f1-f2, f2-f2 of the array and returned value is list of average Tinst's.
	Size of the returned array is smaller by 1 than the size of the input freqs array
	Each pair of frequencies should be either above 100 MHz or below 100 MHz
	"""

	if type(freqs) in [float, int]:
		freqs = [(freqs - dv, freqs + dv)]

	if band.upper() == 'HBA':
		flow=110
		fhigh=250
		fstep=5
		# polynomial coefficients
		T_inst_poly = [6.64031379234e-08, -6.27815750717e-05, 0.0246844426766, -5.16281033712, 605.474082663, -37730.3913315, 975867.990312]
	else:
		print(f"Unknown band {band.upper()}. Exiting.")
		return None
	dpoly = len(T_inst_poly)

	tinsts=[]
	for flower, fupper in freqs:
		tot = 0
		df = fupper - flower
		for ii in range(101):
			freq = flower + ii*(df)/100.
			tinst = 0.0
			for jj in range(dpoly): 
				tinst += T_inst_poly[jj]*(freq)**(dpoly-jj-1)
			tot += tinst
		tot /= 100.
		tinsts.append(tot)

	return tinsts

# 2 tiles out of action -> 94
# Kondratiev et al.
def get_lofar_aeff_max(freqs, nelem=N_TILES, SEPTON = False):
	"""
	Calculate the Aeff using given frequency and EL
	"""

	wavelen = 300.0 / np.array(freqs)
	# HBA
	if np.max(freqs) >= 100.:
		aeff = nelem * (16. if not SEPTON else 1.) * np.minimum((wavelen * wavelen)/3., 1.5625)
	# LBA (LBA_OUTER)
	else:
		aeff = nelem * (wavelen * wavelen)/3.
	return aeff

def calculateBrightness(snr, aeff, beamcorrection, tsys, tsky, tobs, bandwidth = 5, rfiflagged = 0.):
	return snr * (1 * 2 * 1380 *(tsys + tsky) / (aeff / beamcorrection)) / np.sqrt(2 * bandwidth * 1e6 * (1 - rfiflagged) * tobs)


def getCoordinateGrid(centre, frequencies, sampling = 64, nfwhm = 2):
	grids = {frequency: (np.meshgrid(np.linspace(-nfwhm * hwhm(frequency), nfwhm * hwhm(frequency), sampling), np.linspace(-nfwhm * hwhm(frequency), nfwhm * hwhm(frequency), sampling))) for frequency in frequencies}

	centre = centre.galactic

	coords = {}
	for frequency, (gridL, gridB) in grids.items():
		coords[frequency] = centre.spherical_offsets_by(gridL * u.deg, gridB * u.deg)

	return grids, coords

def getSkyRegion(model, coords):
	regions = {}

	for frequency, coord in coords.items():
		regions[frequency] = model.get_sky_temperature(coord, frequency)

	return regions

def applyBeamGuassian(grids, temps, plot = False):
	gaussians = {frequency: gauss2d(hwhm(frequency), grid[0], grid[1]) for frequency, grid in grids.items()}


	convTemp = {}
	for (frequency, gaussian), temp, grid in zip(gaussians.items(), temps.values(), grids.values()):
		convTemp[frequency] = np.sum(np.multiply(temp, gaussian)) / np.sum(gaussian)

		if plot:
			plt.title(f"Raw Sky Temperatures @ {frequency:.3g} MHz")
			plt.xlabel("l [deg]")
			plt.ylabel("b [deg]")
			plt.pcolormesh(grid[0], grid[1], temp)
			plt.colorbar(label = "Temperature [K]")
			plt.scatter(0, 0, alpha = 0.2)
			plt.figure()
			plt.title(f"Beam-Convolved Sky Temperatures @ {frequency:.3g} MHz ({convTemp[frequency]:.3g}K)")
			plt.xlabel("l [deg]")
			plt.ylabel("b [deg]")
			plt.pcolormesh(grid[0], grid[1], np.multiply(temp, gaussian))
			plt.colorbar(label = "Contributed Temperature [K]")
			plt.scatter(0, 0, alpha = 0.2)
			plt.show()



	pars, cov = opt.curve_fit(powerl, np.fromiter(convTemp.keys(), dtype = float), np.fromiter(convTemp.values(), dtype = float))

	return pars, convTemp


def getSourceTsky(source, frequencies, model = pygdsm.LowFrequencySkyModel(freq_unit = 'MHz'), sampling = 64, nhwhm = 2, plot = False):
	model.generate(frequencies)

	referenceValues = {frequency: model.get_sky_temperature(source, frequency) for frequency in frequencies}

	grids, coords = getCoordinateGrid(source, frequencies, sampling, nhwhm)
	temps = getSkyRegion(model, coords)
	pars, convTemp = applyBeamGuassian(grids, temps, plot = plot)


	return pars, convTemp, referenceValues

def getSEFD(tskys, bandwidth = 1, tobs = 1e-3, rfiFraction = 0.):
	sefd = {}
	for freq, tsky in tskys.items():
		sefd[freq] = calculateBrightness(1., aeff = get_lofar_aeff_max(freq), beamcorrection = 1.0, tsys = lofar_tinst_range('HBA', freqs = freq, dv = bandwidth), tsky = tsky, tobs = tobs,  bandwidth = bandwidth, rfiflagged = rfiFraction).item()
	return sefd

def getSEFD_bandavg(args, frequencies, vsamp = 5.):
	width = 1e-3
	freqs = np.arange(frequencies[0] + vsamp / 2, frequencies[1] - vsamp / 2 + vsamp, vsamp)[:, np.newaxis]
	bandFreqs = np.hstack([freqs, freqs])
	bandFreqs[:, 0] -= vsamp / 2
	bandFreqs[:, 1] += vsamp / 2


	aeffAvg = np.mean(get_lofar_aeff_max(freqs))
	tsysAvg = np.mean(lofar_tinst_range('HBA', bandFreqs))

	#print(f"aeffAvg: {get_lofar_aeff_max(bandFreqs)} -> {aeffAvg}")
	#print(f"tsysAvg: {lofar_tinst_range('HBA', bandFreqs)} -> {tsysAvg}")

	source = SkyCoord(args.ra, args.dec, unit = 'hourangle, degree')
	res = getSourceTsky(source, freqs[:, 0].tolist(), model = skyModels[args.model](freq_unit = 'MHz'), sampling = args.samples, nhwhm = args.nhwhm, plot = args.plot)
	tskyAvg = np.mean(list(res[1].values()))

	#print(f"tskyAvg: {list(res[1].values())} -> {tskyAvg}")

	sefd = calculateBrightness(1., aeff = aeffAvg, beamcorrection = 1.0, tsys = tsysAvg, tsky = tskyAvg, tobs = width, bandwidth = frequencies[1] - frequencies[0], rfiflagged = args.rfi_frac)
	print(f"SEFD {np.abs(np.diff(frequencies)).items()}MHz {width / 1e-3}ms: {sefd:.3g} [Jy]")

	return sefd

def getSensitivityLimits(sefd, snr, width_ms, bandwidth_MHz):
	sensitivity = {}
	for freq, sefdv in sefd.items():
		sensitivity[freq] = sefdv * snr / np.sqrt(width_ms * bandwidth_MHz)
	return sensitivity


skyModels = {
	'LFSS': pygdsm.LowFrequencySkyModel,
	'GSM2008': pygdsm.GlobalSkyModel,
	'GSM2016': pygdsm.GlobalSkyModel2016,
	'HASLAM': pygdsm.HaslamSkyModel,
}

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = "Generate sky temperatures for a given RA and DEC and frequency")
	source = parser.add_mutually_exclusive_group(required = True)
	source.add_argument("--ra", '-r', type = str, help = "Right Ascension in hh:mm:ss.s format.")
	source.add_argument("--list", '-l', type = str, help = "File containing lines with format \"{srcName} {RA in rad} {Dec in rad}\\n\". USING --LIST WILL ONLY GENERATE AN OUTPUT TSKY FOREACH SOURCE.")

	parser.add_argument("--dec", '-d', type = str, help = "Declination in dd:mm:ss.s format.")

	parser.add_argument("--output", '-o', default = "./tsky_output.pkl", type = str, help = "Path to output pickle'd dictionary of source Tsky variables")

	parser.add_argument("--freqs", '-f', default = [100, 150, 200], nargs = '+', type = float, help = "Frequencies to sample [MHz].")
	parser.add_argument("--plot", '-p', default = False, action = 'store_true', help = "Whether or not to plot the inspected region of the sky.")
	parser.add_argument("--nhwhm", '-n', default = 2, type = float, help = "Width (in approximated HWHM (half of FWHM) of the beam) to be used during the convolution.")
	parser.add_argument("--samples", '-s', default = 64, type = int, help = "Amount of samples per axis to inspect (forms an n,n grid in galactic coordinate space).")
	parser.add_argument("--model", '-m', default = 'LFSS', choices = ['LFSS', 'GSM2008', 'GSM2016', 'HASLAM'], help = "Choice of SkyModel to generate Tsky values from.")
	parser.add_argument("--rfi_frac", '-R', default = 0., type = float, help = "Fraction of bandwidth that is flagged for RFI (for SEFD/Sensitivity calculations).")

	flags = parser.add_mutually_exclusive_group()

	flags.add_argument("--sefd", '-J', default = False, action = 'store_true', help = "Output SEFD values for a 1MHz / 1ms integrated pulse.")
	flags.add_argument("--sefd_bandavg", '-j', default = False, action = 'store_true', help = "Output SEFD values for a 75MHz 1ms band-averaged pulse. Overwrites the SEFD, sensitivity_width and sensitivity_bw flags.")
	parser.add_argument("--sensitivity_snr", '-S', default = None, type = float, help = "Get system sensitivity for a given signal-to-noise ratio.")
	parser.add_argument("--sensitivity_width", '-W', default = 5, type = float, help = "Get system sensitivity with a given pulse width [ms].")
	parser.add_argument("--sensitivity_bw", '-B', default = 10, type = float, help = "Get system sensitivity with a given bandwidth [MHz].")



	parser.add_argument("--ntiles", default = None, type = int, help = f"Number of HBA tiles used for observation (default: {N_TILES}).")

	args = parser.parse_args()
	if args.ntiles is not None:
		N_TILES = args.ntiles

	if args.list:
		sources = {}
		with open(args.list, 'r') as ref:
			for line in ref.readlines():
				l = line.rstrip('\n').split()
				sources[l[0]] = (float(l[1]), float(l[2]))
		results = {}
		model = skyModels[args.model](freq_unit = 'MHz')
		for source, (ra, dec) in sources.items():
			src = SkyCoord(ra, dec, unit = 'rad')
			results[source] = (getSourceTsky(src, args.freqs, model = model, sampling = args.samples, nhwhm = args.nhwhm, plot = args.plot), src)

			print(f"{source}: {np.mean(list(results[source][0][1].values())):.0f}K")
		with open(args.output, 'wb') as ref:
			pickle.dump(results, ref)
		exit()
	else:
		if args.sefd_bandavg:
			freqs = [110, 185]
			value = getSEFD_bandavg(args, freqs)
			exit()

		source = SkyCoord(args.ra, args.dec, unit = 'hourangle, degree')

		res = getSourceTsky(source, args.freqs, model = skyModels[args.model](freq_unit = 'MHz'), sampling = args.samples, nhwhm = args.nhwhm, plot = args.plot)

		print(f"Power law model = {res[0][0]:.5g} * freq **{res[0][1]:.4g}")

		print(f"\n\nFreq [MHz]:\tConv. Temp. [K]\tRaw Temp [K]\tDiff [K]")

		for (freq, val), origval in zip(res[1].items(), res[2].values()):
			print(f"{freq}:\t\t{val:.5g}\t\t{origval:.5g}\t\t{val - origval:.3g}")

		print()

		if args.sefd or args.sensitivity_snr:
			sefd = getSEFD(res[1], rfiFraction = args.rfi_frac)
		if args.sefd:
			print("\n\nFreq [MHz]:\tSEFD [Jy MHz ms]")
			for freq, val in sefd.items():
				print(f"{freq}:\t\t{val:.5g}")
		if args.sensitivity_snr:
			sensitivity = getSensitivityLimits(sefd, args.sensitivity_snr, args.sensitivity_width, args.sensitivity_bw)
			print("\n\nFreq [MHz]:\tSensitivity Limit [Jy]")
			for freq, val in sensitivity.items():
				print(f"{freq}:\t\t{val:.5g}")
