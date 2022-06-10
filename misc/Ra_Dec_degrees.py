#! /usr/bin/python3

import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Convert RA and DEC to Alt and Az')

parser.add_argument("-R", "--input_RA", help="Input as: XXhYYmZZ.ZZZs")
parser.add_argument("-D", "--input_DEC", help="Input as: XXdYYmZZ.ZZZs")
parser.add_argument("-t", "--time", help = "time for conversion to alt,az, MJD")

args = parser.parse_args()
RA = args.input_RA
DEC = args.input_DEC
time_MJD = args.time


#I-LOFAR
lat = 53.09469
lon = -7.92153
z = 75
#LOFAR-SE
#lat = 57.39885
#lon = 11.93029
#z = 18

target_coord = SkyCoord(RA, DEC, frame = 'icrs')
observer_coord = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=z*u.m)
time = Time(time_MJD, format='mjd')

altaz = target_coord.transform_to(AltAz(obstime=time,location=observer_coord))
print(altaz)
