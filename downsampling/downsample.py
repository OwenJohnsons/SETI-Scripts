import numpy as np 
import blimpy 
import argparse
import re 

parser = argparse.ArgumentParser(description='Downsamples the size of filterbanks')
parser.add_argument('-i', '--inputfiles')

