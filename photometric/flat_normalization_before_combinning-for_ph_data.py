import numpy as np
import pandas as pd
import astropy.io.fits as pyfits
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob
from astropy.io import fits
import os
import warnings
warnings.filterwarnings('ignore')
from astropy import units as u
from pathlib import Path
import ccdproc as ccdp
from astropy.nddata import CCDData
from astropy.stats import mad_std
from argparse import ArgumentParser


'''
clear up 
'''


parser = ArgumentParser()
parser.add_argument("-input", help="flat images you want to normalize, e.g. 'flat*.fits'", dest = 'input')
parser.add_argument("-method",  help="method: it can be 'mean' or 'median' ", dest="method", )
args = parser.parse_args()


input_file = args.input
method1 = args.method


'''list the flat images'''
flat_only = ccdp.ImageFileCollection('.', glob_include = input_file)

'''read data and estimate the mean/median ADU'''

for i,v in enumerate(flat_only.files):
    flat = pyfits.open(flat_only.files[i])
    primary_flat = flat[0].data
    fmean = np.mean(primary_flat)
    fmedian = np.median(primary_flat)
    print('mean = ', fmean, 'median = ',fmedian)
    
    if method1 == 'mean':
        flat[0].data = flat[0].data / fmean
        
    if method1 == 'median':
        flat[0].data = flat[0].data / fmedian
        
        
    print('save Nor_'+v)
    flat.writeto('Nor_'+v)       
    


        
    
    
