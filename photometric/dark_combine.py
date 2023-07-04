import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# %matplotlib inline

import glob
from astropy.io import fits
from scipy.ndimage import interpolation as interp

# from skimage.feature.register_translation import (register_translation, _upsampled_dft)

import os
## This turns off warnings: not a great way to code
## But when we show the images, sometimes we're taking the logarithm of zero and it doesn't like that
## Which would matter if we were doing math, but we're just taking a look at images, so we can ignore it. 
import warnings
warnings.filterwarnings('ignore')

from astropy import units as u
from pathlib import Path
import ccdproc as ccdp
from astropy.nddata import CCDData
from astropy.stats import mad_std


from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-input", help="input dark image, e.g. 'Dark*30S.fits'", dest = 'input')
parser.add_argument("-method",  help="method: it can be 'average' or 'median' ", dest="method", )
parser.add_argument("-output", help="set a name for output combined dark image, e.g. Dark_30S.fits ", dest="output", )
args = parser.parse_args()

input_file = args.input
method1 = args.method
output_file = args.output

darks_only = ccdp.ImageFileCollection('.', glob_include = input_file)


'''combine dark images'''
combined_darks = ccdp.combine(darks_only.files,\
                         method = method1,\
                         sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,\
                         sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,\
                         mem_limit=350e6,unit='adu')

'''add header combined'''
combined_darks.meta['combined'] = True

'''save the combined dark image'''
combined_darks.write(output_file)

print(output_file,' saved')