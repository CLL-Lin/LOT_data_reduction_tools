
# coding: utf-8

# In[1]:


import numpy as np
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


# In[ ]:
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-input", help="input bias image, e.g. 'Bias*fits'", dest = 'input')
parser.add_argument("-method",  help="method: it can be 'average' or 'median' ", dest="method", )
parser.add_argument("-output",  help="set a name for output combined bias image, e.g. Zero.fits ", dest="output", )
args = parser.parse_args()

input_file = args.input
method1 = args.method
output_file = args.output

bias_only = ccdp.ImageFileCollection('.', glob_include = input_file)


'''combine bias images'''
combined_bias = ccdp.combine(bias_only.files,method = method1,sigma_clip = True, sigma_clip_low_thresh = 5, sigma_clip_high_thresh = 5,sigma_clip_func = np.ma.median, sigma_clip_dev_func = mad_std,mem_limit=350e6,unit='adu')

'''add header combined'''
combined_bias.meta['combined'] = True

'''save the combined bias image'''
combined_bias.write(output_file)

print(output_file,' saved')

