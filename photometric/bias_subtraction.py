
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
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


# In[1]:


# cd /sciproc/disk5/carey1/LOT/SLT/slt20201030/bias-dark


# In[7]:


parser = ArgumentParser()
parser.add_argument("-input", help="input images files, e.g. 'Dark*30S*fits'",dest="input")
parser.add_argument("-exposure", help="float type, input the exopsure of data which you want to process", dest = "exposure", default = None)
parser.add_argument("-zero", help="The file name of a master bias image, e.g. Zero.fits", dest="zero")
args = parser.parse_args()


input_file = args.input
masterbias = args.zero
exposure = args.exposure

if exposure == None:


    data_only1 = ccdp.ImageFileCollection('.', glob_include = input_file)

else:
    
    exposure = float(exposure)
    
    data_only = ccdp.ImageFileCollection('.', glob_include = input_file)
    
    lights_table = data_only.summary['file', 'imagetyp', 'exptime']

    data_only1 = lights_table[lights_table['exptime'] == exposure].to_pandas()


    
'''Readin the master bias'''
master_bias = CCDData.read(masterbias, unit = 'adu')

'''data images bias subtraction'''

try:
    for i, dd in enumerate(data_only1.files):

        print(str(i)+' ',dd,' processing')
        a_data = CCDData.read(dd, unit='adu')
        '''Subtraction'''
        data_d = ccdp.subtract_bias(ccd=a_data, master = master_bias)

        '''Save Dark subtracted flat image'''

        data_d.write('B_'+dd, overwrite=True)
        print('B_'+dd,' saved')
        
except AttributeError:
    for i, dd in enumerate(data_only1['file']):

        print(str(i)+' ',dd,' processing')
        a_data = CCDData.read(dd, unit='adu')
        '''Subtraction'''
        data_d = ccdp.subtract_bias(ccd=a_data, master = master_bias)

        '''Save Dark subtracted flat image'''

        data_d.write('B_'+dd, overwrite=True)
        print('B_'+dd,' saved')
        


# In[11]:


# combined_bias.write.help()

