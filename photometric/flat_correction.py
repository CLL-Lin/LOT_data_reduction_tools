
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
parser.add_argument("-input", help="input images you want to process, e.g. 'Proxima_Cen*.fits'", dest = 'input')
parser.add_argument("-flat",  help="input a master flat image file name, e.g. Flat.fits", dest="flat", )
parser.add_argument("-exposure", help="float type, input the exopsure of data which you want to process", dest = "exposure", default = None)
args = parser.parse_args()

input_file = args.input
masterflat = args.flat
exposure = args.exposure

if exposure == None:


    light1 = ccdp.ImageFileCollection('.', glob_include = input_file)

else:
    
    exposure = float(exposure)
    
    light = ccdp.ImageFileCollection('.', glob_include = input_file)
    
    lights_table = light.summary['file', 'imagetyp', 'exptime']

    light1 = lights_table[lights_table['exptime'] == exposure].to_pandas()

    
'''Readin the master flat'''
master_flat = CCDData.read(masterflat, unit = 'adu')


'''Flat correction'''

try:
    for i,v in enumerate(light1.files):

        '''Read light data'''
        a_light = CCDData.read(light1.files[i], unit='adu')
        print(str(i)+'. Raw data: ', light1.files[i])

        '''Flat correction'''
        fd_a_light = ccdp.flat_correct(a_light, master_flat)

        '''Save reduecd light data'''
        if os.path.exists('F_'+light1.files[i]) == False:
            fd_a_light.write('F_'+light1.files[i])
            print('F_'+light1.files[i]+' data saved')
        else:
            print('This data already exists')

except AttributeError:

    for i,v in enumerate(light1['file']):

        '''Read light data'''
        a_light = CCDData.read(v, unit='adu')
        print(str(i)+'. Raw data: ', v)

        '''Flat correction'''
        fd_a_light = ccdp.flat_correct(a_light, master_flat)

        '''Save reduecd light data'''
        if os.path.exists('F_'+v) == False:
            fd_a_light.write('F_'+v)
            print('F_'+v+' data saved')
        else:
            print('This data already exists')


# In[11]:


# combined_bias.write.help()

