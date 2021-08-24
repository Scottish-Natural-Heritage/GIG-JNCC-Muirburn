"""
SCOTTISH WILDFIRE AND MUIRBURN MAPPING PROJECT: PHASE 2
This module contains configuration information for use in the code for running on JASMIN to calculate 
burn locations in Scotland, using data held on CEDA
    
Contributors: 
    Duncan Blake, NatureScot
    Alastair Graham, Geoger Ltd, @ajggeoger
"""
# --Main set up--
ARD_WRKDIR = '/neodc/sentinel_ard/data/sentinel_2/2020/04'           # Input ARD data (default is April 2019)  

# OPERATIONAL LOCATION
GWS_DATA = '/gws/nopw/j04/jncc_muirburn/users/abdb2/Phase2_April20' # GWS output location (needs an output folder added) 

#TEST LOCATION for when used as a notebook
#GWS_DATA = '/home/users/abdb2/phase2_test'

# Name for output file for the run
OUT_SHAPE = 'Wildfire_and_muirburn_April2020.shp'

LANDMASK = '/gws/nopw/j04/jncc_muirburn/data/Scot_LandMask/muirburn_mask_phase2.tif' # Landmask shapefile location

# --Other parameters--
# Image thresholding values - the variable names are set in the code, but the values can be changed here.
# The type is a helper variable so that users know how it is being applied. It is not used in the code (global == to all images). 
THRESHOLD = {'threshdsavi': 0.2853, 'threshpostnbr': 0.2395, 'threshdnbr2': 0.8, 'type': 'global'}

# Image thresholding values - the variable names are set in the code, but the values can be changed here.
# The type is a helper variable so that users know how it is being applied. It is not used in the code (global == to all images). 
GROW = {'dsaviq1thresh': 0.206748, 'postnbrq1thresh': 0.173447, 'cloudthresh': 0.8, 'type': 'global'}

# Cloud cover threshold
CLOUD = 0.6

# Scottish granule filter - only process the granules for Scotland. 
#PROC_GRANULES = ['T29UPB']
PROC_GRANULES = ['T29UPB', 'T29VND', 'T29VNE', 'T29VPC', 'T29VPD', 'T29VPE', 'T30UUF', 'T30UUG', 'T30UVF', 'T30UVG', 'T30UWG', 'T30VUH', 'T30VUJ', 'T30VUK', 'T30VUL', 'T30VVH', 'T30VVJ', 'T30VVK', 'T30VVL', 'T30VWH', 'T30VWJ', 'T30VWL', 'T30VWM', 'T30VWN', 'T30VXM', 'T30VXN']

# Dictionary containing the orbit of maximum coverage for each granule.  NB in three instances there are two orbits.
# This cannot account for the odd SPLIT granule.

MAXORBIT = {'T29UPB':['ORB123'],
              'T29VND':['ORB023'],
              'T29VNE':['ORB023'],
              'T29VPC':['ORB123'],
              'T29VPD':['ORB123'],
              'T29VPE':['ORB023'],
              'T30UUF':['ORB080'],
              'T30UUG':['ORB123'],
              'T30UVF':['ORB080'],
              'T30UVG':['ORB080'],
              'T30UWG':['ORB037'],
              'T30VUH':['ORB123'],
              'T30VUJ':['ORB123'],
              'T30VUK':['ORB123'],
              'T30VUL':['ORB023', 'ORB123'],
              'T30VVH':['ORB080'],
              'T30VVJ':['ORB080'],
              'T30VVK':['ORB123'],
              'T30VVL':['ORB123'],
              'T30VWH':['ORB037', 'ORB080'],
              'T30VWJ':['ORB080'],
              'T30VWL':['ORB080'],
              'T30VWM':['ORB080'],
              'T30VWN':['ORB080','ORB123'],
              'T30VXM':['ORB080'],
              'T30VXN':['ORB080'], 
               }

# Date filter for seasonality - do not process 01 Sept - 31 dec inclusive
MONTHS_OUT = ['09', '10', '11', '12']
