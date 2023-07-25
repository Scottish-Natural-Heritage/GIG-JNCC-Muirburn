"""
SCOTTISH WILDFIRE AND MUIRBURN MAPPING PROJECT: PHASE 2
This module contains configuration information for use in the code for running on JASMIN to calculate 
burn locations in Scotland, using data held on CEDA
    
Contributors: 
    Duncan Blake, NatureScot
    Alastair Graham, Geoger Ltd, @ajggeoger
"""
# --Main set up--
ARD_WRKDIR = '/neodc/sentinel_ard/data/sentinel_2'           # Input ARD data  

# OPERATIONAL LOCATION
GWS_DATA = '/gws/nopw/j04/jncc_muirburn/users/abdb2/Run_2022_2023' # GWS output location (needs an output folder added) 

LANDMASK = '/gws/nopw/j04/jncc_muirburn/data/Scot_LandMask/muirburn_mask_phase2.tif' # Masking shapefile location

# --Other parameters--
# Image thresholding values for core burn pixels - the variable names are set in the code, but the values can be changed here.
# The type is a helper variable so that users know how it is being applied. It is not used in the code (global == to all images). 
#THRESHOLD = {'threshdsavi': 0.2853, 'threshpostnbr': 0.2395, 'threshdnbr2': 0.8, 'threshdnbr2_shad': 0.045, 'type': 'global'}
THRESHOLD = {'threshdsavi': 0.3, 'threshpostnbr': 0.3, 'threshdnbr2': 0.8, 'threshdnbr2_shad': 0.045, 'type': 'global'}

# Threshold values for 'growing' the burn areas - the variable names are set in the code, but the values can be changed here.
# The type is a helper variable so that users know how it is being applied. It is not used in the code (global == to all images). 
GROW = {'dsaviq1thresh': 0.206748, 'postnbrq1thresh': 0.173447, 'cloudthresh': 0.8, 'shadowthresh': 0.045,'type': 'global'}

# Cloud cover threshold
CLOUD = 0.6
# set CLOUD2 to 0.0 if two runs are not required
CLOUD2 = 0.3

# Scottish granule filter - only process the granules for Scotland. 
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

# Date filter - include month previous to date range required
MONTHS_LIST = ['/2022/04','/2022/05','/2022/06','/2022/07','/2022/08','/2022/09','/2022/10','/2022/11','/2022/12','/2023/01','/2023/02','/2023/03','/2023/04']

