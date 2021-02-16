"""
This module contains configuration information for use in the code for running on JASMIN to calculate 
burn locations in Scotland, using data held on CEDA
    
Contributors: 
    Alastair Graham, Geoger Ltd, @ajggeoger
    Duncan Blake, NatureScot
"""
# --Main set up--
ARD_WRKDIR = '/neodc/sentinel_ard/data/sentinel_2/2020/04'           # Input ARD data (default is April 2019)  
GWS = '/gws/nopw/j04/jncc_muirburn'                                  # Group workspace
GWS_DATA = '/gws/nopw/j04/jncc_muirburn/users/abdb2/duncan_v2'                       # GWS output location (needs an output folder added) 
LANDMASK = '/gws/nopw/j04/jncc_muirburn/data/Scot_LandMask/landmask_water.tif' # Landmask shapefile location

# --Other parameters--
# Image thresholding values - the variable names are set in the code, but the values can be changed here.
# The type is a helper variable so that users know how it is being applied. It is not used in the code (global == to all images). 
THRESHOLD = {'threshdsavi': 0.2853, 'threshpostnbr': 0.2395, 'threshdnbr2': 0.8, 'type': 'global'}

# Image thresholding values - the variable names are set in the code, but the values can be changed here.
# The type is a helper variable so that users know how it is being applied. It is not used in the code (global == to all images). 
GROW = {'dsaviq1thresh': 0.206748, 'postnbrq1thresh': 0.173447, 'cloudthresh': 0.8, 'type': 'global'}

# Toggle the file count function on and off. Value can be 'off' or 'on'
FILECOUNT = 'on'

# Cloud cover threshold (for use in future versions, not called in current code)
#CLOUD = 0.9 

# Scottish granule filter - only process the granules for Scotland. 
PROC_GRANULES = ['T29UPB', 'T29VND', 'T29VNE', 'T29VPC', 'T29VPD', 'T29VPE', 'T30UUF', 'T30UUG', 'T30UVF', 'T30UVG', 'T30UWG', 'T30VUH', 'T30VUJ', 'T30VUK', 'T30VUL', 'T30VVH', 'T30VVJ', 'T30VVK', 'T30VVL', 'T30VWH', 'T30VWJ', 'T30VWL', 'T30VWM', 'T30VWN', 'T30VXM', 'T30VXN']

# Date filter for seasonality - do not process 01 Sept - 31 dec inclusive - Scottish ARD only starts in Feb 2019
MONTHS_OUT = ['09', '10', '11', '12']
