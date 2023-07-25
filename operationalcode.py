# -*- coding: utf-8 -*-
"""
SCOTTISH WILDFIRE AND MUIRBURN MAPPING PROJECT: PHASE 2
operationalcode.ipynb

'''
Summary: 
This module contains code for calculating burn locations in Scotland, using Sentinel 2 Analysis Ready Data held on CEDA.
Description:
This code works through a set of CEDA directories and searches for images. 
The list of images is sorted by date. A cleaned, sorted list of files is fed into functions that read the imagery as pre-burn 
and post-burn datasets and passes that to functions that calculate NBR, NBR2 and SAVI. These outputs are then thresholded to 
create a seed layer, ready for region growing. 
Finally a vector shapefile of burn pixels and a logfile are exported.
Contributors: 
    Duncan Blake, NatureScot
    Alastair Graham, Geoger Ltd, @ajggeoger
'''
"""

# --- Imports ---
import logging
import os
import sys
import datetime
import re
import glob
import copy
import gdal
import numpy as np
import rasterio
import matplotlib.pyplot as plt
import fiona
from fiona.crs import from_epsg
from rasterio.features import sieve
from scipy.ndimage import label, generate_binary_structure

import config


# --- Functions ---

def directorycheck(wd, od):
    '''
    Checks whether the supplied directories exist. Exits if they do not.
    
    Return:
    Successful check
    Keyword arguments:
    wd -- workingdirectory
    od -- output directory
    '''
    if os.path.isdir(wd) and os.path.isdir(od) == True:
        pass
    else:
        print('--EXITING--')
        print('one or more directories supplied do not exist')
        sys.exit()
        

def getcloudmetadata (granule):
    '''
    This calculates cloud cover on a granule by granule basis based on the proportion of pixels within the mask 
    (i.e. on land, in Scotland, not inland water or lowland agriculture) that are covered by cloud.
    Keyword arguments:
    granule - full pathname to granule
    '''

    gran_root = granule.split('vmsk')[0]
    # create valid pixel filename
    valid_file = gran_root + "valid.tif"

    validpixel = gdal.Open(valid_file)
    s2_pixel = validpixel.GetRasterBand(1).ReadAsArray()

    imagetransform = validpixel.GetGeoTransform()

    cols = validpixel.RasterXSize
    rows = validpixel.RasterYSize

    # open land/sea raster and get transform
    landmask = gdal.Open(config.LANDMASK)
    land_transform = landmask.GetGeoTransform()

    # calculate offsets (upper left corner) for land/sea raster import. NB offsets are in pixels
    xoffset = int((imagetransform[0] - land_transform[0])/(imagetransform[1]))
    yoffset = int((imagetransform[3] - land_transform[3])/(imagetransform[5]))

    # import land/sea area that matches location and size of the valid pixel image
    s2_land = landmask.ReadAsArray(xoff=xoffset, yoff=yoffset, xsize=cols, ysize=rows).astype('uint16')
    
    # create topographic shadow layer
    topo_file = gran_root + "toposhad.tif"
    topomask = gdal.Open(topo_file)
    topoin = topomask.GetRasterBand(1).ReadAsArray()
    s2_topo = np.zeros((topoin.shape)).astype('uint16')
    topo_index = np.nonzero(topoin < 1)
    s2_topo[topo_index[0], topo_index[1]] = 1
    del topo_index
    
    # create and sum array of all pixels that could be assessed for burns if there were no clouds
    valid_land = s2_land * s2_topo * s2_pixel
    sum_valid_land = np.sum(valid_land)
    
    # import cloudmask and convert to zeros/ones where ones are the clouds
    cloud_file = gran_root + "clouds.tif"
    cloudmask = gdal.Open(cloud_file)
    cloudin = cloudmask.GetRasterBand(1).ReadAsArray()
    s2_cloud = np.ones((cloudin.shape)).astype('uint16')
    cloud_index = np.nonzero(cloudin < 1)
    s2_cloud[cloud_index[0], cloud_index[1]] = 0
    del cloud_index
   
    # calculate number of valid pixels that are masked by cloud
    cloud_land = valid_land * s2_cloud
    sum_cloud_land = np.sum(cloud_land)
    
    # calulate cloud proportion checking first for where there is no valid land
    if sum_valid_land == 0.0:
        cloud_prop = 1.0
        # check for where there is no cloud to avoid division by zero
    elif sum_cloud_land == 0.0:
        cloud_prop = 0.0
    else:
        cloud_prop = sum_cloud_land/sum_valid_land
    
    # print/log results
    y = granule.split('/')[8]
    z = y.split('_utm')[0]
    print('Cloud cover prop for granule ' + str(z) + ' is ' + str(cloud_prop))
    logging.info('Cloud cover prop for granule ' + str(z) + ' is ' + str(cloud_prop))
    
    return cloud_prop


def getdatalist(wd, ymlist, grantoproc, cloud_thresh1, cloud_thresh2):
    '''
    Walks the supplied directories and finds the required imagery
    Return:
    Cleaned list of images to be processed on this run
    Keyword arguments:
    wd -- working directory
    ymlist -- list of year/month combinations to be processed
    grantoproc -- list of granules to be processed
    cloud_thresh1 -- cloud threshold for first run
    cloud_thresh2 -- cloud threshold for 2nd run (0.0 if not required)
    '''
    
    inputfiles = []
    inputfiles2 = []
    inputfilessort2 = []
    
    for yearmonth in ymlist:
        dirstr = wd + yearmonth
        # r = root, d = directories, f = filenames
        for r, d, f in os.walk(dirstr, followlinks=True):
            for name in glob.fnmatch.filter(f, '*osgb_vmsk_sharp_rad_srefdem_stdsref.tif'): 
                # check granule is a Scottish one to be processed
                if name.split('_')[3] in grantoproc:
                    # get image orbit
                    orbit = name.split('_')[4]
                    # check granule falls within cloud threshold
                    cloud_cover = getcloudmetadata(os.path.join(r, name))
                    if cloud_cover < cloud_thresh1:
                        # create [imagename, imagepath, granule, orbit, date]
                        paramlist = [name, r, name.split('_')[3], orbit, name.split('_')[1]] 
                        inputfiles.append(paramlist)
                    # check a second run is required
                    if cloud_thresh2 >0.0:
                        if cloud_cover < cloud_thresh2:
                            # create [imagename, imagepath, granule, orbit, date]
                            paramlist2 = [name, r, name.split('_')[3], orbit, name.split('_')[1]] 
                            inputfiles2.append(paramlist2)
                      
    # sort by granule and date
    inputfilessort = sorted(sorted(inputfiles, key = lambda x : x[4], reverse = True), key = lambda x : x[2], reverse = False)
    if cloud_thresh2 > 0.0:
        inputfilessort2 = sorted(sorted(inputfiles2, key = lambda x : x[4], reverse = True), key = lambda x : x[2], reverse = False)
   
    return inputfilessort, inputfilessort2

def maskimage(imagename, imagetransform, cols, rows, cloudname, toponame):
    '''
    Masks the image dataset by the land/sea and cloud and topographic shadow masks
    
    Return:
    An array of the chosen bands with cloud, shadow and sea areas as zeroes
    Profile of the input image

    Keyword arguments:
    imagename -- path and filename of image to be processed
    imagetransform, cols, rows - parameters for importing area of land raster to match the granule extent
    cloudname -- path and filename of cloud cover dataset to be processed
    toponame -- path and filename of topographic shadow dataset to be processed
    ''' 
    # open land/sea raster and get transform
    landmask = gdal.Open(config.LANDMASK)
    land_transform = landmask.GetGeoTransform()

    # calculate offsets (upper left corner) for land/sea raster import. NB offsets are in pixels
    xoffset = int((imagetransform[0] - land_transform[0])/(imagetransform[1]))
    yoffset = int((imagetransform[3] - land_transform[3])/(imagetransform[5]))

    # import land/sea area that matches location and size of the S2 image
    s2_land = landmask.ReadAsArray(xoff=xoffset, yoff=yoffset, xsize=cols, ysize=rows).astype('uint16')

    # import cloudmask and convert to zeros/ones
    cloudmask = gdal.Open(cloudname)
    cloudin = cloudmask.GetRasterBand(1).ReadAsArray()
    s2_cloud = np.zeros((cloudin.shape)).astype('uint16')
    cloud_index = np.nonzero(cloudin < 1)
    s2_cloud[cloud_index[0], cloud_index[1]] = 1
    
    # import topomask and convert to zeros/ones
    topomask = gdal.Open(toponame)
    topoin = topomask.GetRasterBand(1).ReadAsArray()
    s2_topo = np.zeros((topoin.shape)).astype('uint16')
    topo_index = np.nonzero(topoin < 1)
    s2_topo[topo_index[0], topo_index[1]] = 1
    
    # read in bands required from image using rasterio
    with rasterio.open(imagename, 'r') as s2_image:
        s2_profile = s2_image.profile
        s2_array = s2_image.read((3,7,9,10)).astype('uint16') # double bracket required for 3D array
        # multiply 3D array of bands by two mask arrays to make all area of sea and cloud to be zero.
        s2_array = s2_array * s2_cloud * s2_land * s2_topo  

    del s2_cloud, s2_land, cloudin, cloud_index, topoin, topo_index, s2_topo

    return s2_array, s2_profile

def nbr(swir1, nir):
    '''
    Calculates Normalised Burn Ratio (NBR).
    The NBR equation is the opposite way round to normal but this is as recommended by Filiponi: to make the behavior of the indices consistent i.e. the values of all the indices increase if fire has occurred. 
    Return:
    NBR image
    Keyword arguements:
    swir1 -- short wave infra-red 1 band 
    nir -- near infra-red band
    '''
    nbr = np.zeros(nir.shape).astype(rasterio.float32)
    mask = nir*swir1!=0 # Deal with divide by zero
    nbr[mask] = ((swir1[mask] - nir[mask])/(swir1[mask] + nir[mask]))

    return nbr
    

def nbr2(swir2, swir1):
    '''
    Calculates NBR2
    Return:
    NBR2 image
    Keyword arguements:
    swir2 -- short wave infra-red 2 band 
    swir1 -- short wave infra-red 1 band 
    '''
    nbr2 = np.zeros(swir1.shape).astype(rasterio.float32)
    mask = swir1*swir2!=0 # Deal with divide by zero
    nbr2[mask] = ((swir2[mask] - swir1[mask])/(swir2[mask] + swir1[mask]))

    return nbr2

def savi(nir, red, L = 0.5):
    '''
    Calculates Soil Adjusted Vegetation Index (SAVI) 
    The default version of paramer L is set to 0.5. Specify a value in the function call to over-ride this.
    
    Return:
    SAVI image
    Keyword arguements:
    nir -- near infra-red band 
    red -- red band
    L -- equation parameter
    '''
    savi = np.zeros(nir.shape).astype(rasterio.float32)
    mask = nir*red!=0 # Deal with divide by zero
    savi[mask] = (-1 * (1.5 * ((nir[mask] - red[mask]) / (nir[mask] + red[mask] + L))))
  
    return savi

def threshold_imgs(dsavi, postnbr, dnbr2, thresholds):
    '''
    Uses a specified dictionary of thresholds applied to three input images to create a layer of seed areas ready for region growing. 
    
    Return:
    Seed area array
    Keyword arguements:
    dsavi -- pre/post difference in SAVI  
    postnbr -- post NBR 
    dnbr2 -- pre/post difference in NBR2
    thresholds -- dictionary of thresholds
    '''

    # Create an array of zeros the same shape as the input raster
    reclassArray = np.zeros(postnbr.shape)

    # where dsavi is greater than the median and post fire image NBR is greater than the mean 
    reclassArray[np.where((dsavi>=(thresholds['threshdsavi'])) & (postnbr>=(thresholds['threshpostnbr'])))] = 1 
    
    # attempt to solve issue of edges of clouds being falsely identified which have high values in dnbr2
    reclassArray[np.where(dnbr2>=(thresholds['threshdnbr2']))] = 0
    reclassArray[np.where(dnbr2<=(thresholds['threshdnbr2_shad']))] = 0
    
    # rasterio function to exclude clumps of pixels smaller than 3.  Diagonally joined pixels are allowed.
    sievedArray = rasterio.features.sieve(reclassArray.astype(rasterio.uint8), size=3, connectivity=8)

    return sievedArray

def grow_burn(dsavi, postnbr, dnbr2, thresholds):
    '''
    Uses a specified dictionary of thresholds applied to three input images to `grow` burn areas from the seed areas. 
    
    Return:
    Burn area array
    Keyword arguements:
    dsavi -- pre/post difference in SAVI  
    postnbr -- post NBR 
    dnbr2 -- pre/post difference in NBR2
    thresholds -- dictionary of thresholds
    '''
    # Create an array of zeros the same shape as the input raster
    extendArray = np.zeros(postnbr.shape)

    # reclassify the second array to contain extended burn pixels
    extendArray[np.where((dsavi>=(thresholds['dsaviq1thresh'])) & (postnbr>=(thresholds['postnbrq1thresh'])))] = 1 
    # solve issue of edges of clouds being falsely identified
    extendArray[np.where(dnbr2>=(thresholds['cloudthresh']))] = 0
    extendArray[np.where(dnbr2<=(thresholds['shadowthresh']))] = 0

    # rasterio function to exclude clumps of pixels smaller than 3.  Diagonally joined pixels are allowed.
    burnedArray = rasterio.features.sieve(extendArray.astype(rasterio.uint8), size=3, connectivity=8)

    return burnedArray

def burn_intersect(seed, burn):
    '''
    Extracts the full extent of burns where the 'grown' burn areas intersect with seed burn areas. 
    
    Return:
    Burn extents array

    Keyword arguements:
    seed -- array of seed/core burn patches - output from def threshold_imgs
    burn -- array of wider possible burn areas - output from def grow_burn
    '''
    # create structure to allow diagonally joined pixels to be included in part of same patch
    s = generate_binary_structure(2,2) 
    # create numberd 'zones' of contiguous burn pixels
    burn_labeled, num_features = label(burn, structure=s) 

    # create deep copy as burn_labeled is required later
    seed_labeled = copy.deepcopy(burn_labeled)

    # convert burn zones back to zero if they are not in the core dataset
    seed_labeled[np.where(seed == 0)] = 0 
    
    # create 1D array of unique zone numbers
    unique_array = np.unique(seed_labeled) 

    # create boolean array which is TRUE where an element in burn_labeled is in the unique array
    ix = np.isin(burn_labeled, unique_array)

    # multiplying arrays converts false returns (areas not in core burns) to zeros)
    intersectArray = burn_labeled * ix

    return intersectArray


def gran_process (toprocess):

    '''
    Main function for carrying out comparisons between a list of granules. 
    
    Return:
    Nothing - writes outputs to an open shapefile

    Keyword arguements:
    toprocess -- list of granules to process
    '''     
    # get the last month in the list so that these are not used as post fire images later
    last_gran = toprocess[-1]
    last_date = last_gran[4]
    # gets year and month for year long processing runs
    last_month = last_date[2:6]

    x = 0
    y = 1
    newpost = True
    while x < (len(toprocess)-1):
        # get new cropped post-fire image
        if newpost:
            postlist = toprocess[x]
            
            # check not in last month, in which case break out the while loop
            date = postlist[4]
            month = date[2:6]
            if month == last_month:
                break

            # gets cloud name
            names = postlist[0].split('vmsk')[0]
            cloudname_post = names + "clouds.tif"
            
            # gets topographic shadow
            tnames = postlist[0].split('vmsk')[0]
            toponame_post = tnames + "toposhad.tif"

            # open Sentinel 2 image to get transform
            postimage = gdal.Open(os.path.join(postlist[1], postlist[0]))
            postim_transform = postimage.GetGeoTransform()

            no_cols_post = postimage.RasterXSize
            no_rows_post = postimage.RasterYSize

            post_array, post_profile = maskimage(os.path.join(postlist[1], postlist[0]), postim_transform, no_cols_post, no_rows_post, os.path.join(postlist[1], cloudname_post), os.path.join(postlist[1], toponame_post))

        # PRE FIRE IMAGE
        prelist = toprocess[x+y]

        # gets cloud name
        names = prelist[0].split('vmsk')[0]
        cloudname_pre = names + "clouds.tif"
     
        # gets topographic shadow
        tnames = prelist[0].split('vmsk')[0]
        toponame_pre = tnames + "toposhad.tif"

        # open Sentinel 2 image to get transform
        preimage = gdal.Open(os.path.join(prelist[1], prelist[0]))
        preim_transform = preimage.GetGeoTransform()
        no_cols_pre = preimage.RasterXSize
        no_rows_pre = preimage.RasterYSize
                 
        pre_array, pre_profile = maskimage(os.path.join(prelist[1], prelist[0]), preim_transform, no_cols_pre, no_rows_pre, os.path.join(prelist[1], cloudname_pre), os.path.join(prelist[1], toponame_pre))
        
        # check if pre and post fire images have the same granule size for 2022-2023 change in particular
        if (no_cols_post == no_cols_pre) and (no_rows_post == no_rows_pre):
            print("Dimensions of images are equal")
            equalDims = True
            # go ahead with rest of code
        else:
            print("Dimensions are NOT equal")
            equalDims = False
            # find each image's bounding box 
            # r1 has left, top, right, bottom of dataset's bounds in geospatial coordinates. 
            r1 = [preim_transform[0], preim_transform[3], preim_transform[0] + (preim_transform[1] * no_cols_pre), preim_transform[3] + (preim_transform[5] * no_rows_pre)] 
            r2 = [postim_transform[0], postim_transform[3], postim_transform[0] + (postim_transform[1] * no_cols_post), postim_transform[3] + (postim_transform[5] * no_rows_post)] 
            print ('bounding box: %s' % str(r1))
            print ('bounding box: %s' % str(r2))

            # find intersection between bounding boxes (max left, min top, min right, max bottom)
            intersection = [max(r1[0], r2[0]), min(r1[1], r2[1]), min(r1[2], r2[2]), max(r1[3], r2[3])] 
            print ('intersection coords: %s' % str(intersection))

            # left offset in number of pixels
            colStart1 = abs(int(round((intersection[0]-r1[0])/preim_transform[1])))
            # top offset in number of pixels
            rowStart1 = abs(int(round((intersection[1]-r1[1])/preim_transform[5])))
            # right intersection - left intersection divided by pixel width and add on x offset
            colEnd1 = abs(int(round((intersection[2]- (intersection[0]))/preim_transform[1])) + colStart1)
            # bottom intersection - top intersection divided by pixel width and add on y offset
            rowEnd1 = abs(int(round((intersection[3]- (intersection[1]))/preim_transform[5])) + rowStart1)

            colStart2 = abs(int(round((intersection[0]-r2[0])/postim_transform[1])))
            rowStart2 = abs(int(round((intersection[1]-r2[1])/postim_transform[5])))
            colEnd2 = abs(int(round((intersection[2]- (intersection[0]))/postim_transform[1])) + colStart2)
            rowEnd2 = abs(int(round((intersection[3]- (intersection[1]))/postim_transform[5])) + rowStart2)

            print("Pre fire image array: column start: %s , column end: %s, row start: %s, row end: %s" % (str(colStart1),str(colEnd1),str(rowStart1),str(rowEnd1)))
            print("Post fire image array: column start: %s , column end: %s, row start: %s, row end: %s" % (str(colStart2),str(colEnd2),str(rowStart2),str(rowEnd2)))        
            
            # slice arrays
            pre_array_isect = pre_array[:, rowStart1:rowEnd1, colStart1:colEnd1]
            post_array_isect = post_array[:, rowStart2:rowEnd2, colStart2:colEnd2]
            
        message = ('Processing post-fire granule:', postlist[2], postlist[4], postlist[3], 'against pre-fire granule:', prelist[2], prelist[4], prelist[3])
        print(message)
        logging.info(message)

        # if equal dimensions use the full arrays, if not use the area of overlap or intersect arrays
        if equalDims:
            postr, postnir, postswir1, postswir2 = post_array.astype(float)
            prer, prenir, preswir1, preswir2 = pre_array.astype(float)
        else:
            postr, postnir, postswir1, postswir2 = post_array_isect.astype(float)
            prer, prenir, preswir1, preswir2 = pre_array_isect.astype(float)
            
        print('--CALCULATING postNBR--')
        postnbr = nbr(postswir1, postnir)

        print('--CALCULATING dNBR2--')
        postnbr2 = nbr2(postswir2, postswir1)
        prenbr2 = nbr2(preswir2, preswir1)
        # combine so that nodata values (zero) in either index layer stay nodata in the output
        mask = postnbr2*prenbr2!=0
        # Pre/post NBR2 difference
        dnbr2 = np.zeros(postnbr2.shape).astype(rasterio.float32)
        dnbr2[mask] = postnbr2[mask] - prenbr2[mask]

        print('--CALCULATING dSAVI--')
        postsavi = savi(postnir, postr)
        presavi = savi(prenir, prer)
        # combine so that nodata values (zero) in either index layer stay nodata in the output
        mask = postsavi*presavi!=0
        # Pre/post SAVI difference
        dsavi = np.zeros(postsavi.shape).astype(rasterio.float32)
        dsavi[mask] = postsavi[mask] - presavi[mask]

        # Thresholding
        print('--CALCULATING THRESHOLDING--')
        thresholds = config.THRESHOLD 
        # print('Thresholds used: ', thresholds)
        burnseed = threshold_imgs(dsavi, postnbr, dnbr2, thresholds)

        # Region growing
        print('--CALCULATING BURN REGIONS--')
        thresholds = config.GROW 
        # print('Thresholds used: ', thresholds)
        burnarray = grow_burn(dsavi, postnbr, dnbr2, thresholds)

        # Seed/burn intersection
        print('--CALCULATING ESTIMATED BURN EXTENTS --')
        burnextents = burn_intersect(burnseed, burnarray)
        # REMOVE IN 2024:convert first and last rows and columns to zeros to account for S2 ARD issue where pixel values are spurious
        #lhs
        burnextents[0:, 0] = 0
        #rhs
        burnextents[0:,-1] = 0
        #bottom
        burnextents[-1, 0:] = 0
        #top
        burnextents[0, 0:] = 0

        # Convert raster features to shapes and write to the shapefile
        print('--SAVING DATA--')
        #saveraster(od, postnbr, pre_profile, 'postnbr', prelist[0], postlist[0])
        #saveVector(od, burnextents_sub, pre_profile, preim_transform, prelist[0], postlist[0])
        prename1 = prelist[0]
        postname1 = postlist[0]

        #create output base name
        prename = prename1.split('_')
        postname = postname1.split('_')

        # need to transform GDAL transform into rasterio affine transformation which is structured differently
        im_affine = rasterio.Affine.from_gdal(*preim_transform)
        # make a copy of the array to use as a mask in rasterio.features.shape to avoid all the zeros being turned into polygons
        mask = burnextents.astype(rasterio.uint8)

        shapes = (
                    {'properties': {'raster_val': v, 'pre': prename1, 'post': postname1, 'predate': prename[1], 'postdate': postname[1],'granule': prename[3]}, 'geometry': s}
                    for i, (s, v) 
                    in enumerate(
                        rasterio.features.shapes(burnextents, mask=mask, transform=im_affine)))

        for rec in shapes:
            c.write(rec)

        # if pre-fire image is partial and it is not the last in the sublist then next image becomes pre-fire image and post-fire image stays the same 
        # get orbit number(s) for maximum coverage of granule being processed
        orbit_list = config.MAXORBIT.get(prelist[2], 'granule not in dictionary')

        if (prelist[3] not in orbit_list) and ((len(toprocess) - (x+y)>1)):
            y = y+1
            newpost = False
        # else the prefire granule was a full granule and y is reset to 1, x is incremented, and a new granule becomes the post-fire image
        else: 
            y = 1
            x = x+1
            newpost = True
    print('-----')

# def saveraster(od, datafile, profile, name, prename, postname):
#     '''
#     Saves spatial data to tif file 
    
#     Return:
#     NA
    
#     Keyword arguements:
#     od -- output directory  
#     datafile -- data to be saved 
#     profile -- spatial data profile for rasterio
#     name -- identifier for the dataset allowing decision making in the code
#     prename -- name of the preburn input image
#     postname -- name of the postburn input image
#     '''
    
#     #create output base name
#     prename = prename.split('_')
#     postname = postname.split('_')

#     kwds = profile
    
#     if name == 'burnseed':
#         #Export the thresholded raster
#         kwds.update(dtype=rasterio.uint8,
#             count=1,
#             compress='lzw')
#         outname = prename[0] + prename[1] + prename[3] + prename[4] + '_' + postname[0] + postname[1] + postname[3] + postname[4] + '_' + name + '.tif'

#         with rasterio.open(os.path.join(od, outname), 'w', **kwds) as dst_dataset:
#             dst_dataset.write(datafile, 1)


#     elif name == 'burnarea':
#         #Export the thresholded raster
#         kwds.update(dtype=rasterio.uint8,
#             count=1,
#             compress='lzw')
#         outname = prename[0] + prename[1] + prename[3] + prename[4] + '_' + postname[0] + postname[1] + postname[3] + postname[4] + '_' + name + '.tif'

#         with rasterio.open(os.path.join(od, outname), 'w', **kwds) as dst_dataset:
#             dst_dataset.write(datafile, 1)

#     else:
#         # Change the format driver for the destination dataset to
#         #kwds['driver'] = 'GTiff'
#         kwds['dtype'] = 'float32'
#         kwds['count'] = 1

#         outname = prename[0] + prename[1] + prename[3] + prename[4] + '_' + postname[0] + postname[1] + postname[3] + postname[4] + '_' + name + '.tif'

#         with rasterio.open((os.path.join(od, outname)), 'w', **kwds) as dst_dataset:
#         # Write data to the destination dataset.
#             dst_dataset.write(datafile, 1)    

#======================================================================    
# ========================== MAIN CODE ======================================
# ======================================================================

if __name__ == "__main__":
    
    # Set working directory
    wd = config.ARD_WRKDIR
    
    # Set output directory
    od = config.GWS_DATA
    
    # Set logfile 
    logfile = os.path.join(od, (datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")+'-processing.log'))
    logging.basicConfig(filename=logfile, level=logging.INFO, format='Date-Time : %(asctime)s : Line No. : %(lineno)d - %(message)s')

    # Check directory validity
    directorycheck(wd, od)
    logging.info('Directories validated')
    
    # Start timer
    starttime1 = datetime.datetime.now()
    print('--STARTING PROCESSING--')

    # Get list of granules to process that match granule name, date and cloud coverage criteria   
    # Change to process each unique granule name separately into a separate shapefile.  This means whole process does not have to be 
    # re-run if errors occur
    for gran in config.PROC_GRANULES:
        toprocess, toprocess2 = getdatalist(wd, config.MONTHS_LIST , gran, config.CLOUD, config.CLOUD2)

        print('Processing list constructed' + ' - number of files to process: ' + str(len(toprocess)))
        logging.info('Processing list constructed' + ' - number of files to process: ' + str(len(toprocess)))
    
    
        # If too few images for comparison, exit the program
        if len(toprocess) < 2:
            print('--EXITING--')
            print('Too few images to process in test')
            logging.error('Too few images supplied for processing')
            sys.exit()
    
        if toprocess2:
            print('Second processing list constructed' + ' - number of files to process: ' + str(len(toprocess2)))
            logging.info('Second processing list constructed' + ' - number of files to process: ' + str(len(toprocess2)))
            if len(toprocess2) < 2:
                print('--EXITING--')
                print('Too few images to process in test')
                logging.error('Too few images supplied for processing')
                sys.exit()
            
        # as this list contains granules selected using higher cloud threshold we do not need to worry about second list.
        checklist = copy.deepcopy(toprocess)
    
        # set up shapefile and open to write outputs  during processing
        outname = gran + '.shp'
        crs = from_epsg(27700)
        driver='ESRI Shapefile'
        schema = {'geometry': 'Polygon',
        'properties': {'raster_val': 'int',
                    'pre': 'str',
                    'post': 'str',
                    'predate':'str',
                    'postdate':'str',
                    'granule':'str'}}

        with fiona.open(os.path.join(od,outname),
                         'w',
                          driver=driver,
                          crs=crs,
                          schema=schema) as c:
        
            print('Starting processing list 1')
            logging.info('Starting processing list 1')
            gran_process(toprocess)

            if toprocess2:
                print('Starting processing list 2')
                logging.info('Starting processing list 2')
                gran_process(toprocess2)
            else:
                print('No second run')
                logging.info('No second run')
        
            c.close()
        
        print('--WRITING OUTPUT FILES--')
        logging.info('Writing output file')

    # Stop timer
    endtime1=datetime.datetime.now()
    deltatime1=endtime1-starttime1
    print(("Time to process:  {0}  hr:min:sec".format(deltatime1)))
    logging.info("Time to process:  {0}  hr:min:sec".format(deltatime1))
