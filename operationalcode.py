# -*- coding: utf-8 -*-
"""operationalcode.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1mqL5s_fbdsUjcW-uS5U4aNgjrQXiKKby

'''
Summary: 
This module contains code for testing the scaling up on JASMIN of calculating burn locations in Scotland, using Sentinel 2 Analysis Ready Data held on CEDA.
Description:
This code takes in a working and output directory which is validated, and searches for new images (checking the images found against a saved list of previously processed images). The list of processed images is retrieved from a file that is saved in pickle format. The working directory is generated by a crawl of the working directory filesystem. 
The list of images is sorted by date. A check is made that images are larger than 1GB so that only full granules are used (partially covered granules throw up errors when the division is applied). A cleaned, sorted list of files is fed into functions that read the imagery (having checked they are for the same granule) as pre-burn and post-burn datasets and passes that to functions that calculate NBR, NBR2 and SAVI. These outputs are then thresholded to create a seed layer, ready for region growing. 
Finally the imagery datasets and logfiles are exported.
Contributors: 
    Alastair Graham, Geoger Ltd, @ajggeoger
    Duncan Blake, NatureScot
Possible things to do:
- look at cloud masking
- improve code around partial image cover granules 
'''
"""

# --- Imports ---
import logging
import os
import sys
import datetime
import glob
import pickle
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
        
def picklecheck(od):
    '''
    Check if the pickle file exists i.e. the code has been previously run
    
    Return: 
    A list of previously processed images (which may be empty)
    Keyword arguments:
    od -- output directory
    '''
    if os.path.isfile(os.path.join(od, 'imagelist.pkl')):
        # if it exists, load it as proc_list
        infile = open(os.path.join(od, 'imagelist.pkl'),'rb')
        proc_list = pickle.load(infile)
        infile.close()
    else:
        # create proc_list if doesn't exist
        proc_list = []
    return proc_list

def cleanlistfunc(inputfiles, proc_list):
    '''
    Cleans the file list to make sure that duplication does not occur
    
    Return:
    Cleaned list of data to process
    Keyword arguments:
    inputfiles -- list of crawled input files
    proc_list -- list of previously processed images
    '''

    count = 0
    rem, index, x = [], [], []
    for f in inputfiles:
        count +=1
        index.append(count)
        for e in proc_list:
            if e[0]==f[0]:
                rem.append(count)
    
    p = [x for x in index if x not in rem]
    res_list = [inputfiles[i-1] for i in (p)] 
    return res_list #[::-1]



def getdatalist(wd, proc_list, grantoproc, months_out):
    '''
    Walks the supplied directory and finds the required imagery
    Return:
    Cleaned list of images to be processed on this run
    Keyword arguments:
    wd -- working directory
    proc_list -- list of previously processed images
    '''
    
    inputfiles = []

    for r, d, f in os.walk(wd, followlinks=True):
        for name in glob.fnmatch.filter(f, '*osgb_vmsk_sharp_rad_srefdem_stdsref.tif'):
            # create [imagename, imagepath, granule, size, date]
            size = (os.stat(os.path.join(r, name)).st_size)/(1024*1024*1024)
            if name.split('_')[1][4:6] not in months_out:
                if name.split('_')[3] in grantoproc:
                    paramlist = [name, r, name.split('_')[3], size, name.split('_')[1]] # os.path.join(r, name)
                    #print(paramlist)
                    inputfiles.append(paramlist)
    
    # sort by granule and date
    inputfiles2 = sorted(sorted(inputfiles, key = lambda x : x[4]), key = lambda x : x[2], reverse = False)

    # call cleaning function
    res_list = cleanlistfunc(inputfiles2, proc_list)
    # print(res_list)
    return res_list


def countfiles(wd):
    '''
    Walks the supplied directory and counts files to be processed in each folder
    Return:
    Number of images to be processed in each folder in wd
    Keyword arguments:
    wd -- working directory
    '''
    fileno = []
    for r, d, f in os.walk(wd, followlinks=True):
        for name in glob.fnmatch.filter(f, '*osgb_vmsk_sharp_rad_srefdem_stdsref.tif'):
        # os.walk method is used to travel throught the wd.
            fileno.append([r, len(f)])
    
    unique = [] 
    templist = []

    for item in fileno:
        if item[0] in templist:
            continue
        else:
            unique.append(item)
            templist.append(item[0])
                    
    return fileno


def maskimage(imagename, imagetransform, cols, rows, cloudname):
    '''
    Masks the image dataset by the land/sea and cloud masks
    
    Return:
    An array of the chosen bands with cloud and sea areas as zeroes
    Profile of the input image

    Keyword arguments:
    imagename -- path and filename of image to be processed
    imagetransform, cols, rows - parameters for importing area of land raster to match the granule extent
    cloudname -- path and filename of cloud cover dataset to be processed
    ''' 
    # open land/sea raster and get transform
    landmask = gdal.Open(config.LANDMASK)
    land_transform = landmask.GetGeoTransform()

    # calculate offsets for land/sea raster import. NB offsets are in pixels
    xoffset = int((imagetransform[0] - land_transform[0])/(imagetransform[1]))
    yoffset = int((imagetransform[3] - land_transform[3])/(imagetransform[5]))

    # import land/sea area that matches location and size of the S2 image
    s2_land = landmask.ReadAsArray(xoff=xoffset, yoff=yoffset, xsize=cols, ysize=rows).astype('uint16')

    # import cloudmask and convert to zeros/ones
    cloudmask = gdal.Open(cloudname)
    cloudin = cloudmask.GetRasterBand(1).ReadAsArray().astype('uint16')
    s2_cloud = np.zeros((cloudin.shape))
    cloud_index = np.nonzero(cloudin < 1)
    s2_cloud[cloud_index[0], cloud_index[1]] = 1
    
    # read in bands required from image using rasterio
    with rasterio.open(imagename, 'r') as s2_image:
      s2_profile = s2_image.profile
      s2_array = s2_image.read((3,7,9,10)).astype('uint16') # double bracket required for 3D array
      # multiply 3D array of bands by two mask arrays to make all area of sea and cloud to be zero.
      s2_array = s2_array * s2_cloud * s2_land   

    del s2_cloud, s2_land, cloudin, cloud_index

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

    logging.debug('NBR calculated')
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

    logging.debug('NBR2 calculated')
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
  
    logging.debug('SAVI calculated')
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

def saveraster(od, datafile, profile, name, prename, postname):
    '''
    Saves spatial data to tif file 
    
    Return:
    NA
    
    Keyword arguements:
    od -- output directory  
    datafile -- data to be saved 
    profile -- spatial data profile for rasterio
    name -- identifier for the dataset allowing decision making in the code
    prename -- name of the preburn input image
    postname -- name of the postburn input image
    '''
    
    #create output base name
    prename = prename.split('_')
    postname = postname.split('_')

    kwds = profile
    
    if name == 'burnseed':
        #Export the thresholded raster
        kwds.update(dtype=rasterio.uint8,
            count=1,
            compress='lzw')
        outname = prename[0] + prename[1] + prename[3] + prename[4] + '_' + postname[0] + postname[1] + postname[3] + postname[4] + '_' + name + '.tif'

        with rasterio.open(os.path.join(od, outname), 'w', **kwds) as dst_dataset:
            dst_dataset.write(datafile, 1)


    elif name == 'burnarea':
        #Export the thresholded raster
        kwds.update(dtype=rasterio.uint8,
            count=1,
            compress='lzw')
        outname = prename[0] + prename[1] + prename[3] + prename[4] + '_' + postname[0] + postname[1] + postname[3] + postname[4] + '_' + name + '.tif'

        with rasterio.open(os.path.join(od, outname), 'w', **kwds) as dst_dataset:
            dst_dataset.write(datafile, 1)

    else:
        # Change the format driver for the destination dataset to
        #kwds['driver'] = 'GTiff'
        kwds['dtype'] = 'float32'
        kwds['count'] = 1

        outname = prename[0] + prename[1] + prename[3] + prename[4] + '_' + postname[0] + postname[1] + postname[3] + postname[4] + '_' + name + '.tif'

        with rasterio.open((os.path.join(od, outname)), 'w', **kwds) as dst_dataset:
        # Write data to the destination dataset.
            dst_dataset.write(datafile, 1)    


def saveVector(od, burnedArray, profile, transform, prename, postname):
    '''
    Creates and saves the vector layer  
    
    Return:
    NA
    
    Keyword arguements:
    od -- output directory  
    burnedArray -- estimated burn extents data to be saved 
    profile -- spatial data profile for rasterio
    transform -- transform data for writing files
    prename -- name of the preburn input image
    postname -- name of the postburn input image
    '''
    prename1 = prename
    postname1 = postname

    #create output base name
    prename = prename1.split('_')
    postname = postname1.split('_')
    # satname 
    outname = prename[0] + prename[1] + prename[3] + prename[4] + '_' + postname[0] + postname[1] + postname[3] + postname[4] + '.shp'
    # need to transform GDAL transform into rasterio affine transformation which is structured differently
    im_affine = rasterio.Affine.from_gdal(*transform)
    # make a copy of the array to use as a mask in rasterio.features.shape to avoid all the zeros being turned into polygons
    mask = burnedArray.astype(rasterio.uint8)

    shapes = (
                {'properties': {'raster_val': v, 'pre': prename1, 'post': postname1, 'predate': prename[1], 'postdate': postname[1],'granule': prename[3]}, 'geometry': s}
                for i, (s, v) 
                in enumerate(
                    rasterio.features.shapes(burnedArray, mask=mask, transform=im_affine)))
    
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
      for rec in shapes:
        c.write(rec)
      c.close()

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
    logging.basicConfig(filename=logfile, level=logging.DEBUG, format='%(asctime)s %(message)s')

    # Check directory validity
    directorycheck(wd, od)
    logging.debug('Directories validated')

    # Get count of files (toggle on-off set in config file)
    if config.FILECOUNT == 'on':
        file_count = countfiles(wd)

    # Get data and list of processed files
    # First call in any file names that have been processed. Then get unprocessed files, for the granules in PROC_GRANULES, ignoring certain months listed in MONTHS_OUT
    
    proc_list = picklecheck(od)
    toprocess = getdatalist(wd, proc_list, config.PROC_GRANULES, config.MONTHS_OUT)

    print('Processing list constructed' + ' - number of files to process: ' + str(len(toprocess)))
    # print(*toprocess, sep = "\n")
    logging.debug('Processing list constructed')
    # logging.debug(toprocess)

    # Start timer
    starttime1 = datetime.datetime.now()
    print('--STARTING PROCESSING--')

    cleanlist = []
    # Look for full scenes: remove to process all images (what is effect of null data?)
    for j in toprocess:
        if j[3] > 1:
            cleanlist.append(j)

    # Get total number of files to process
    tot2process = len(cleanlist)
    print('Number of full scenes to process: ' + str(tot2process))

    # If too few images for comparison, exit the program
    if len(cleanlist) < 2:
            print('--EXITING--')
            print('Too few images to process in test')
            logging.error('Too few images supplied for processing')

            sys.exit()

    checklist = copy.deepcopy(cleanlist)

    count = 1
    runno = 0
    while len(cleanlist) > 0:
        print('--GETTING DATA--')
        runno = runno+1
        if count == 1:
            postlist = cleanlist.pop()

            # post-fire image
            # gets cloud name
            names = postlist[0].split('_')[:7]
            names.append('clouds.tif')
            s = '_'
            cloudname = s.join(names)
            # print(cloudname)
            
            # open Sentinel 2 image to get transform
            postimage = gdal.Open(os.path.join(postlist[1], postlist[0]))
            print('POST BURN IMAGE')
            print('Name: ', postimage.GetDescription())
            print('Bands: ', postimage.RasterCount)
            print('Width: ', postimage.RasterXSize)
            print('Height: ', postimage.RasterYSize)
            # print('CRS: ', postimage.GetProjection())
            postim_transform = postimage.GetGeoTransform()
            no_cols = postimage.RasterXSize
            no_rows = postimage.RasterYSize

            print('Cropping to land and cloud masks')

            post_array, post_profile = maskimage(os.path.join(postlist[1], postlist[0]), postim_transform, no_cols, no_rows, os.path.join(postlist[1], cloudname))

            logging.debug('POST image data read')
        
        count = 2
        prelist = cleanlist.pop()

        # PRE FIRE IMAGE
        # create associated cloud image name
        names = prelist[0].split('_')[:7]
        names.append('clouds.tif')
        s = '_'
        cloudname = s.join(names)

        # open Sentinel 2 image to get transform
        preimage = gdal.Open(os.path.join(prelist[1], prelist[0]))
        print('PRE BURN IMAGE')
        print('Name: ', preimage.GetDescription())
        print('Bands: ', preimage.RasterCount)
        print('Width: ', preimage.RasterXSize)
        print('Height: ', preimage.RasterYSize)
        # print('CRS: ', preimage.GetProjection())
        preim_transform = preimage.GetGeoTransform()
        no_cols = preimage.RasterXSize
        no_rows = preimage.RasterYSize

        pre_array, pre_profile = maskimage(os.path.join(prelist[1], prelist[0]), preim_transform, no_cols, no_rows, os.path.join(prelist[1], cloudname))

        if prelist[2]==postlist[2]:

            postr, postnir, postswir1, postswir2 = post_array.astype(float)
            #postnir, postswir1 = post_array.astype(float)
            prer, prenir, preswir1, preswir2 = pre_array.astype(float)
            #prenir, preswir1 = pre_array.astype(float)

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


            # Save data
            print('--SAVING DATA--')
            #saveraster(od, postnbr, pre_profile, 'postnbr', prelist[0], postlist[0])
            #saveraster(od, dnbr2, pre_profile, 'dnbr2', prelist[0], postlist[0])
            #saveraster(od, dsavi, pre_profile, 'dsavi', prelist[0], postlist[0])
            #saveraster(od, burnseed, pre_profile, 'burnseed', prelist[0], postlist[0])
            #saveraster(od, burnarray, pre_profile, 'burnarea', prelist[0], postlist[0])
            #saveraster(od, burnextents, pre_profile, 'burnextent', prelist[0], postlist[0])

            saveVector(od, burnextents, pre_profile, preim_transform, prelist[0], postlist[0])

            print('Processed', runno, 'of', tot2process, 'files')
        
        if len(cleanlist) >= 1:
            postlist = prelist
            post_array, post_profile = pre_array, pre_profile
  
    print('--WRITING OUTPUT--')
    logging.debug('Writing output file')


    # pickle file
    checklist.pop()
    with open(os.path.join(od, 'imagelist.pkl'),'wb') as outfile:
        pickle.dump(checklist,outfile)
    
    # text file
    with open(os.path.join(od, 'imagelist.txt'), 'w') as outfiletxt:
        outfiletxt.writelines("%s" % line for line in checklist)

    # Stop timer
    endtime1=datetime.datetime.now()
    deltatime1=endtime1-starttime1
    print(("Time to process:  {0}  hr:min:sec".format(deltatime1)))
    logging.debug("Time to process:  {0}  hr:min:sec".format(deltatime1))



