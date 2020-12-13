"""
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

"""

# --- Imports ---
import logging
import os
import sys
import datetime
import glob
import pickle
import copy

import numpy as np
import rasterio
import fiona
from rasterio.features import sieve
from rasterio.mask import mask
import geopandas as gpd

import config # config.py configuration parameters


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
        for name in glob.fnmatch.filter(f, '*vmsk_sharp_rad_srefdem_stdsref.tif'):
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
        for name in glob.fnmatch.filter(f, '*vmsk_sharp_rad_srefdem_stdsref.tif'):
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
            print(item)
                    
    return fileno


def getlandmask(landmaskpath):
    '''
    Opens and reads the land mask dataset. 
    Takes in a path to a .shp file that holds only polygons for areas of land for a given place of interest e.g. Scotland.
    
    Return:
    Polygons of the land mass.

    Keyword arguments:
    landmaskpath -- the path to the shapefile to be processed
    '''

    # landmask is path to file from CONFIG
    with fiona.open(landmaskpath, "r") as landmaskshp:
        shapes = [feature["geometry"] for feature in landmaskshp]

    
        logging.debug('landmask data read')
        return shapes


def getcloudmask(cloudname):
    '''
    Opens and reads the cloud mask dataset. 
    Takes in a path to a .tif file that holds the cloudmask for the granule being processed.
    
    Return:
    A mask where non-cloud equals 1

    Keyword arguments:
    cloudname -- the path to the cloudmask to be processed
    '''
    with rasterio.open(cloudname) as clouddataset:
        cloudin = clouddataset.read(1)

        new_mask = np.zeros((cloudin.shape))
        cloud_index = np.nonzero(cloudin < 1)
        new_mask[cloud_index[0], cloud_index[1]] = 1
        # copy the profile and update to integer
        profile = clouddataset.profile.copy()
        profile.update(dtype=rasterio.uint8,
            count=1,
            compress='lzw')  

    with rasterio.open(os.path.join(od, 'cloudmask.tif'), "w", **profile) as dest:
        dest.write(new_mask, 1)



def maskify(image, cloudmask):#, landmask):
    '''
    Masks the band being processed by the cloud mask
    
    Return:
    A dataset containing valid data 

    Keyword arguments:
    image -- the image band to be processed
    cloudmask -- the cloudmask    
    '''    
    maskedimage = image * cloudmask
    logging.debug('Cloud masking complete')
    return maskedimage


def masktheland(dataset, outimagename):
    '''
    Masks the image dataset by the land mask
    
    Return:
    A dataset on file containing all data over land

    Keyword arguments:
    dataset -- the rasterio object to be processed
    '''    
    
    out_image, out_transform = rasterio.mask.mask(dataset, landmask, crop=True)

    out_meta = dataset.meta    
    out_meta.update({"driver": "GTiff",
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        "transform": out_transform})

    # Write to temporary file
    with rasterio.open(os.path.join(od, outimagename), "w", **out_meta) as dest:
        dest.write(out_image)


def pre(imagename, cloudname):
    '''
    Opens and reads the pre fire image. 
    Takes in a path to a .tif file.
    
    Return:
    Individual bands and the image profile.

    Keyword arguments:
    imagename -- the path to the image to be processed
    cloudname -- name of the associated cloud mask
    '''
    with rasterio.open(imagename) as dataset:
        print('PRE BURN IMAGE')
        print('Name: ', dataset.name)
        print('Bands: ', dataset.count)
        print('Width: ', dataset.width)
        print('Height: ', dataset.height)
        print('CRS: ', dataset.crs)

        print('Cropping to land mask')
        masktheland(dataset)
        
    with rasterio.open(os.path.join(od, 'temp.tif')) as dataset:
        profile = dataset.profile.copy()
        transform = dataset.transform

        cloudmask = getcloudmask(cloudname) 
        
        print('Masking for cloud')
        red = maskify(dataset.read(3), cloudmask)
        nir = maskify(dataset.read(7), cloudmask)
        swir1 = maskify(dataset.read(9), cloudmask)
        swir2 = maskify(dataset.read(10), cloudmask)
        
        logging.debug('PRE image data read')
        return red, nir, swir1, swir2, profile, transform


def post(imagename, cloudname):
    '''
    Opens and reads the post fire image. 
    Takes in a path to a .tif file.
    
    Return:
    Individual bands.

    Keyword arguments:
    imagename -- the path to the image to be processed
    cloudname -- name of the associated cloud mask
    '''

    with rasterio.open(imagename) as dataset:
        print('POST BURN IMAGE')
        print('Name: ', dataset.name)
        print('Bands: ', dataset.count)
        print('Width: ', dataset.width)
        print('Height: ', dataset.height)
        print('CRS: ', dataset.crs)

        print('Cropping to land mask')
        masktheland(dataset, 'postimagecrop.tif')

    getcloudmask(cloudname)

    with rasterio.open(os.path.join(od, 'cloudmask.tif')) as dataset:
        print('Cropping cloud to land mask')
        masktheland(dataset, 'postimagecrop_cloud.tif')
        
    with rasterio.open(os.path.join(od, 'postimagecrop.tif')) as postdataset:
        profile = postdataset.profile.copy()
         
        
        print('Masking for cloud')
        red = maskify(postdataset.read(3), dataset)
        nir = maskify(postdataset.read(7), dataset)
        swir1 = maskify(postdataset.read(9), dataset)
        swir2 = maskify(postdataset.read(10), dataset)

        logging.debug('POST image data read')
        return red, nir, swir1, swir2, profile


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
    nbr = ((swir1 - nir)/(swir1 + nir)).astype(rasterio.float32)
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
    nbr2 = ((swir2 - swir1)/(swir2 + swir1)).astype(rasterio.float32)
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
    savi = (-1 * (1.5 * ((nir - red) / (nir + red + L)))).astype(rasterio.float32)
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

    # Create a copy of the postnbr image and reset all values in output raster to 0
    reclassArray = postnbr.copy()
    reclassArray[np.where(reclassArray != 0)] = 0

    # where dsavi is greater than the median and post fire image NBR is greater than the mean 
    reclassArray[np.where((dsavi>=(thresholds['threshdsavi'])) & (postnbr>=(thresholds['threshpostnbr'])))] = 1 
    
    # attempt to solve issue of edges of clouds being falsely identified which have high values in dnbr2
    # rasterio function to exclude clumps of pixels smaller than 3.  Diagonally joined pixels are allowed.
    
    reclassArray[np.where(dnbr2>=(thresholds['threshdnbr2']))] = 0
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

    # read in a single band as a template for the output 
    extendArray = postnbr.copy()

    # reset all values in output raster to 0 
    extendArray[np.where(extendArray != 0)] = 0

    # reclassify the second array to contain extended burn pixels
    extendArray[np.where((dsavi>=(thresholds['dsaviq1thresh'])) & (postnbr>=(thresholds['postnbrq1thresh'])))] = 1 
    # solve issue of edges of clouds being falsely identified
    extendArray[np.where(dnbr2>=(thresholds['cloudthresh']))] = 0

    # rasterio function to exclude clumps of pixels smaller than 3.  Diagonally joined pixels are allowed.
    burnedArray = rasterio.features.sieve(extendArray.astype(rasterio.uint8), size=3, connectivity=8)

    return burnedArray

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


def saveVector(od, sievedArray, burnedArray, profile, transform, prename, postname):
    '''
    Saves spatial data to shp file. Calculates the vector layer at the same time 
    
    Return:
    NA
    
    Keyword arguements:
    od -- output directory  
    sieveArray -- seed burn data to be saved
    burnedArray -- region grown data to be saved 
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


    shapes = (
                {'properties': {'raster_val': v, 'pre': prename1, 'post': postname1, 'predate': prename[1], 'postdate': postname[1],'granule': prename[3]}, 'geometry': s}
                for i, (s, v) 
                in enumerate(
                    rasterio.features.shapes(sievedArray, transform=transform)))


    coreShapesGeoms = list(shapes)
    # convert geoJSON objects to a geopandas data frame
    gpd_coreShapes  = gpd.GeoDataFrame.from_features(coreShapesGeoms)
    gpd_coreShapes = gpd_coreShapes.set_crs(epsg=27700)

    # do the same for extended burn areas
    extendShapes = (
                {'properties': {'raster_val': v}, 'geometry': s}
                for i, (s, v) 
                in enumerate(
                    rasterio.features.shapes(burnedArray, transform=transform)))

    extendShapesGeoms = list(extendShapes)

    # convert geoJSON objects to a geopandas data frame
    gpd_extendShapes  = gpd.GeoDataFrame.from_features(extendShapesGeoms)
    gpd_extendShapes = gpd_extendShapes.set_crs(epsg=27700)

    #Use spatial join to detect polygons in extended burn areas that intersect core burn pixels

    # subset the geodataframes to only include rows that are burns
    gpd_coreBurnShapes = gpd_coreShapes[gpd_coreShapes["raster_val"] == 1.0]
    gpd_extendBurnShapes = gpd_extendShapes[gpd_extendShapes["raster_val"] == 1.0]

    # indexing needs to be reset to enable spatial join to work properly
    gpd_coreBurnShapes = gpd_coreBurnShapes.reset_index(drop=True)
    gpd_extendBurnShapes = gpd_extendBurnShapes.reset_index(drop=True)

    # carry out spatial join
    gpd_spatialJoin = gpd.sjoin(gpd_extendBurnShapes, gpd_coreBurnShapes, how="inner", op='intersects')

    # get rid of attribute columns which are not required
    gpd_spatialJoin = gpd_spatialJoin.drop(columns=['index_right','raster_val_left','raster_val_right'])

    # drop duplicate geometries
    gpd_finalShapes = gpd_spatialJoin.drop_duplicates(subset = 'geometry', keep = 'first')

    # export to shapefile
    gpd_finalShapes.to_file(os.path.join(od,outname), driver='ESRI Shapefile')

    


# ======================================================================    
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
    
    landmask = getlandmask(config.LANDMASK)
    proc_list = picklecheck(od)
    toprocess = getdatalist(wd, proc_list, config.PROC_GRANULES, config.MONTHS_OUT)

    print('Processing list constructed')
    logging.debug('Processing list constructed')
    logging.debug(toprocess)


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
            # create associated cloud image name
            names = postlist[0].split('_')[:7]
            names.append('clouds.tif')
            s = '_'
            cloudname = s.join(names)

            postred, postnir, postswir1, postswir2, postprofile = post(os.path.join(postlist[1], postlist[0]), os.path.join(postlist[1], cloudname))
        
        
        count = 2
        prelist = cleanlist.pop()

        # pre-fire image
        # create associated cloud image name
        names = prelist[0].split('_')[:7]
        names.append('clouds.tif')
        s = '_'
        cloudname = s.join(names)
        prered, prenir, preswir1, preswir2, preprofile, pretransform = pre(os.path.join(prelist[1], prelist[0]), os.path.join(prelist[1], cloudname))
                
        if prelist[2]==postlist[2]:

            #PROCESSING

            print('--CALCULATING postNBR--')
            #prenbr = nbr(preswir1, prenir)
            postnbr = nbr(postswir1, postnir)
            #dnbr = postnbr - prenbr


            print('--CALCULATING dNBR2--')
            # Pre/post NBR2 difference
            dnbr2 = nbr2(postswir2, postswir1) - nbr2(preswir2, preswir1)


            print('--CALCULATING dSAVI--')
            # Pre/post SAVI difference
            dsavi = savi(postnir, postred) - savi(prenir, prered)


            # Thresholding
            print('--CALCULATING THRESHOLDING--')
            thresholds = config.THRESHOLD 
            print('Thresholds used: ', thresholds)
            burnseed = threshold_imgs(dsavi, postnbr, dnbr2, thresholds)

            # Region growing
            print('--CALCULATING BURN REGIONS--')
            thresholds = config.GROW 
            print('Thresholds used: ', thresholds)
            burnarray = grow_burn(dsavi, postnbr, dnbr2, thresholds)

            # Save data
            print('--SAVING DATA--')
            #saveraster(od, postnbr, preprofile, 'postnbr', prelist[0], postlist[0])
            #saveraster(od, dnbr2, preprofile, 'dnbr2', prelist[0], postlist[0])
            #saveraster(od, dsavi, preprofile, 'dsavi', prelist[0], postlist[0])
            saveraster(od, burnseed, preprofile, 'burnseed', prelist[0], postlist[0])
            saveraster(od, burnarray, preprofile, 'burnarea', prelist[0], postlist[0])

            saveVector(od, burnseed, burnarray, preprofile, pretransform, prelist[0], postlist[0])

            print('Processed', runno, 'of', tot2process, 'files')
        


        if len(cleanlist) >= 1:
            postlist = prelist
            postred, postnir, postswir1, postswir2, postprofile = prered, prenir, preswir1, preswir2, preprofile
    


    print('--WRITING OUTPUT--')
    logging.debug('Writing output file')


    # pickle file
    checklist.pop()
    with open(os.path.join(od, 'imagelist.pkl'),'wb') as outfile:
        pickle.dump(checklist,outfile)
    
    # text file
    with open(os.path.join(od, 'imagelist.txt'), 'w') as outfiletxt:
        outfiletxt.writelines("%s" % line for line in checklist)

    # clean up temp file
    if os.path.exists(os.path.join(od, 'temp.tif')):
        os.remove(os.path.join(od, 'temp.tif'))
    else:
        print("The file does not exist")
        pass

    # Stop timer
    endtime1=datetime.datetime.now()
    deltatime1=endtime1-starttime1
    print(("Time to process:  {0}  hr:min:sec".format(deltatime1)))
    logging.debug("Time to process:  {0}  hr:min:sec".format(deltatime1))
