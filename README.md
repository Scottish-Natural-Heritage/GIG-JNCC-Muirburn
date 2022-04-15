# GIG-JNCC-Muirburn
Code linked to the NatureScot/JNCC wildfire and muirburn monitoring project, designed to run on JASMIN using analysis ready Sentinel-2 data from the CEDA archive.

## Summary
This code is designed to run on a JASMIN Science Server (provided by JNCC for this particular project). It is an attempt to scale up the use of Sentinel-2 (S2) Analysis Ready Data (ARD) - created by the JNCC and hosted on CEDA - for the purposes of mapping wildfire and muirburn areas in Scotland. The methodology for the assessment of burn detection was developed initially by NatureScot EO analysts using trial sites in the Cairngorms and the Isle of Skye with some aspects amended when the code was run nationally.

The code runs through the following steps:
* Obtains a list of granules to process that match a Scottish granule name and fall within a date range and cloud coverage criteria given in the config.py file.  Note cloud coverage percentages are refined first by combining the cloud mask with a land cover mask giving areas of interest (i.e. area of Scotland excluding sea, water bodies and agricultural areas)
* Obtains a second list of granules if a second cloud cover threshold has been defined  (this increases the probability of burn detection though takes longer to process)
* Loops through each granule name in turn (i.e. 'T29UPB' to 'T30VXN') and creates a sublist of granules to process.
* Identifies an image pair to process (termed pre and post fire image) and masks out areas of water ,agriculture, sea, cloud and cloud shadow and topographic shadow on each.
* Calculates indices on remaining areas of the granules (NBR, NBR2, SAVI) and the differences in NBR2 and SAVI.
* Thresholds the index layers using the values in the config.py file to identify core or seed burn pixels.
* 'Grows' these burn areas using a second set of threshold values.
* Exports possible burn areas to a shapefile

## How To
The code is presented as a single script (operationalcode.py) and configuration file (config.py). To run the code, log into the JASMIN Science Server and clone this repository. Navigate into the repository folder, change the details in the configuration file to match the system set up and run 'python operationalcode.py'.
