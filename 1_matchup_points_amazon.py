# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:16:42 2019

@author: john

"""

### This script combined with GEE_pull_functions_rivers extracts 
### Landsat surface reflectance at a point over a water mask and
### caculates reach medians over all bands

### load libraries
import time
import ee
import os
#import numpy
import pandas
#import feather

ee.Initialize()

# Note: This script uses python 3 which has syntax differences from python 2.7
# Source necessary functions.
exec(open("D:/Dropbox/projects/tss_amazon/GEE_matchup_point_functions.py").read())    

#### LOAD IN IN_SITU POINTS
points = ee.FeatureCollection("users/johngardner87/amazon_matchup_sites_v1")

wrs = ee.FeatureCollection("users/johngardner87/WRS2_descending")\
    .filterBounds(points)\
    .filterMetadata('MODE', 'equals', 'D')

pr = wrs.sort('PR').aggregate_array('PR').getInfo()

dummyBands = ee.Image(-99).rename('Null_CS')\
    .addBands(ee.Image(-99).rename('Null_TIR2'))

#
def addDummy(i):
    return i.addBands(dummyBands)

#
l8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
l7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR').map(addDummy)
l5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR').map(addDummy)

# Standardize band names between the various collections and aggregate 
# them into one image collection
bn8 = ['B1','B2','B3', 'B4', 'B5','B6','B7', 'B10','B11','pixel_qa']
bn57 = ['Null_CS', 'B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6','Null_TIR2', 'pixel_qa']
bns = ['Aerosol','Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2', 'TIR1','TIR2','pixel_qa']

ls5 = l5.select(bn57, bns)
ls7 = l7.select(bn57, bns)
ls8 = l8.select(bn8, bns)

# set cloud threshold. Should turn down to 50 since I filter to that later
ls = ee.ImageCollection(ls5.merge(ls7).merge(ls8))\
    .filter(ee.Filter.lt('CLOUD_COVER', 50))\
    .filterBounds(wrs)

## Set up a counter and a list to keep track of what's been done already
#counter = 0
#done = []    
#pr = [i for i in pr if i not in done]


# make a folder in your google drive manually to output data
#dlDir = 'D:/GoogleDrive/Ohio_SR_Points' 
#filesDown = os.listdir(dlDir)  # -->
#filesDown = [str(i.replace(".csv", "")) for i in filesDown]

#pr  = [i for i in pr if i not in filesDown]

print(len(pr))
print(pr)

radius = 500

#for tiles in pr:
for x in range(0,len(pr)):
    
    print(pr[x])
    
    tile = wrs.filterMetadata('PR', 'equals', pr[x])
    
    river = points.filterBounds(tile.geometry())\
        .map(Buff)
    
    stack = ls.filterBounds(tile.geometry().centroid())
    
    out = stack.map(RefPull).flatten() 
    
    #print(out.limit(1).getInfo())
    
    dataOut = ee.batch.Export.table.toDrive(collection = out,\
                                            description = str(pr[x]),\
                                            folder = 'matchup_amazon',\
                                            fileFormat = 'csv',\
                                            selectors = ['Aerosol','Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2', 'TIR1','TIR2', 'sd_NirSD', 'pixel_qa', 'clouds', 'dswe', 'hillShadow', 'pCount_dswe3','pCount_dswe1', 'pCount_shadow','system:index', 'time', 'CLOUD_COVER', 'SOLAR_ZENITH_ANGLE', 'SOLAR_AZIMUTH_ANGLE', 'siteID'])
    #Check how many existing tasks are running and take a break if it's >15  
    maximum_no_of_tasks(15, 60)
  # Send next task.
    dataOut.start()
    print('done')

# Make sure all Earth engine tasks are completed prior to moving on.  
maximum_no_of_tasks(1,300)
print('done')

