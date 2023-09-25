# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:16:42 2019

@author: john

"""

### This script combined with GEE_pull_functions_rivers extracts 
### Landsat surface reflectance over a channel water mask and
### caculates reach medians over all bands.
### Inputs: centerlines (SWORD), reach polygons (made in previous script),
###         and a DEM (MERIT)

### load libraries
import time
import ee
import os
#import numpy
#import pandas
#import feather

ee.Initialize()

# Note: This script uses python 3 which has syntax differences from python 2.7

# Source necessary functions.
exec(open("D:/Dropbox/projects/globalTSS/GEE_pull_functions_rivers.py").read())    



# load reach polygons with same ID and footprint as NHD reach, but wider.
poly = ee.FeatureCollection("users/johngardner87/reach_polygons_amazon");

# load amazon basin
hb3 = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_3")

def simplify500(feature): 
    return feature.simplify(maxError= 500)

def simplify100(feature): 
    return feature.simplify(maxError =100)

hb3_amazon = hb3.filter(ee.Filter.eq('PFAF_ID', 622))

centerline = ee.FeatureCollection("users/johngardner87/SWORD_SA")\
    .filterBounds(hb3_amazon.geometry())\
    .filter(ee.Filter.lt('type', 4))

centerline_simple = centerline.map(simplify100)

wrs = ee.FeatureCollection("users/johngardner87/WRS2_descending")\
    .filterBounds(hb3_amazon)\
    .filterBounds(centerline_simple)\
    .filterMetadata('MODE', 'equals', 'D')

pr = wrs.aggregate_array('PR').getInfo()

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
    .filter(ee.Filter.lt('CLOUD_COVER', 60))\
    .filterBounds(hb3_amazon)

## Set up a counter and a list to keep track of what's been done already
#counter = 0
#done = []    

#pr = [i for i in pr if i not in done]

# load reach polygons with same ID and footprint as NHD reach, but wider.

# make a folder in your google drive manually to output data
dlDir = 'D:/GoogleDrive/AmazonSR' 
filesDown = os.listdir(dlDir)  # -->
filesDown = [str(i.replace(".csv", "")) for i in filesDown]

pr  = [i for i in pr if i not in filesDown]

print(len(pr))

#for tiles in pr:
for x in range(0,len(pr)):

    print(pr[x])
    
    tile = wrs.filterMetadata('PR', 'equals', pr[x])
    
    reach_poly = poly.filterBounds(tile.geometry())
    reach = centerline.filterBounds(tile.geometry())
    
    stack = ls.filterBounds(tile.geometry().centroid())
    
    out = stack.map(RefPull).flatten() 
    
    dataOut = ee.batch.Export.table.toDrive(collection = out,\
                                            description = str(pr[x]),\
                                            folder = 'AmazonSR',\
                                            fileFormat = 'csv',\
                                            selectors = ['Aerosol','Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2', 'TIR1','TIR2', 'sd_NirSD', 'pixel_qa', 'clouds', 'dswe', 'hillShadow', 'pCount_dswe3','pCount_dswe1', 'pCount_shadow','system:index', 'CLOUD_COVER', 'SOLAR_ZENITH_ANGLE', 'SOLAR_AZIMUTH_ANGLE', 'reach_ID'])
    #Check how many existing tasks are running and take a break if it's >15  
    maximum_no_of_tasks(15, 60)
  # Send next task.
    dataOut.start()
    print('done')

# Make sure all Earth engine tasks are completed prior to moving on.  
maximum_no_of_tasks(1,300)
print('done')

