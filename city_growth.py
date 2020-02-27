#!/usr/bin/python

import requests 
import urllib
import json
import pandas as pd
import numpy as np
import geopandas as gpd
import shapely
from shapely.geometry import Point, Polygon, MultiPoint
from shapely import wkt
import shapely.speedups
from shapely.ops import transform, nearest_points
import plotly.express as px
import plotly.graph_objects as go
import os
import sys
import gdal
import rasterio
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
import glob
from functools import partial
import pyproj
import osmnx as ox
from IPython.display import Image
import make_fishnet
import cbd_osm
import geocoder
from pandana.loaders import osm
import pandana
import pylab as pl
ox.config(log_console=True, use_cache=True)


===================================================================================================================================


def get_city_proj_crs(to_crs):
    
    """
    Function to indentify local projection for cities dynamically
    
    Input:
    to_crs : name of city / country; epsg if known
    
    Returns:
    Local epsg (in string)
    """

    if isinstance(to_crs, int):
        to_crs = to_crs
    elif isinstance(to_crs, str):
        city, country = to_crs.split(',')
        url = "http://epsg.io/?q={}&format=json&trans=1&callback=jsonpFunction".format(city)
        r = requests.get(url)
        if r.status_code == 200:
            js = json.loads(r.text[14:-1])
            
            if js['number_result'] != 0:
                lis = []
                for i in js['results']:
                    res = i
                    if (res['unit'] == 'metre') and (res['accuracy'] == 1.0):
                        lis.append(res['code'])
                if len(lis) == 0:
                    for i in js['results']:
                        res = i
                        if res['unit'] == 'metre':
                            lis.append(res['code'])
                    return lis[0]
                else:
                    return lis[0]
   
            else:
                url = "http://epsg.io/?q={}&format=json&trans=1&callback=jsonpFunction".format(country)
                r = requests.get(url)
                if r.status_code == 200:
                    js = json.loads(r.text[14:-1])

                    if js['number_result'] != 0:
                        lis = []
                        for i in js['results']:
                            res = i
                            if (res['unit'] == 'metre') and (res['accuracy'] == 1.0):
                                lis.append(res['code'])
                        if len(lis) == 0:
                            for i in js['results']:
                                res = i
                                if res['unit'] == 'metre':
                                    lis.append(res['code'])
                            return lis[0]
                        else:
                            return lis[0]  
                        
===================================================================================================================================


def reproject_geom(geom, in_epsg=None, out_epsg=None):
    
    """
    Function to transform shapely geometry from in_epsg to out_epsg
    
    Input:
    geom: (required) : Shapely object (point/polygon/polyline)
    in_epsg : Input EPSG for geometry; defaults to 4326
    out_epsg: Output local EPSG for geometry; defaults to 3857
    
    Returns:
    Transformed Shapely object
    """
    
    if not in_epsg:
        in_epsg = 'epsg:3857' 
    if not out_epsg:
        out_epsg = 'epsg:4326'
    
    project = partial(
        pyproj.transform,
        pyproj.Proj(init= in_epsg), # source coordinate system
        pyproj.Proj(init= out_epsg)) # destination coordinate system

    geom = transform(project, geom)
    
    return geom


===================================================================================================================================


def polygonize_raster(ras_path, shp_path, string):
    """
    Function to polygonize a raster based on the pixel size of base raster.
    
    Inputs:
    ras_path: path to base raster location that is to be polygonized
    shp_path: path to where the shapefile will be saved
    string: name of the city
    
    Returns:
    Geodataframe with polygons equivalent to raster pixels.
    """
    
    print("Polygonizing Raster!!")
    
    import polygonize as pz
    
    path = os.getcwd()
       
    outSHPfn = path+"\\shapefiles\\{}".format(shp_path)
    lat, lon = pz.main(ras_path,outSHPfn)

    sh = gpd.read_file(path+"\\shapefiles\\{}".format(shp_path))
    sh.crs = {'init':'epsg:4326'}

    rio = rasterio.open(ras_path)
    shp_arr = np.array(sh.geometry).reshape(rio.shape[0], rio.shape[1])
    
    ### The following code is creating a 2x2 point window in a 2D array to use the four points of pixel and creates a polygon
    
    pols = []
    for row in range(shp_arr.shape[0]-1):
        for col in range(shp_arr.shape[1]-1):
            pols.append(shapely.geometry.box(shp_arr[row+1][col].x, shp_arr[row+1][col].y, shp_arr[row][col+1].x, shp_arr[row][col+1].y ))

    gdf = gpd.GeoDataFrame()
    gdf['ID'] = [i for i in range(len(pols))]
    gdf['geometry'] = pols
    gdf.set_geometry('geometry', inplace=True)
    gdf.crs = {'init':'epsg:4326'}

    print("Populating avearge height!!")

    av_h = []
    for i in gdf.geometry:
        coords = getFeatures(convert_geom_to_shp(i, string))
        out_img, out_transform = mask(dataset=rio, shapes=coords, crop=True)
        av_h.append(out_img.sum()/out_img.shape[2])

    gdf['avg_height'] = av_h
    gdf['Lon'] = [i.centroid.x for i in gdf.geometry]
    gdf['Lat'] = [i.centroid.y for i in gdf.geometry]

    
    return gdf


===================================================================================================================================


def nodes_from_osm(gdf, ras_path):
    """
    Function to count the number of intersections in each polygon
    
    Input:
    gdf : A geodataframe (Ideally a grid file for the city)
    ras_path : raster path of base data
    
    Returns:
    Geodataframe with 'node_count' column added to it.
    """
    
    print("Populating nodes from OSM!!")
    
    rio = rasterio.open(ras_path)
    pol = shapely.geometry.box(rio.bounds[0], rio.bounds[1], rio.bounds[2], rio.bounds[3])
    
    G = ox.graph_from_polygon(pol)
    nodes = ox.graph_to_gdfs(G, nodes=True, edges=False)

    
    gdf.crs = {'init':'epsg:4326'}
    gdf_copy = gdf.copy()
    gdf_copy.to_crs(epsg=out_proj, inplace=True) 
    nodes.to_crs(epsg=out_proj, inplace=True)
    
    nodes['geomType'] = nodes.geom_type
    nodes = nodes[nodes['geomType'] != 'GeometryCollection']
    merged = gpd.sjoin( gdf_copy, nodes, how='left', op='intersects')
    grp = merged.groupby('ID').count()
    
    gdf['node_count'] = grp.Lon.tolist()
    
    return gdf

===================================================================================================================================


def get_distance_cbd(gdf):
    """
    Function to calculate distance to CBD from each pixel center
    
    Input:
    gdf: Geodataframe (ideally a grid file for the city)
    
    Returns:
    Geodataframe with 'dis_to_cbd' column added to it
    """
    
    print("Populating distance to CBD!!")
    
    gdf_copy = gdf.copy()
    gdf_copy.to_crs(epsg=out_proj, inplace=True) 
    cbd = gdf_copy[gdf_copy.node_count == gdf_copy.node_count.max()].geometry.iloc[0].centroid
    gdf['dis_to_cbd'] = [i.centroid.distance(cbd) for i in gdf_copy.geometry]
    
    return gdf


===================================================================================================================================


def get_distance_water(gdf, gdf_water):
    """
    Function to calculate distance to nearest water body
    
    Input:
    gdf: Geodataframe (Ideally a grid network of the city)
    gdf_water : Water mask geodataframe
    
    Returns:
    Geodataframe with column 'dist_to_water' added to it
    """
    
    print("Populating distance to water!!")
    gdf_land = gdf.copy()
    gdf_water.to_crs(epsg=out_proj, inplace=True)
    gdf_land.to_crs(epsg=out_proj, inplace=True)
    
    water_dist = []
    
    dest_water = MultiPoint([i.centroid for i in gdf_water.geometry])
    
    for i in gdf_land.index:
        if i % 1000 == 0:
            print("{0} of {1} rows processed" .format(i, len(gdf_land)))

        temp_cent = gdf_land.geometry[i].centroid
        nearest_geoms = nearest_points(temp_cent, dest_water)
        water_dist.append(nearest_geoms[0].distance(nearest_geoms[1]))
    
    gdf['dist_to_water'] = water_dist
    
    return gdf


===================================================================================================================================

def get_iso(city):
    """
    Function to get ISO-3 codes for countries
    
    Input:
    city: city name (Ideally in (city, country) format)
    
    Returns:
    ISO-3 code for the country
    """
    
    try:
        country = city.split(',')[1].strip().lower()
        if country == 'south korea':  ### incorrect output for South Korea's ISO code with API
            return 'kor'
        else:
            url = "https://restcountries.eu/rest/v2/name/{}".format(country)
            r = requests.get(url)
            return r.json()[0]['alpha3Code'].lower()
    except IndexError:
        url = "https://restcountries.eu/rest/v2/capital/{}".format(city)
        r = requests.get(url)
        return r.json()[0]['alpha3Code'].lower()
    
    
===================================================================================================================================


def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio accepts them"""
    import json
    return [json.loads(gdf.to_json())['features'][0]['geometry']]


===================================================================================================================================


def get_population(gdf, city, ras_path):
    """
    Function to estimate population for each pixel
    
    gdf: Geodataframe (Ideally a grid network of the city)
    city: city name (Ideally in (city, country) format)
    ras_path: path to population raster
    
    Returns:
    Geodataframe with aggregated population in 'population' column
    """
    
    print("Populating population from Facebook!!")
    
    iso = get_iso(city)
    
    fbras_lis = glob.glob(os.getcwd()+"\\*{}*tif".format(iso))
    
    rio = rasterio.open(ras_path)
    
    pol = shapely.geometry.box(rio.bounds[0], rio.bounds[1], rio.bounds[2], rio.bounds[3])
    
    ## Parsing through multiple rasters to check which one intersects with input data
    
    if len(fbras_lis)>1:
        for path_ in fbras_lis:
            ras = rasterio.open(path_)
            if pol.intersects(shapely.geometry.box(ras.bounds[0], ras.bounds[1], ras.bounds[2], ras.bounds[3])):
                pop_path = path_
    else:
        pop_path = fbras_lis[0]
        
        
    pop = rasterio.open(pop_path)
    
    fb_pop = []
    
    for i in gdf.index:
        _gdf = gdf[gdf.index == i]

        _coords = getFeatures(_gdf)

        _out_img, _out_transform = mask(dataset=pop, shapes=_coords, crop=True)

        outimg = np.nan_to_num(_out_img)
        outimg = outimg.reshape(outimg.shape[1], outimg.shape[2])

        fb_pop.append(outimg.sum())
        
    gdf['population'] = fb_pop
    
    return gdf


===================================================================================================================================



def get_subway_data(city):
    """
    Function seacrches for subway data for cities on either local files or queries data from OSM.
    
    city: city name (Ideally in (city, country) format)
    
    Returns:
    Geodataframe with subway stations as point features or None if OSM query contains no data
    """
    
    string = city.split(",")[0]
    
    str_ = string.replace(' ','')
     
    sub_path = os.getcwd()+"\Subway and growth\subway_census_v1\station_points"
    
    sub_data = gpd.read_file(sub_path+"\\subway_stations2.shp")
    
    if str_ in sub_data.CITY1.unique():
        city_data = sub_data[sub_data.CITY1 == str_]
        return pois
    else:
        bbox = geocoder.arcgis(city).geojson['features'][0]['bbox']
        amenities = ['subway', 'light_rail', 'metro', 'underground', 'monorail', 'tram']
        osm_tags = '"railway"~"{}"'.format('|'.join(amenities))
        try:
            pois = osm.node_query(bbox[1], bbox[0], bbox[3], bbox[2],tags=osm_tags) ##lat_min, lng_min, lat_max, lng_max
            return pois
    
        except RuntimeError as e:
            if e.args[0] == "OSM query results contain no data.":
                return None


===================================================================================================================================



def get_distance_subway(gdf, city):
    """
    Function computes distance to nearest metro station for each pixel
    
    gdf: Geodataframe (Ideally a grid network of the city)
    city : city name (Ideally in (city, country) format)
    
    Returns:
    Geodataframe with 'dis_to_subway' column added to it
    """
    
    print("Populating distance to Subway stations!!")
     
    if not gdf.crs:
        gdf.crs = {'init':'epsg:4326'}
        
    gdf_copy = gdf.copy()
    
    gdf_copy.to_crs(epsg=out_proj, inplace=True)
    
    data =  get_subway_data(city)       
    
    if data:
        metro_dist = []

        dest = MultiPoint([i for i in data.geometry])

        for i in gdf_copy.index:
            if i % 1000 == 0:
                print("{0} of {1} rows processed" .format(i, len(gdf_copy)))

            temp_cent = gdf_copy.geometry.iloc[i].centroid

            nearest_geoms = nearest_points(temp_cent, dest)
            metro_dist.append(nearest_geoms[0].distance(nearest_geoms[1]))
    else:
        metro_dist = [0 for i in range(len(gdf_copy))]
        
    gdf['dis_to_subway'] = metro_dist
    
    return gdf



===================================================================================================================================


def get_distance_hwy(gdf, city):
    
    print("Populating distance to highways!!")
        
    if not gdf.crs:
        gdf.crs = {'init':'epsg:4326'}
    
    gdf_copy = gdf.copy()
    
    gdf_copy.to_crs(epsg=out_proj, inplace=True)
    
    bbox = geocoder.arcgis("{}".format(city)).geojson['features'][0]['properties']['raw']['extent']
    
    highway = ['motorway', 'trunk', 'primary', 'secondary', 'tertiary', 'residential']
    osm_tags = '"highway"~"{}"'.format('|'.join(highway))
    try:
        highway_pois = osm.node_query(bbox['ymin'],bbox['xmin'],bbox['ymax'],bbox['xmax'],tags=osm_tags) ##lat_min, lng_min, lat_max, lng_max

        highway_pois = highway_pois[['lat', 'lon', 'highway', 'name']]
        ## Adding geometry to the dataset
        highway_pois['geometry'] = (list(zip(highway_pois.lon,highway_pois.lat)))
        highway_pois['geometry'] = highway_pois.geometry.apply(lambda x: Point(x))
        highway_pois = gpd.GeoDataFrame(highway_pois, geometry='geometry')
        highway_pois.crs = {'init':'epsg:4326'}

        highway_pois.to_crs(epsg=out_proj, inplace=True)

        hwy_dist = []

        dest_hwy = MultiPoint([i for i in highway_pois.geometry])

        for i in gdf_copy.index:
            if i % 1000 == 0:
                print("{0} of {1} rows processed" .format(i, len(gdf_copy)))

            temp_cent = gdf_copy.geometry.iloc[i].centroid

            nearest_geoms = nearest_points(temp_cent, dest_hwy)
            hwy_dist.append(nearest_geoms[0].distance(nearest_geoms[1]))
    except RuntimeError as e:
            if e.args[0] == "OSM query results contain no data.":
                 hwy_dist = [0 for i in range(len(gdf))]
            else:
                print(e)
    
    gdf['dis_to_hwy'] = hwy_dist
    
    return gdf


===================================================================================================================================


def get_built_year(gdf, string):
    """
    Function to identify minimum and maximum of year for each of the pixels
    
    Input:
    gdf: Geodataframe (Ideally a grid network of the city)
    string : city name 
    
    Returns:
    Geodataframe with 'Yr_min_built' and 'Yr_max_built' columns added to it
    """
    
    print("Populating builup year!!")
    
    ras_path = glob.glob(os.getcwd()+"\Data\DLR Data\\*{}*_WSFevolution.tif".format(string))[0]
    
    ras = rasterio.open(ras_path)
    
    min_, max_ = [], []
    for i in gdf.index:
        if i % 1000 == 0:
            print("{0} of {1} rows processed" .format(i, len(gdf)))

        gdf_ = convert_geom_to_shp(gdf.geometry[i], 'Auckland')
        coords = getFeatures(gdf_)
        out_img, out_transform = mask(dataset=ras, shapes=coords, crop=True)
        un = np.unique(out_img)
        if (un[0] == 0) and (len(un) > 1):
            min_.append(un[1])
        else:
            min_.append(un[0])

        max_.append(un[-1])
            
    gdf['Yr_min_built'] = min_
    gdf['Yr_max_built'] = max_
    
    return gdf


===================================================================================================================================


def get_time(t1, t2):
    """
    Function to return difference between two timestamps
    
    t1: initial time
    t2: final time
    """
    diff = t2 - t1
    
    c = round(diff.total_seconds() / 60, 2)
    
    return c


===================================================================================================================================


def main(city):
    """
    Main function to compute all the metrics required for analysis of city growth
    
    Input:
    city: city name (Ideally in (city, country) format)
    
    Returns:
    None; Writes a shapefile at "pixel_level_files" location
    """
    
    import datetime
    global out_proj
    
    string = city.split(',')[0]
    path = os.getcwd()
    
    ras_path = glob.glob(os.getcwd()+"\Data\DLR Data\\*{}*WSF3D_AW3D30.tif".format(string))[0]
    shp_path = "{}_grid.shp".format(string)
    
    out_proj = get_city_proj_crs(city)
    
    t1 = datetime.datetime.now()
    
    gdf = polygonize_raster(ras_path, shp_path, string)

    gdf_water = gdf[gdf.avg_height == 0.0]
    gdf_water.reset_index(drop=True, inplace=True)
    
    gdf_land = gdf[gdf.avg_height != 0.0]
    gdf_land.reset_index(drop=True, inplace=True)
    
    gdf_water.crs = {'init':'epsg:4326'}
    
    gdf_land.crs = {'init':'epsg:4326'}
    
    gdf_land = nodes_from_osm(gdf_land, ras_path)

    gdf_land = get_distance_cbd(gdf_land)

    gdf_land = get_distance_water(gdf_land, gdf_water)

    gdf_land = get_population(gdf_land, city, ras_path)
    
    gdf_land = get_distance_subway(gdf_land, city)
    
    gdf_land = get_distance_hwy(gdf_land, city)
    
    gdf_land = get_built_year(gdf_land, string)
    
    print("Writing data to shapefile!!")
    
    gdf_land.to_file(path+"\\pixel_level_files\\{}_pixel_level_data.shp".format(string))
    
    t2 = datetime.datetime.now()
    
    print('Total time taken to run the analysis: {}'.format(get_time(t1, t2)))


===================================================================================================================================


if __name__ == "__main__":
    if len(sys.argv) == 2:
        main(sys.argv[1])



