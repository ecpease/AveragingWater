"""
Steps:

1. Load packages (some you'll need to get from Github or find alternatives)
2. Functions to call later (convert polyline to points, calculate volume)
3. Load shapefiles (rivers and gages)
4. Load raster and create shapefile that contains just the extent of the raster tif
5. Clip TCEQ river segments to bounding box of raster
6. Merge TCEQ river segments into one single feature and clip USGS gages and create USGS site number object
7. Create new raster the extent of the study area (same nrow, ncol), set all values to zero and set nan value
   (here, it is -1000 ... we can change this)
8. Set dates to iterate through and naming convention (yyyymm)
9. Call get_USGS_data.py script to run function getting all volumes for USGS gages
10. Create points from polylines to overlay in raster 

"""

import pandas as pd
pd.options.display.max_columns = 20
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon
import os
import rasterio
from earthpy.clip import clip_shp, _clip_line_poly # this is a good clipping library
import numpy as np
import skspatial # from Ross K. github
import matplotlib.pyplot as plt
import get_USGS_data # other Python script

# Convert lines to points
def polyline_to_pts(line_geom):
	if type(line_geom) == LineString:
		x, y = line_geom.coords.xy
	else: # If multistring
		x, y = [], []
		for line in line_geom:
			xi, yi = line.coords.xy
			x += xi
			y += yi 
	return x,y

# Calculate volume for each USGS site
def get_vol(vol_dict, yyyymm, site):
	if site in vol_dict[yyyymm].keys():
		val = vol_dict[yyyymm][site]
	else:
		val = np.nan
	return val


def main():

	# Load data
	TCEQ_lines = gpd.read_file(os.path.join("..", "Data", "TCEQ_classified.shp")) # streams
	crs = TCEQ_lines.crs
	usgs_sites = gpd.read_file(os.path.join('..','Data','usgs_tceq_merged.shp')) # usgs gages

	# Make extent of master raster into geodataframe for clipping
	r = rasterio.open(os.path.join('..','Data','hagm_grid','MasterRaster_Houston.tif'))
	transform = r.transform # raster info on size, how its made, etc.
	bounds = r.bounds # corners
	pts = [(bounds.left, bounds.top), (bounds.right, bounds.top), (bounds.right, bounds.bottom), (bounds.left, bounds.bottom)]
	pts = Polygon(pts) # make the points a polygon shp of raster boundary
	extent_gdf = gpd.GeoDataFrame({'geometry':[pts]})
	extent_gdf.crs = crs # set the new gdf crs equal to the tceq crs
	ulc, lrc = (bounds.left, bounds.top), (bounds.right, bounds.bottom) # upper left corner, lower right corner

	# Clip the TCEQ streams to the bounding box gdf of Houston
	TCEQ_lines = TCEQ_lines.explode() # get rid of multigeometry
	TCEQ_lines.reset_index(inplace=True, drop=True)
	TCEQ_lines = clip_shp(TCEQ_lines, extent_gdf)
	print(TCEQ_lines.head())
	print(TCEQ_lines.columns)
	print(len(TCEQ_lines))
	
	# Make TCEQ clipped shapefile into one big shapefile based on Basin Name - not individual line segments
	# Clip USGS sites shapefile to the extent of the geodataframe
	# Assign 'sites' object as each unique site number
	TCEQ_lines = TCEQ_lines.dissolve('BASIN_NAME')
	print(len(TCEQ_lines))
	usgs_sites = clip_shp(usgs_sites, extent_gdf) # earthpy function
	sites = usgs_sites['site_no'].unique()

	# Create a raster the same size/shape as original but with all zeros
	# Then if values are less than -1000, set them to NaN
	array = r.read(1)
	nrow, ncol = array.shape
	active_mask = np.zeros((nrow, ncol))
	active_mask[array <= -1000.] = np.nan # -1000 is set as nan value - change this if not the case

	# Iterate through years and months to set naming convention (year yyyy month mm)
	dates = []
	start, end = 1989,2009
	for y in range(start, end+1):
		for m in range(12):
			dates.append(str(y)+'{0:02}'.format(m+1)) 

	# Get volumes for each site iterated through time (89-09) from other script
	vol_dict = get_USGS_data.USGS_vols_to_dict(sites, start=start, end=end)	
	for yyyymm in vol_dict.keys():
		print(yyyymm)
		final_array = np.zeros((nrow, ncol)) # set a raster to zeros initially
		for itr, dfrow in TCEQ_lines.iterrows():
			x, y = polyline_to_pts(dfrow['geometry']) # get x,y from TCEQ clipped shp and convert to points
			temp_df = pd.DataFrame({'UTMx':x,'UTMy':y}) # create dataframe using x and y values of the TCEQ clipped shp
			temp_df['geometry'] = temp_df.apply(lambda i: Point(i['UTMx'], i['UTMy']), axis=1)
			temp_df = gpd.GeoDataFrame(temp_df, geometry='geometry')
			temp_df.crs = crs

			# Spatial join on the sites using geometries to get the site number
			temp_df = gpd.sjoin(temp_df, usgs_sites[['site_no','geometry']], how='left')
			temp_df['acre_ft'] = temp_df.apply(lambda i: get_vol(vol_dict ,yyyymm, i['site_no']), axis=1)
			# temp_df['acre_ft'] = temp_df['acre_ft'].interpolate()

			if not pd.isnull(temp_df['acre_ft']).all():
				ml = skspatial.interp2d(temp_df,'acre_ft', res = 1609.35, ulc = ulc, lrc = lrc) # set up skpatial model object
				hgrid = ml.points_to_grid()
				hgrid[np.isnan(hgrid)] = 0 # look at this more! (check out type)
				final_array += np.array(hgrid)

		# Set all nan or 0 values to nan in the raster
		final_array[np.isnan(active_mask)] = np.nan
		final_array[final_array == 0] = np.nan # sometimes there are relevant zeros - check this out more! 

		# Export dataset (monthly) to raster
		if not os.path.exists(os.path.join('..','output_rasters')): os.mkdir(os.path.join('..','output_rasters'))
		new_dataset = rasterio.open(os.path.join('..','output_rasters',f'{yyyymm}.tif'), 'w', driver='GTiff',
                            height = final_array.shape[0], width = final_array.shape[1],
                            count=1, dtype=str(final_array.dtype),
                            crs=crs,
                            transform=transform, nodata= np.nan)
		new_dataset.write(final_array, 1)
		new_dataset.close()
	print(temp_df.head())

if __name__ == '__main__':
	main()
