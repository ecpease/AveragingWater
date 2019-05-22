import pandas as pd
pd.options.display.max_columns = 10
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon
import os
import rasterio
from earthpy.clip import _clip_line_poly
import numpy as np
import skspatial
import matplotlib.pyplot as plt

# cubic feet to acre feet
def cft2acreft(vol):
	return vol * 2.29569e-5

def USGS_vols_to_dict(sites, start=2000, end=2010): # default dates, but these are changed in the other script
	data_dir = os.path.join('..','Data','SurfaceWater','USGS_Data')

	# Set naming convention (yyyymm)
	vol_dict = {}
	for y in range(start,end+1):
		for m in range(12):
			vol_dict[str(y)+'{0:02}'.format(m+1)] = {}


	for site in sites:
		df = pd.read_csv(os.path.join(data_dir,site+'.csv'))
		if len(df) > 1: # calls gages that have data 
			columns = df.columns
			usecol = [item for item in columns if item.endswith('00060_00003')] # discharge column (mean = 00003, max=1, min=2)
			if len(usecol) == 1: # 1 because sometimes csvs have multiple 00060 columns for some reason
				df = df[['datetime', usecol[0]]] # df is only datetime and discharge columns
				df['datetime'] = pd.to_datetime(df['datetime']) 
				df['cubic_ft'] = df[usecol] * (60*60*24)
				df.set_index('datetime', inplace=True) # reset index after resampling
				df = pd.DataFrame(df.resample('D')['cubic_ft'].mean()) # mean daily
				# df['cubic_ft'] = df['cubic_ft'].interpolate()
				df = pd.DataFrame(df.resample('M')['cubic_ft'].sum()) # sum for the month
				df['acre_ft'] = cft2acreft(df['cubic_ft']) # convert to acrefeet
				df.reset_index(inplace=True) # no more datetime index

				# set column format to yyyymm for each site through time
				df['yyyymm'] = df.apply(lambda i: str(i['datetime'].year) + '{0:02}'.format(i['datetime'].month), axis=1)
				df.set_index('yyyymm', inplace=True)
				available_months = df.index.tolist() # list of the months where there is available data (sometimes gages break)
				for yyyymm in vol_dict.keys():
					if yyyymm in available_months:
						vol_dict[yyyymm][site] = df.loc[yyyymm]['acre_ft'] # find the date key, find site key, then adds the acre-ft value to the dictionary
	return vol_dict


if __name__ == '__main__':
	sites = ['8052700','8095300','8164370','8449100','8052745']
	d = USGS_vols_to_dict(sites, 2000, 2015)
	print(d['201505']['8052700'])





