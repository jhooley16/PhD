## Calculate the geoid height using EIGEN6c4

from netCDF4 import Dataset
import os
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
import numpy as np

# load the monthly along track data

for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if os.path.isfile('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + '/' + year + month + '_raw.nc'):
            print(year, month)

            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + '/' + year + month + '_raw.nc', 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            nc.close()
            
            # Construct a two column ascii table for the location of points
            input_file = open('/Users/jmh2g09/Documents/PhD/Data/Geoid/INPUT.txt', 'w')
            for i in range(0, len(lat)):
                print(lon[i], lat[i], file=input_file)
            input_file.close()

            # insert this into gmt function 'grdtrack' and interpolate the points
            os.system('gmt grdtrack /Users/jmh2g09/Documents/PhD/Data/Geoid/INPUT.txt -f0x,1y -G/Users/jmh2g09/Documents/PhD/Data/Geoid/EIGEN6c4/EIGEN6c4_GRID.nc > /Users/jmh2g09/Documents/PhD/Data/Geoid/OUTPUT.txt')

            os.system('rm /Users/jmh2g09/Documents/PhD/Data/Geoid/INPUT.txt')
            # open the output geoid data
            lon_new = []
            lat_new = []
            geoid = []
            f = open('/Users/jmh2g09/Documents/PhD/Data/Geoid/OUTPUT.txt', 'r')
            for line in f:
                columns = line.strip().split()
                lon_new.append(float(columns[0]))
                lat_new.append(float(columns[1]))
                geoid.append(float(columns[2]))
            f.close()
            os.system('rm /Users/jmh2g09/Documents/PhD/Data/Geoid/OUTPUT.txt')
            
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + '/' + year + month + '_geoid_track.nc', 'w')
            nc.createDimension('station', np.size(lat_new))
            latitudes = nc.createVariable('latitude', float, ('station',))
            longitudes = nc.createVariable('longitude', float, ('station',))
            geoid_height = nc.createVariable('geoid_height', float, ('station',))

            latitudes.long_name = 'latitude'
            latitudes.standard_name = 'latitude'
            latitudes.units = 'degrees_north'
            longitudes.long_name = 'longitude'
            longitudes.standard_name = 'longitude'
            longitudes.units = 'degrees_east'
            geoid_height.long_name = 'EIGEN6c4_geoid_height'
            geoid_height.standard_name = 'geoid_height_above_WGS84_ellipsoid'
            geoid_height.units = 'm'

            latitudes[:] = lat_new
            longitudes[:] = lon_new
            geoid_height[:] = geoid

            nc.close()