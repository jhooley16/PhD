import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

year = input('What Year? (xxxx): ')
lat_resolution = input('What latitude resolution?: ')
lon_resolution = input ('What longitude resolution?: ')

# Look at the raw data files
os.chdir('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + '/DOT_track')

# Cycle through each raw file
for file in os.listdir():
    if file[-15:] == '09_DOT_track.nc':
        print('Gridding: ' + file)
        nc = Dataset(file, 'r')
    
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        dot = nc.variables['dynamic_ocean_topography'][:]
        ice_conc = nc.variables['sea_ice_concentration'][:]
    
        nc.close()
    
        # Grid dynamic ocean topography to a 0.5-degree grid
        # grid05(data, lon, lat, lat_resolution, lon_resolution)
        data = funct.grid05(dot, lon, lat, float(lat_resolution), float(lon_resolution))
        grid_dot = data['Grid']
        grid_lon = data['Lon']
        grid_lat = data['Lat']
        
        data = funct.grid05(ice_conc, lon, lat,float(lat_resolution), float(lon_resolution))
        grid_ice = data['Grid']

        month = file[4:6]
        pl.figure()
        pl.clf()
        m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
        m.drawmapboundary()
        m.drawcoastlines(zorder=10)
        m.fillcontinents(zorder=10)
        m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        
        grid_lats, grid_lons = np.meshgrid(grid_lat, grid_lon)
        stereo_x, stereo_y = m(grid_lons, grid_lats)
        
        m.pcolor(stereo_x, stereo_y, grid_dot)
        m.colorbar()
        pl.clim(np.mean(grid_dot) + 3*np.std(grid_dot), np.mean(grid_dot) - 3*np.std(grid_dot))
        m.contour(stereo_x, stereo_y, grid_ice, [60,])
        pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/'+ year +'/DOT/Figures/' 
            + year + month + '_DOT_gridded_' + lon_resolution + 'x' 
            + lat_resolution + '.png', format='png', dpi=1200)
        pl.close()
    
        # Put the data in a .nc file in /Users/jmh2g09/Documents/PhD/Data/Gridded     
        nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year + '/DOT' + year + month + '_DOT.nc', 'w', format='NETCDF4_CLASSIC')

        nc.createDimension('lat', np.size(grid_lat))
        nc.createDimension('lon', np.size(grid_lon))

        latitudes = nc.createVariable('Latitude', float, ('lat',))
        longitudes = nc.createVariable('Longitude', float, ('lon',))
        gridded_dot = nc.createVariable('dynamic_ocean_topography', float, ('lon','lat'))
        gridded_ice = nc.createVariable('sea_ice_concentration', float, ('lon','lat'))

        latitudes.long_name = 'latitude'
        latitudes.standard_name = 'latitude'
        latitudes.units = 'degrees_north'
        longitudes.long_name = 'longitude'
        longitudes.standard_name = 'longitude'
        longitudes.units = 'degrees_east'
        gridded_dot.long_name = 'dynamic_ocean_topography'
        gridded_dot.standard_name = 'sea_surface_height_above_GOCO05s_geoid'
        gridded_dot.units = 'm'
        gridded_ice.long_name = 'sea_ice_concentration'
        gridded_ice.standard_name = 'sea_ice_concentration'
        gridded_ice.units = '%'

        latitudes[:] = grid_lat
        longitudes[:] = grid_lon
        gridded_dot[:] = grid_dot
        gridded_ice[:] = grid_ice

        nc.close()
        print('Complete!')