import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    lat_resolution = '0.5'    #input('What latitude resolution?: ')
    lon_resolution = '1.0'      #input ('What longitude resolution?: ')

    # Look at the raw data files
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year)

    # Cycle through each raw file
    for file in os.listdir():
        if file[-7:] == '_raw.nc':
            print('Gridding: ' + file)
            nc = Dataset(file, 'r')

            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            ssh = nc.variables['sea_surface_height'][:]
            ice_conc = nc.variables['ice_concentration'][:]

            nc.close()

            # Grid dynamic ocean topography to a 1x0.5-degree grid
            # grid05(data, lon, lat, lat_resolution, lon_resolution)
            data = funct.grid05(ssh, lon, lat, float(lat_resolution), float(lon_resolution))
            grid_ssh = data['Grid']
            grid_lon = data['Lon']
            grid_lat = data['Lat']

            data = funct.grid05(ice_conc, lon, lat, float(lat_resolution), float(lon_resolution))
            grid_ice = data['Grid']

            # Make the longitudes between 0 and 360
            grid_ssh, grid_lon = funct.lon_convert(grid_lon, grid_ssh)
            grid_ice, grid_lon = funct.lon_convert(grid_lon, grid_ice)
            #grid_lon[grid_lon < 0] += 361

            #grid_ssh = grid_ssh[np.argsort(grid_lon), :]
            #grid_ice = grid_ice[np.argsort(grid_lon), :]
            #grid_lon = grid_lon[np.argsort(grid_lon)]
        
            
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
        
            m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(grid_ssh), cmap='RdBu_r')
            m.colorbar()
            #pl.clim(0, -2)
            #pl.clim(np.mean(np.ma.masked_invalid(grid_dot)) - 3*np.std(np.ma.masked_invalid(grid_dot)), np.mean(np.ma.masked_invalid(grid_dot)) + 3*np.std(np.ma.masked_invalid(grid_dot)))
            m.contour(stereo_x, stereo_y, np.ma.masked_invalid(grid_ice), [40,])
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/SSH/'+ year +'/Figures/' 
                + year + month + '_SSH_gridded_' + lon_resolution + 'x' 
                + lat_resolution + '.png', format='png', transparent=True, dpi=300)
            pl.close()
    
            # Put the data in a .nc file in /Users/jmh2g09/Documents/PhD/Data/Gridded     
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/SSH/' + year + '/' + year + month + '_SSH.nc', 'w')

            nc.createDimension('lat', np.size(grid_lat))
            nc.createDimension('lon', np.size(grid_lon))

            latitudes = nc.createVariable('latitude', float, ('lat',))
            longitudes = nc.createVariable('longitude', float, ('lon',))
            gridded_ssh = nc.createVariable('sea_surface_height', float, ('lat','lon'))
            gridded_ice = nc.createVariable('sea_ice_concentration', float, ('lat','lon'))

            latitudes.long_name = 'latitude'
            latitudes.standard_name = 'latitude'
            latitudes.units = 'degrees_north'
            longitudes.long_name = 'longitude'
            longitudes.standard_name = 'longitude'
            longitudes.units = 'degrees_east'
            gridded_ssh.long_name = 'sea_surface_height'
            gridded_ssh.standard_name = 'sea_surface_height_above_WGS84_ellipsoid'
            gridded_ssh.units = 'm'
            gridded_ice.long_name = 'sea_ice_concentration'
            gridded_ice.standard_name = 'sea_ice_concentration'
            gridded_ice.units = '%'

            latitudes[:] = grid_lat
            longitudes[:] = grid_lon
            gridded_ssh[:] = np.transpose(grid_ssh)
            gridded_ice[:] = np.transpose(grid_ice)

            nc.close()
            print('Complete!')