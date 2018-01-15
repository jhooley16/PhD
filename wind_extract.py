import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap, shiftgrid

for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
    wind_file = '/Users/jmh2g09/Documents/PhD/Data/Wind/ERA_Interim_' + year + '_10m_wind.nc'
    
    nc = Dataset(wind_file, 'r')
    lat = nc.variables['latitude'][:] # shape = 81
    pre_lon = nc.variables['longitude'][:] # shape = 720
    
    pre_u10 = nc.variables['u10'][:] # (12, 81, 720)
    pre_v10 = nc.variables['v10'][:] # (12, 81, 720)
    nc.close()
    
    lon = np.append(pre_lon, -pre_lon[0])
    u10 = np.dstack((pre_u10, pre_u10[:, :, 0]))
    v10 = np.dstack((pre_v10, pre_v10[:, :, 0]))
    
    for imonth in range(12):
        months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
        month = months[imonth]
        
        nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Wind/INPUT_u.nc', 'w')
        nc.createDimension('lat', np.size(lat))
        nc.createDimension('lon', np.size(lon))
        latitude = nc.createVariable('lat', float, ('lat',))
        longitude = nc.createVariable('lon', float, ('lon',))
        wind_save = nc.createVariable('wind', float, ('lat','lon',))
        latitude[:] = lat
        longitude[:] = lon
        wind_save[:] = u10[imonth, :, :]
        nc.close()
        
        nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Wind/INPUT_v.nc', 'w')
        nc.createDimension('lat', np.size(lat))
        nc.createDimension('lon', np.size(lon))
        latitude = nc.createVariable('lat', float, ('lat',))
        longitude = nc.createVariable('lon', float, ('lon',))
        wind_save = nc.createVariable('wind', float, ('lat','lon',))
        latitude[:] = lat
        longitude[:] = lon
        wind_save[:] = v10[imonth, :, :]
        nc.close()
        
        os.system('gmt grdsample /Users/jmh2g09/Documents/PhD/Data/Wind/INPUT_u.nc \
            -G/Users/jmh2g09/Documents/PhD/Data/Wind/OUTPUT_u.nc -I1.0/0.5 -R-180/180/-79/-50')
        os.system('rm /Users/jmh2g09/Documents/PhD/Data/Wind/INPUT_u.nc')
        
        os.system('gmt grdsample /Users/jmh2g09/Documents/PhD/Data/Wind/INPUT_v.nc \
            -G/Users/jmh2g09/Documents/PhD/Data/Wind/OUTPUT_v.nc -I1.0/0.5 -R-180/180/-79/-50')
        os.system('rm /Users/jmh2g09/Documents/PhD/Data/Wind/INPUT_v.nc')

        nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Wind/OUTPUT_u.nc', 'r')
        lat_wind = nc.variables['y'][:]
        lon_wind = nc.variables['x'][:]
        u_wind_data = nc.variables['z'][:]
        nc.close()
        os.system('rm /Users/jmh2g09/Documents/PhD/Data/Wind/OUTPUT_u.nc')
        
        nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Wind/OUTPUT_v.nc', 'r')
        v_wind_data = nc.variables['z'][:]
        nc.close()
        os.system('rm /Users/jmh2g09/Documents/PhD/Data/Wind/OUTPUT_v.nc')
        
        nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Wind/' + year + month + '_wind.nc', 'w')
        nc.createDimension('lat', np.size(lat_wind))
        nc.createDimension('lon', np.size(lon_wind))
        latitudes = nc.createVariable('latitude', float, ('lat',))
        longitudes = nc.createVariable('longitude', float, ('lon',))
        u_wind = nc.createVariable('u_wind', float, ('lat','lon',))
        v_wind = nc.createVariable('v_wind', float, ('lat','lon',))
        latitudes.long_name = 'latitude'
        latitudes.standard_name = 'latitude'
        latitudes.units = 'degrees_north'
        longitudes.long_name = 'longitude'
        longitudes.standard_name = 'longitude'
        longitudes.units = 'degrees_east'
        u_wind.standard_name = 'u_wind'
        u_wind.long_name = '10m_u_wind_component'
        u_wind.units = 'meters per second'
        v_wind.standard_name = 'v_wind'
        v_wind.long_name = '10m_v_wind_component'
        v_wind.units = 'meters per second'
        latitudes[:] = lat_wind
        longitudes[:] = lon_wind
        u_wind[:] = u_wind_data
        v_wind[:] = v_wind_data
        nc.close()

        pl.figure()
        pl.clf()
        m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
        m.drawmapboundary()
        m.drawcoastlines(zorder=10)
        m.fillcontinents(zorder=10)
        m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        grid_lats, grid_lons = np.meshgrid(lat_wind, lon_wind)
        stereo_x, stereo_y = m(grid_lons, grid_lats)
        m.pcolor(stereo_x, stereo_y, np.transpose(u_wind_data))
        m.colorbar()
        #pl.clim(5, -5)
        pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Wind/Figures/'+ year + month + '_10m_u_wind.png', 
            format='png')
        pl.close()
        
        pl.figure()
        pl.clf()
        m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
        m.drawmapboundary()
        m.drawcoastlines(zorder=10)
        m.fillcontinents(zorder=10)
        m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        grid_lats, grid_lons = np.meshgrid(lat_wind, lon_wind)
        stereo_x, stereo_y = m(grid_lons, grid_lats)
        m.pcolor(stereo_x, stereo_y, np.transpose(v_wind_data))
        m.colorbar()
        #pl.clim(5, -5)
        pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Wind/Figures/'+ year + month + '_10m_v_wind.png', 
            format='png')
        pl.close()