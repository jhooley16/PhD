import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap, shiftgrid

nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Wind/ERA_Interim_10m_UV_wind_201011-201602.nc', 'r')

lat = nc.variables['latitude'][:] # shape = 81
pre_lon = nc.variables['longitude'][:] # shape = 720
time = nc.variables['time'][:] # shape = 64 # hours since 1900-01-01 00:00:0.0 

pre_u10 = nc.variables['u10'][:] # (time, lat, lon)
pre_v10 = nc.variables['v10'][:] # (time, lat, lon)

nc.close()
print(pre_lon)
lon = np.append(pre_lon, -pre_lon[0])
print(lon)
pause
u10 = np.dstack((pre_u10, pre_u10[:, :, 0]))
v10 = np.dstack((pre_v10, pre_v10[:, :, 0]))

# Cycle through the data, extract each month separately
#new_time = 1900 + time / (24 * 365.25)
count = 0
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if os.path.isdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + year + month + '_elev'):
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Wind/INPUT_u.nc', 'w')
            nc.createDimension('lat', np.size(lat))
            nc.createDimension('lon', np.size(lon))
            latitude = nc.createVariable('lat', float, ('lat',))
            longitude = nc.createVariable('lon', float, ('lon',))
            wind_save = nc.createVariable('wind', float, ('lat','lon',))
            latitude[:] = lat
            longitude[:] = lon
            wind_save[:] = np.squeeze(u10[count, :, :])
            nc.close()
    
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Wind/INPUT_v.nc', 'w')
            nc.createDimension('lat', np.size(lat))
            nc.createDimension('lon', np.size(lon))
            latitude = nc.createVariable('lat', float, ('lat',))
            longitude = nc.createVariable('lon', float, ('lon',))
            wind_save = nc.createVariable('wind', float, ('lat','lon',))
            latitude[:] = lat
            longitude[:] = lon
            wind_save[:] = np.squeeze(v10[count, :, :])
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
    
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Wind/OUTPUT_v.nc', 'r')
            v_wind_data = nc.variables['z'][:]
            nc.close()
    
            os.system('rm /Users/jmh2g09/Documents/PhD/Data/Wind/OUTPUT_u.nc')
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
            #m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        
            grid_lons, grid_lats = np.meshgrid(lon_wind, lat_wind)
            stereo_x, stereo_y = m(grid_lons, grid_lats)
            
            m.pcolor(stereo_x, stereo_y, u_wind_data)
            
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Wind/Figures/'+ year + month + '_10m_u_wind.png', 
                format='png')
            pl.close()
            
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            #m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
            
            m.pcolor(stereo_x, stereo_y, v_wind_data)
            
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Wind/Figures/'+ year + month + '_10m_v_wind.png', 
                format='png')
            pl.close()
            
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            #m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])

            ugrid, newlons = shiftgrid(180.,u_wind_data,lon_wind,start=False)
            vgrid, newlons = shiftgrid(180.,v_wind_data,lon_wind,start=False)
            
            uproj,vproj,xx,yy = m.transform_vector(ugrid,vgrid,newlons,lat_wind,70, 70,returnxy=True,masked=True)

            m.quiver(xx, yy,  -uproj, -vproj)
            
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Wind/Figures/'+ year + month + '_wind_quiver.png', 
                format='png')
            pl.close()
            
            count += 1