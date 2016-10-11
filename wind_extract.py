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

lon = np.append(pre_lon, -pre_lon[0])

u10 = np.dstack((pre_u10, pre_u10[:, :, 0]))
v10 = np.dstack((pre_v10, pre_v10[:, :, 0]))

# Cycle through the data, extract each month separately
#new_time = 1900 + time / (24 * 365.25)
count = 0
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if os.path.isdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + year + month + '_elev'):
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Wind/' + year + month + '_wind.nc', 'w')
    
            nc.createDimension('lat', np.size(lat))
            nc.createDimension('lon', np.size(lon))

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

            latitudes[:] = lat
            longitudes[:] = lon
            u_wind[:] = np.squeeze(u10[count, :, :])
            v_wind[:] = np.squeeze(v10[count, :, :])

            nc.close()
            
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            #m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        
            grid_lons, grid_lats = np.meshgrid(lon, lat)
            stereo_x, stereo_y = m(grid_lons, grid_lats)
            
            m.pcolor(stereo_x, stereo_y, np.squeeze(u10[count, :, :]))
            
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
        
            grid_lons, grid_lats = np.meshgrid(lon, lat)
            stereo_x, stereo_y = m(grid_lons, grid_lats)
            
            m.pcolor(stereo_x, stereo_y, np.squeeze(v10[count, :, :]))
            
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

            ugrid, newlons = shiftgrid(180.,np.squeeze(np.flipud(u10[count, :, :])),lon,start=False)
            vgrid, newlons = shiftgrid(180.,np.squeeze(np.flipud(v10[count, :, :])),lon,start=False)
            
            uproj,vproj,xx,yy = m.transform_vector(ugrid,vgrid,newlons,np.flipud(lat),70, 70,returnxy=True,masked=True)

            m.quiver(xx, yy,  -uproj, -vproj)
            
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Wind/Figures/'+ year + month + '_wind_quiver.png', 
                format='png')
            pl.close()
            
            count += 1