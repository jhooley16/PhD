import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

monthly_average = np.full((12, 7, 59, 361), fill_value=np.NaN)
monthly_average_2 = np.full((12, 7, 59, 361), fill_value=np.NaN)
monthly_average_3 = np.full((12, 7, 59, 361), fill_value=np.NaN)
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year)
    
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:

        if os.path.isfile(year + month + '_DOT_ice_only.nc'):
            print(year, month)
            # Load the lat and lon data for each file
            nc = Dataset(year + month + '_DOT_ice_only.nc', 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            nc.close()

            # Filter the DOT with a 500km gaussian filter
            os.system('gmt grdfilter ' + year + month + '_DOT_ice_only.nc?"dynamic_ocean_topography_seasonal_offset" -D4 -Fg500 -Nr -f0y -f1x -GDOT_filt.nc')
            os.system('gmt grdfilter ' + year + month + '_DOT_ice_only.nc?"dynamic_ocean_topography_no_offset" -D4 -Fg500 -Nr -f0y -f1x -GDOT_2_filt.nc')
            os.system('gmt grdfilter ' + year + month + '_DOT_ice_only.nc?"dynamic_ocean_topography_constant_offset" -D4 -Fg500 -Nr -f0y -f1x -GDOT_3_filt.nc')
            
            # Open the filtered data
            nc = Dataset('DOT_filt.nc', 'r')
            dot_filt = np.array(nc.variables['z'][:])
            nc.close()
            
            nc = Dataset('DOT_2_filt.nc', 'r')
            dot_2_filt = np.array(nc.variables['z'][:])
            nc.close()
            
            nc = Dataset('DOT_3_filt.nc', 'r')
            dot_3_filt = np.array(nc.variables['z'][:])
            nc.close()

            os.system('gmt grdfilter ' + year + month + '_DOT.nc?"sea_ice_concentration" -D4 -Fg1000 -Ni -f0y -f1x -GICE_filt.nc')

            nc = Dataset('ICE_filt.nc', 'r')
            ice_filt = np.array(nc.variables['z'][:])
            nc.close()

            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        
            grid_lats, grid_lons = np.meshgrid(lat, lon)
            stereo_x, stereo_y = m(grid_lons, grid_lats)
            pl.rcParams['contour.negative_linestyle'] = 'solid'
            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_filt)), cmap='RdBu_r')
            m.colorbar()
            #m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_filt)), 15, colors='k')
            #pl.clim([0, -2.25])
            #pl.clim(np.mean(np.ma.masked_invalid(dot_filt)) + 3*np.std(np.ma.masked_invalid(dot_filt)), np.mean(np.ma.masked_invalid(dot_filt)) - 3*np.std(np.ma.masked_invalid(dot_filt)))
            #m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_filt)), [20,])
            pl.savefig('Figures/' + year + month + '_DOT_filt_ice_only.png', format='png', transparent=True, dpi=300)
            pl.close()

            nc = Dataset(year + month + '_DOT_filt_ice.nc', 'w')

            nc.createDimension('lat', np.size(lat))
            nc.createDimension('lon', np.size(lon))

            latitudes = nc.createVariable('latitude', float, ('lat',))
            longitudes = nc.createVariable('longitude', float, ('lon',))
            sea_surface_height = nc.createVariable('dynamic_ocean_topography_seasonal_offset', float, ('lat','lon'))
            sea_surface_height_2 = nc.createVariable('dynamic_ocean_topography_no_offset', float, ('lat','lon'))
            sea_surface_height_3 = nc.createVariable('dynamic_ocean_topography_constant_offset', float, ('lat','lon'))
            sea_ice_concentration = nc.createVariable('sea_ice_concentration', float, ('lat', 'lon'))
        
            latitudes[:] = lat
            longitudes[:] = lon
            sea_surface_height[:] = dot_filt
            sea_surface_height_2[:] = dot_2_filt
            sea_surface_height_3[:] = dot_3_filt
            sea_ice_concentration[:] = ice_filt
            nc.close()
            
            os.system('rm DOT_filt.nc DOT_2_filt.nc DOT_3_filt.nc ICE_filt.nc')
            
            monthly_average[int(month)-1, int(year)-2010, :, :] = dot_filt
            monthly_average_2[int(month)-1, int(year)-2010, :, :] = dot_2_filt
            monthly_average_3[int(month)-1, int(year)-2010, :, :] = dot_3_filt

month_average = np.mean(monthly_average[:, 1:-1, :, :], axis=1)
month_average_2 = np.mean(monthly_average_2[:, 1:-1, :, :], axis=1)
month_average_3 = np.mean(monthly_average_3[:, 1:-1, :, :], axis=1)
annual_mean = np.mean(month_average, axis=0)
annual_mean_2 = np.mean(month_average_2, axis=0)
annual_mean_3 = np.mean(month_average_3, axis=0)

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
grid_lats, grid_lons = np.meshgrid(lat, lon)
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(annual_mean)), cmap='RdBu_r')
m.colorbar()
pl.savefig('../Figures/mean_DOT_filt_ice_only.png', format='png', transparent=True, dpi=300)
pl.close()

nc = Dataset('../MDT_mean_ice_only.nc', 'w')
nc.createDimension('lat', np.size(lat))
nc.createDimension('lon', np.size(lon))
latitudes = nc.createVariable('latitude', float, ('lat',))
longitudes = nc.createVariable('longitude', float, ('lon',))
sea_surface_height = nc.createVariable('mean_dynamic_topography_seasonal_offset', float, ('lat','lon'))
sea_surface_height_2 = nc.createVariable('mean_dynamic_topography_no_offset', float, ('lat','lon'))
sea_surface_height_3 = nc.createVariable('mean_dynamic_topography_constant_offset', float, ('lat','lon'))
latitudes[:] = lat
longitudes[:] = lon
sea_surface_height[:] = annual_mean
sea_surface_height_2[:] = annual_mean_2
sea_surface_height_3[:] = annual_mean_3
nc.close()

for imnth in range(12):
    pl.figure()
    pl.clf()
    m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])

    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(month_average[imnth, :, :])), cmap='RdBu_r')
    m.colorbar()
    pl.savefig('../Figures/'+ str(imnth + 1) + '_DOT_filt_ice_only.png', format='png', transparent=True, dpi=300)
    pl.close()

    nc = Dataset('../' + str(imnth + 1) + '_mean_DOT_filt_ice_only.nc', 'w')
    nc.createDimension('lat', np.size(lat))
    nc.createDimension('lon', np.size(lon))
    latitudes = nc.createVariable('latitude', float, ('lat',))
    longitudes = nc.createVariable('longitude', float, ('lon',))
    sea_surface_height = nc.createVariable('dynamic_ocean_topography_seasonal_offset', float, ('lat','lon'))
    sea_surface_height_2 = nc.createVariable('dynamic_ocean_topography_no_offset', float, ('lat','lon'))
    sea_surface_height_3 = nc.createVariable('dynamic_ocean_topography_constant_offset', float, ('lat','lon'))
    
    latitudes[:] = lat
    longitudes[:] = lon
    sea_surface_height[:] = month_average[imnth, :, :]
    sea_surface_height_2[:] = month_average_2[imnth, :, :]
    sea_surface_height_3[:] = month_average_3[imnth, :, :]
    nc.close()