import os
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as pl
import functions as funct
# 
# ## Load the AVISO MSS
# nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/AVISO/AVISO_MSS/mss_cnes_cls2015.nc', 'r')
# latitude = nc.variables['NbLatitudes'][:]
# longitude = nc.variables['NbLongitudes'][:]
# mss = np.array(np.transpose(nc.variables['mss'][:]))
# nc.close()
# print('LAT: ', np.shape(latitude))
# print('LON: ', np.shape(longitude))
# print(np.shape(mss))
# 
# # Force mask to NaN
# mss[abs(mss) > 200] = np.NaN
# 
# # Save into dummy .nc file for filtering/resampling
# nc_test = Dataset('test_mss.nc', 'w')
# nc_test.createDimension('lat', np.size(latitude))
# nc_test.createDimension('lon', np.size(longitude))
# latitudes = nc_test.createVariable('lat', float, ('lat',))
# longitudes = nc_test.createVariable('lon', float, ('lon',))
# mss_test = nc_test.createVariable('z', float, ('lat','lon'))
# latitudes[:] = latitude
# longitudes[:] = longitude
# mss_test[:] = mss
# nc_test.close()
# 
# # Filter and resample
# os.system('gmt grdsample test_mss.nc -Gtest_mss_subsampled.nc -I1.0/0.5 -R0/360/-79/-50 -fg')
# os.system('gmt grdfilter test_mss_subsampled.nc -D4 -Fg600 -Nr -fg -Gtest_mss_filt.nc')
# os.system('rm test_mss.nc test_mss_subsampled.nc')
# 
# # Open filtered and resampled mss data
# nc = Dataset('test_mss_filt.nc', 'r')
# latitude = nc.variables['lat'][:]
# longitude = nc.variables['lon'][:]
# mss_sampled = np.squeeze(nc.variables['z'][:])
# nc.close()
# os.system('rm test_mss_filt.nc')

anomaly_anomaly = np.full((59, 361, 72), fill_value = np.NaN)
ix = 0
for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        print(year, month)
        ## Load the AVISO altimetry ssha data
        aviso_file = '/Users/jmh2g09/Documents/PhD/Data/AVISO/AVISO_Data/dt_global_allsat_msla_h_y' + year + '_m' + month + '.nc' 
    
        nc = Dataset(aviso_file, 'r')
        latitude = nc.variables['lat'][:]
        longitude = nc.variables['lon'][:]
        ssha = np.array(np.squeeze(nc.variables['sla'][:]))
        nc.close()
        
        # Force mask to NaN
        ssha[abs(ssha) > 100] = np.NaN
        
        # Save into a dummy .nc file for filtering and subsampling
        nc_test = Dataset('aviso.nc', 'w')
        nc_test.createDimension('lat', np.size(latitude))
        nc_test.createDimension('lon', np.size(longitude))
        latitudes = nc_test.createVariable('lat', float, ('lat',))
        longitudes = nc_test.createVariable('lon', float, ('lon',))
        ssha_test = nc_test.createVariable('ssha', float, ('lat','lon'))
        latitudes[:] = latitude
        longitudes[:] = longitude
        ssha_test[:] = ssha
        nc_test.close()
    
        # Filter and subsample
        os.system('gmt grdsample aviso.nc -Gaviso_subsampled.nc -I1.0/0.5 -R0/360/-79/-50 -fg')
        os.system('gmt grdfilter aviso_subsampled.nc -D4 -Fg600 -Nr -fg -Gaviso_filt.nc')
        os.system('rm aviso.nc aviso_subsampled.nc')
        
        # Open filtered and subsampled AVISO data
        nc = Dataset('aviso_filt.nc', 'r')
        latitude = nc.variables['lat'][:]
        longitude = nc.variables['lon'][:]
        ssha_sampled = np.squeeze(nc.variables['z'][:])
        nc.close()
        os.system('rm aviso_filt.nc')
    
        # Pad out empty meridional sections at each end
        ssha_sampled[:, 0] = ssha_sampled[:, 1]
        ssha_sampled[:, -1] = ssha_sampled[:, -2]
        
        ### Load the CryoSat-2 altimetry ssha
        cs2_file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year + month + '_grid.nc'
        
        nc = Dataset(cs2_file, 'r')
        longitude = nc.variables['lon'][:]
        latitude = nc.variables['lat'][:]
        cs2_ssh = nc.variables['sea_surface_height_anomaly_seasonal_offset'][:]
        if year == '2015' and month == '09':
            ice = nc.variables['filtered_sea_ice_concentration'][:]
        nc.close()
        
        anomaly_anomaly[:, :, ix] = cs2_ssh - ssha_sampled
        
        ix += 1

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
grid_lats, grid_lons = np.meshgrid(latitude, longitude)
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(np.nanmean(anomaly_anomaly, axis=2))), cmap='RdBu_r')
c = m.colorbar()
pl.clim(-0.1, 0.1)
m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice)), [20,], colors='k')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/AVISO/mean_CS2-AVISO.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()