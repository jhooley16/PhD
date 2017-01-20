import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year)
    # Run through the months and calculate the DOT_track for each month
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        file_exists = False
        for file in os.listdir():
            if file == year + month + '_geoid_track.nc':
                file_exists = True
                nc_geoid = Dataset(file, 'r')
                lat = nc_geoid.variables['latitude'][:]
                lon = nc_geoid.variables['longitude'][:]
                geoid_height = nc_geoid.variables['geoid_height'][:]
                nc_geoid.close()
                
            if file == year + month + '_raw.nc':
                file_exists = True
                nc_ssh = Dataset(file, 'r')
                sea_surface_height = nc_ssh.variables['sea_surface_height_seasonal_offset'][:]
                sea_surface_height_2 = nc_ssh.variables['sea_surface_height_no_offset'][:]
                sea_surface_height_3 = nc_ssh.variables['sea_surface_height_constant_offset'][:]
                
                ice_conc = nc_ssh.variables['ice_concentration'][:]
                surface = nc_ssh.variables['surface_type'][:]
                nc_ssh.close()
                
        if file_exists == True:
            print(year, month)
            dot = sea_surface_height[surface == 2] - geoid_height[surface == 2]
            dot_2 = sea_surface_height_2[surface == 2] - geoid_height[surface == 2]
            dot_3 = sea_surface_height_3[surface == 2] - geoid_height[surface == 2]
            
            lat_2 = lat[surface == 2]
            lon_2 = lon[surface == 2]
            ice_2 = ice_conc[surface == 2]
            surface_2 = surface[surface == 2]
            
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
            stereo_x, stereo_y = m(lon_2, lat_2)
            m.scatter(stereo_x, stereo_y, c=dot, marker='.', edgecolors='none')
            m.colorbar()
            pl.clim(np.mean(dot) + 3*np.std(dot), np.mean(dot) - 3*np.std(dot))
            pl.savefig('DOT_track/Figures/' + year + month + '_DOT_ice_only.png', 
                format='png', dpi=300)
            pl.close()

            nc_dot = Dataset('DOT_track/' + year + month + '_DOT_track_ice_only.nc', 'w')
            nc_dot.createDimension('station', np.size(lat_2))
            latitudes = nc_dot.createVariable('latitude', float, ('station',))
            longitudes = nc_dot.createVariable('longitude', float, ('station',))
            dynamic_ocean_topography = nc_dot.createVariable('dynamic_ocean_topography_seasonal_offset', float, ('station',))
            dynamic_ocean_topography_2 = nc_dot.createVariable('dynamic_ocean_topography_no_offset', float, ('station',))
            dynamic_ocean_topography_3 = nc_dot.createVariable('dynamic_ocean_topography_constant_offset', float, ('station',))
            sea_ice_concentration = nc_dot.createVariable('sea_ice_concentration', float, ('station',))
            surface_type = nc_dot.createVariable('surface_type', float, ('station',))
            
            latitudes.long_name = 'latitude'
            latitudes.standard_name = 'latitude'
            latitudes.units = 'degrees_north'
            longitudes.long_name = 'longitude'
            longitudes.standard_name = 'longitude'
            longitudes.units = 'degrees_east'
            dynamic_ocean_topography.long_name = 'dynamic_ocean_topography_seasonal_offset'
            dynamic_ocean_topography.standard_name = 'sea_surface_height_above_EGM08_geoid'
            dynamic_ocean_topography.units = 'm'
            dynamic_ocean_topography_2.long_name = 'dynamic_ocean_topography_no_offset'
            dynamic_ocean_topography_2.standard_name = 'sea_surface_height_above_EGM08_geoid_no_retracker_offset'
            dynamic_ocean_topography_2.units = 'm'
            dynamic_ocean_topography_3.long_name = 'dynamic_ocean_topography_constant_offset'
            dynamic_ocean_topography_3.standard_name = 'sea_surface_height_above_EGM08_geoid_no_retracker_offset'
            dynamic_ocean_topography_3.units = 'm'
            surface_type.standard_name = 'surface_type'
            surface_type.long_name = 'ocean_(1)_or_lead_(2)_surface'
            sea_ice_concentration.long_name = 'sea_ice_concentration'
            sea_ice_concentration.standard_name = 'sea_ice_concentration'
            sea_ice_concentration.units = '%'
    
            latitudes[:] = lat_2
            longitudes[:] = lon_2
            dynamic_ocean_topography[:] = dot
            dynamic_ocean_topography_2[:] = dot_2
            dynamic_ocean_topography_3[:] = dot_3
            surface_type[:] = surface_2
            sea_ice_concentration[:] = ice_2
    
            nc_dot.close()