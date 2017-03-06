import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
import functions as funct

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
                sea_surface_height = nc_ssh.variables['sea_surface_height'][:]
                surface = nc_ssh.variables['surface'][:]
                mode = nc_ssh.variables['mode'][:]
                ice_conc = nc_ssh.variables['ice_concentration'][:]
                nc_ssh.close()

        if file_exists == True:
            print(year, month)
            sea_surface_height_1 = []
            sea_surface_height_2 = []
            sea_surface_height_3 = []
            # Apply the offset
            for iel in range(len(lat)):
                # If no offset is required (LRM mode)
                if mode[iel] == 0:
                    sea_surface_height_1.append(sea_surface_height[iel])
                    sea_surface_height_2.append(sea_surface_height[iel])
                    sea_surface_height_3.append(sea_surface_height[iel])
                
                if mode[iel] != 0:
                    # If the point is SAR ocean (ocean, SAR)
                    if surface[iel] == 1:
                        sea_surface_height_1.append(sea_surface_height[iel] + funct.apply_offset(month, 'SAR_ocean'))
                        sea_surface_height_2.append(sea_surface_height[iel] + funct.apply_offset('constant', 'SAR_ocean'))
                        sea_surface_height_3.append(sea_surface_height[iel])
                    # If the point is SAR ice (ice, SAR)
                    if surface[iel] == 2:
                        sea_surface_height_1.append(sea_surface_height[iel] + funct.apply_offset(month, 'ice'))
                        sea_surface_height_2.append(sea_surface_height[iel] + funct.apply_offset('constant', 'ice'))
                        sea_surface_height_3.append(sea_surface_height[iel])
            
            sea_surface_height_1 = np.array(sea_surface_height_1)
            sea_surface_height_2 = np.array(sea_surface_height_2)
            sea_surface_height_3 = np.array(sea_surface_height_3)
            
            # Calculate DOT
            dot_1 = sea_surface_height_1 - geoid_height
            dot_2 = sea_surface_height_2 - geoid_height
            dot_3 = sea_surface_height_3 - geoid_height
            
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
            stereo_x, stereo_y = m(lon, lat)
            m.scatter(stereo_x, stereo_y, c=dot_1, marker='.', edgecolors='none', s=2)
            c = m.colorbar()
            c.set_label('DOT (m)')
            pl.clim(0, -3)
            pl.savefig('DOT_track/Figures/' + year + month + '_DOT.png', 
                format='png', transparent=True, dpi=300, bbox_inches='tight')
            pl.close()

            nc_dot = Dataset('DOT_track/' + year + month + '_DOT_track.nc', 'w')
            nc_dot.createDimension('station', np.size(lat))
            latitudes = nc_dot.createVariable('latitude', float, ('station',))
            longitudes = nc_dot.createVariable('longitude', float, ('station',))
            dynamic_ocean_topography = nc_dot.createVariable('dynamic_ocean_topography_seasonal_offset', float, ('station',))
            dynamic_ocean_topography_2 = nc_dot.createVariable('dynamic_ocean_topography_constant_offset', float, ('station',))
            dynamic_ocean_topography_3 = nc_dot.createVariable('dynamic_ocean_topography_no_offset', float, ('station',))
            sea_ice_concentration = nc_dot.createVariable('sea_ice_concentration', float, ('station',))

            latitudes.long_name = 'latitude'
            latitudes.standard_name = 'latitude'
            latitudes.units = 'degrees_north'
            longitudes.long_name = 'longitude'
            longitudes.standard_name = 'longitude'
            longitudes.units = 'degrees_east'
            dynamic_ocean_topography.long_name = 'dynamic_ocean_topography_seasonal_offset'
            dynamic_ocean_topography.standard_name = 'sea_surface_height_above_EIGEN6c4_geoid_seasonal_retracker_offset'
            dynamic_ocean_topography.units = 'm'
            dynamic_ocean_topography_2.long_name = 'dynamic_ocean_topography_constant_offset'
            dynamic_ocean_topography_2.standard_name = 'sea_surface_height_above_EIGEN6c4_geoid_constant_retracker_offset'
            dynamic_ocean_topography_2.units = 'm'
            dynamic_ocean_topography_3.long_name = 'dynamic_ocean_topography_no_offset'
            dynamic_ocean_topography_3.standard_name = 'sea_surface_height_above_EIGEN6c4_geoid_no_retracker_offset'
            dynamic_ocean_topography_3.units = 'm'
            sea_ice_concentration.long_name = 'sea_ice_concentration'
            sea_ice_concentration.standard_name = 'sea_ice_concentration'
            sea_ice_concentration.units = '%'
    
            latitudes[:] = lat
            longitudes[:] = lon
            dynamic_ocean_topography[:] = dot_1
            dynamic_ocean_topography_2[:] = dot_2
            dynamic_ocean_topography_3[:] = dot_3
            sea_ice_concentration[:] = ice_conc
    
            nc_dot.close()