from netCDF4 import Dataset
import os
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
import numpy as np

# load the monthly along track data

for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if os.path.isfile('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + '/' + year + month + '_raw.nc'):
            print(year, month)

            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + month + '.nc', 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            sea_surface_height = nc.variables['sea_surface_height'][:]
            surface = nc.variables['surface'][:]
            time = nc.variables['time'][:]
            mode = nc.variables['mode'][:]
            ice_conc = nc.variables['ice_concentration'][:]
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
                geoid.append(float(columns[2]))
            f.close()
            os.system('rm /Users/jmh2g09/Documents/PhD/Data/Geoid/OUTPUT.txt')

            # Apply the offset
            sea_surface_height_1 = []
            sea_surface_height_2 = []
            for iel in range(len(lat)):
                # If no offset is required (LRM mode)
                if mode[iel] == 0:
                    sea_surface_height_1.append(sea_surface_height[iel])
                    sea_surface_height_2.append(sea_surface_height[iel])
                if mode[iel] != 0:
                    # If the point is SAR ocean (ocean, SAR)
                    if surface[iel] == 1:
                        sea_surface_height_1.append(sea_surface_height[iel] + funct.apply_offset(month, 'SAR_ocean'))
                        sea_surface_height_2.append(sea_surface_height[iel] + funct.apply_offset('constant', 'SAR_ocean'))
                    # If the point is SAR ice (ice, SAR)
                    if surface[iel] == 2:
                        sea_surface_height_1.append(sea_surface_height[iel] + funct.apply_offset(month, 'ice'))
                        sea_surface_height_2.append(sea_surface_height[iel] + funct.apply_offset('constant', 'ice'))
            sea_surface_height_1 = np.array(sea_surface_height_1)
            sea_surface_height_2 = np.array(sea_surface_height_2)
            
            
            # Calculate DOT
            dot_1 = sea_surface_height_1 - geoid
            dot_2 = sea_surface_height_2 - geoid
            
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
            stereo_x, stereo_y = m(lon[surface == 2], lat[surface == 2])
            m.scatter(stereo_x, stereo_y, c=dot_1[surface == 2], marker='.', edgecolors='none', s=2)
            c = m.colorbar()
            c.set_label('DOT (m)')
            pl.clim(0, -3)
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + '/DOT_track/Figures/' + year + month + '_DOT_ice.png', 
                format='png', transparent=True, dpi=300, bbox_inches='tight')
            pl.close()
                      
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + month + '_ice.nc.nc', 'w')
            nc.createDimension('station', np.size(lat[surface == 2]))

            latitudes = nc.createVariable('latitude', float, ('station',))
            longitudes = nc.createVariable('longitude', float, ('station',))
            sea_surface_height = nc.createVariable('sea_surface_height', float, ('station',))
            surface_type = nc.createVariable('surface', 'i', ('station',))
            sea_ice_concentration = nc.createVariable('ice_concentration', float, ('station',))
            time_var = nc.createVariable('time', float, ('station',))
            mode_var = nc.createVariable('mode', 'i', ('station',))
            geoid_height = nc.createVariable('geoid_height', float, ('station',))
            dynamic_ocean_topography = nc.createVariable('dynamic_ocean_topography_seasonal_offset', float, ('station',))
            dynamic_ocean_topography_2 = nc.createVariable('dynamic_ocean_topography_constant_offset', float, ('station',))
            sea_surface_heights = nc.createVariable('sea_surface_height_seasonal_offset', float, ('station',))
            sea_surface_heights_2 = nc.createVariable('sea_surface_height_constant_offset', float, ('station',))
            
            latitudes.standard_name = 'latitude'
            latitudes.units = 'degrees_north'
            longitudes.standard_name = 'longitude'
            longitudes.units = 'degrees_east'
            surface_type.standard_name = 'surface_type'
            surface_type.units = 'ocean_(1)_or_lead_(2)_surface'
            sea_ice_concentration.standard_name = 'sea_ice_concentration'
            sea_ice_concentration.units = '%'
            sea_surface_height.standard_name = 'sea_surface_height_above_WGS84'
            sea_surface_height.units = 'm'
            time_var.standard_name = 'time'
            time_var.units = 'days_since_1/1/1950'
            mode_var.standard_name = 'satellite_mode_type'
            mode_var.units = 'LRM_(0)_or_SAR_(1)_or_SARIn_(2)'
            geoid_height.standard_name = 'EIGEN6c4_geoid_height_above_WGS84'
            geoid_height.units = 'm'
            dynamic_ocean_topography.standard_name = 'sea_surface_height_above_EIGEN6c4_seasonal_offset'
            dynamic_ocean_topography.units = 'm'
            dynamic_ocean_topography_2.standard_name = 'sea_surface_height_above_EIGEN6c4_constant_offset'
            dynamic_ocean_topography_2.units = 'm'
            sea_surface_heights.standard_name = 'sea_surface_height_above_WGS84_seasonal_offset'
            sea_surface_heights.units = 'm'
            sea_surface_heights_2.standard_name = 'sea_surface_height_above_WGS84_constant_offset'
            sea_surface_heights_2.units = 'm'
            
            latitudes[:] = lat[surface == 2]
            longitudes[:] = lon[surface == 2]
            sea_surface_height[:] = sea_surface_height[surface == 2]
            sea_ice_concentration[:] = ice_conc[surface == 2]
            surface_type[:] = surface[surface == 2]
            time_var[:] = time[surface == 2]
            mode_var[:] = mode[surface == 2]
            geoid_height[:] = geoid[surface == 2]
            dynamic_ocean_topography[:] = dot_1[surface == 2]
            dynamic_ocean_topography_2[:] = dot_2[surface == 2]
            sea_surface_heights[:] = sea_surface_height_1[surface == 2]
            sea_surface_heights_2[:] = sea_surface_height_2[surface == 2]
            
            nc.close()