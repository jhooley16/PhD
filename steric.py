import numpy as np
from netCDF4 import Dataset
from datetime import date
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as pl
import functions as funct
from scipy import signal
import os

steric_ts = []
dot_ts = []
baristatic_ts = []

dates = []

# Load the SAM index file
sam_file = '/Users/jmh2g09/Documents/PhD/Data/Wind/SAM_2010-2016.txt'
sam_index = []
f = open(sam_file, 'r')
for line in f:
    line = line.strip()
    sam_index.append(float(line))
f.close()
sam_index = np.array(sam_index)

steric_seasonal = np.full((12, 5), fill_value=np.NaN)
dot_seasonal = np.full((12, 5), fill_value=np.NaN)
baristatic_seasonal = np.full((12, 5), fill_value=np.NaN)
ice_extent = np.full((12, 5), fill_value=np.NaN)
wind_seasonal = np.full((12, 5), fill_value=np.NaN)
SAM_seasonal  = np.full((12, 5), fill_value=np.NaN)

steric_months = np.full((59, 361, 12, 7), fill_value=np.NaN)
ice_months = np.full((59, 361, 12, 7), fill_value=np.NaN)
trigger = False
A = 0
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    # Cycle through the months
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        DOT_file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year + '/Anomalies/' + year + month + '_DOT_anomaly.nc'
        if os.path.isfile(DOT_file):
            dates.append(date(int(year), int(month), 15))
            
            # Load DOT
            nc = Dataset(DOT_file, 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            ice = nc.variables['sea_ice_concentration'][:]
            nc.close()
            
            ## Filter the dot file to match the resolution of GRACE
            os.system('gmt grdfilter ' + DOT_file + '?"dynamic_ocean_topography_anomaly_seasonal_offset" -D4 -Fg6000 -Nr -f0y -f1x -GFILT.nc')
            
            nc = Dataset('FILT.nc', 'r')
            dot = nc.variables['z'][:]
            nc.close()
            
            os.system('rm FILT.nc')
            
#             pl.figure()
#             pl.clf()
#             m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
#             m.drawmapboundary()
#             m.drawcoastlines(zorder=10)
#             m.fillcontinents(zorder=10)
#             m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
#             m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
# 
#             grid_lats, grid_lons = np.meshgrid(lat, lon)
#             stereo_x, stereo_y = m(grid_lons, grid_lats)
# 
#             m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot)), cmap='RdBu_r')
#             m.colorbar()
#             pl.clim([-0.1, 0.1])
#             pl.show()
#             pl.close()
#             pause
            
            dot[np.isnan(dot)] = 999
            
            ice_months[:, :, int(month)-1, int(year)-2010] = ice
            
            grid_lats, grid_lons = np.meshgrid(lat, lon)
            
            GRACE_file = '/Users/jmh2g09/Documents/PhD/Data/GRACE/' + year + month + '_GRACE.nc'
            # Load GRACE
            nc = Dataset(GRACE_file, 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            baristatic = nc.variables['GRACE'][:] # (lat, lon)
            nc.close()
            
            wind_file = '/Users/jmh2g09/Documents/PhD/Data/Wind/' + year + month + '_wind.nc'
            # Load GRACE
            nc = Dataset(wind_file, 'r')
            wind_lat = nc.variables['latitude'][:]
            wind = nc.variables['u_wind'][:] # (lat, lon)
            nc.close()
            
            # Apply the GMT land mask
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/mask.nc', 'r')
            # Load the mask (ocean == 1)
            ocean_mask = nc.variables['z'][:]
            ocean_mask[ocean_mask != 1] = np.NaN
            ice_mask = nc.variables['z'][:]
            ice_mask[ice_mask != 1] = np.NaN
            nc.close()
            
            ocean_mask[np.isnan(baristatic)] = np.NaN
            ocean_mask[np.isnan(dot)] = np.NaN
            
            ice_mask[np.isnan(baristatic)] = np.NaN
            ice_mask[dot > 100] = np.NaN
            ice_mask[ice < 20] = np.NaN
            
            dot[dot > 100] = np.NaN
            
            # Calculate the steric
            steric = dot - baristatic
                        
            steric_months[:, :, int(month)-1, int(year)-2010] = steric
            
            ## Calculate the surface area of each cell
            # Mesh the lat and lon together calculate the surface area for each cell
            grid_lon, grid_lat = np.meshgrid(lon, lat)
            # Calculate the surface area of each cell
            S = funct.surface_area(grid_lat, grid_lon, 0.5, 1.0)
            # Get total area for the ocean region
            total_area = np.nansum(ocean_mask * S)

            total_ice_area = np.nansum(ice_mask * S)

            # Calculate weighed mean for all the ocean
            steric_ts.append(np.nansum(steric * S) / total_area)
            dot_ts.append(np.nansum(dot * S) / total_area)
            baristatic_ts.append(np.nansum(baristatic * S) / total_area)

            if int(year) >= 2011 and int(year) < 2016:
                it = int(year) - 2011
                # Calculate the monthly contribution for all the ocean
                steric_seasonal[int(month) - 1, it] = np.nansum(steric * S) / total_area
                dot_seasonal[int(month) - 1, it] = np.nansum(dot * S) / total_area
                baristatic_seasonal[int(month) - 1, it] = np.nansum(baristatic * S) / total_area
                ice_extent[int(month) - 1, it] = total_ice_area
                wind_seasonal[int(month) - 1, it] = np.nanmean(wind[0, :])
                SAM_seasonal[int(month) - 1, it] = sam_index[A]
            
            A += 1

fig = pl.figure()
pl.plot(dates, signal.detrend(np.array(dot_ts)), label='DOT')
pl.plot(dates, signal.detrend(np.array(baristatic_ts)) + .1, label='GRACE')
pl.plot(dates, signal.detrend(np.array(steric_ts)) - .1, label='Steric DOT')
#pl.plot(dates, signal.detrend(np.array(sam_index)) - .2, label='SAM index')
pl.legend(loc='best', prop={'size':8})
pl.xticks(rotation='vertical')
fig.autofmt_xdate()
pl.ylabel('DOT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Steric/steric_timeseries.png', transparent=True, doi=300, bbox_inches='tight')
pl.close()

steric_season = np.nanmean(steric_seasonal, axis=1)
dot_season = np.nanmean(dot_seasonal, axis=1)
baristatic_season = np.nanmean(baristatic_seasonal, axis=1)
ice_extent_season = np.nanmean(ice_extent, axis=1)
wind_season = np.nanmean(wind_seasonal, axis=1)
SAM_season = np.nanmean(SAM_seasonal, axis=1)

fig = pl.figure()
pl.plot(range(1, 13), np.array(dot_season) * 100, label='DOT')
pl.plot(range(1, 13), np.array(baristatic_season) * 100, label='GRACE')
pl.plot(range(1, 13), np.array(steric_season) * 100, label='Steric DOT')
pl.plot(range(1, 13), np.array(ice_extent_season) / 5000000, label='Ice Extent (*5e6)')
#pl.plot(range(1, 13), wind_season - 7., label='wind 10m u-component (+ 7)')
#pl.plot(range(1, 13), SAM_season, label='SAM Index')
pl.xticks(range(1, 13))
pl.xlim([1, 12])
pl.legend(loc='best', prop={'size':8})
pl.ylabel('DOT (cm)')
pl.xlabel('Month')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Steric/steric_seasonal.png', transparent=True, doi=300, bbox_inches='tight')
pl.close()

steric_months_avg = np.nanmean(steric_months, axis=3)
steric_months_avg[steric_months_avg > 100] = np.NaN

ice_months_avg = np.nanmean(ice_months, axis=3)

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
    
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(steric_months_avg[:, :, imnth])), cmap='RdBu_r')
    m.colorbar()
    pl.clim([-0.1, 0.1])
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_months_avg[:, :, imnth])), colors='k', levels=[20])

    pl.title('Steric height estimated from (Altimetry - GRACE), (m)')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Steric/steric_' + str(imnth + 1) + '.png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()