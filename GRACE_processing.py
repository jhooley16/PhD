from netCDF4 import Dataset
import numpy as np
from scipy import interpolate, signal
from datetime import date
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as pl
import os
import functions as funct

### Sea Level Pressure
SSP_file = '/Users/jmh2g09/Documents/PhD/Data/SSPressure/ERA_Interim_mean_sea_level_pressure.nc'

# Load SSP data
nc = Dataset(SSP_file, 'r')
SSP_lat = nc.variables['latitude'][:]
SSP_lon = nc.variables['longitude'][:]
SSP_time = nc.variables['time'][:] / 24 + date.toordinal(date(1900, 1, 1))# Convert to days
SSP_dates = []
for it in range(len(SSP_time)):
    SSP_dates.append(date.fromordinal(int(SSP_time[it])))
# SSP (Pa)
SSP_data = nc.variables['msl'][:]# (time, lat, lon)
nc.close()

# Apply the GMT land mask
nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/mask_full_05degree.nc', 'r')
# Load the mask (ocean == 1)
ocean_mask_05 = np.flipud(nc.variables['z'][:])
ocean_mask_05[ocean_mask_05 != 1] = np.NaN
nc.close()

SSP_data_anom_masked = np.full(np.shape(SSP_data), fill_value=np.NaN)
for it in range(len(SSP_time)):
    SSP_data_anom_masked[it, :, :] = SSP_data[it, :, :] * ocean_mask_05

# Calculate the surface area of each cell
grid_lon, grid_lat = np.meshgrid(SSP_lon, SSP_lat)
S = funct.surface_area(grid_lat, grid_lon, 0.5, 0.5)
# Get total area for the ocean region
total_area = np.nansum(ocean_mask_05 * S)

# Calculate the area-weighted mean of the data
SSPA_ts = np.full(len(SSP_time), fill_value=np.NaN)
for it in range(len(SSP_time)):
    SSPA_ts[it] = np.nansum(SSP_data_anom_masked[it, :, :] * S) / total_area

SSPA_ts_mean = np.nanmean(SSPA_ts)

SSPA_ts_anom = (SSPA_ts - SSPA_ts_mean) / (1025*9.81)

fig = pl.figure()
pl.plot(SSP_dates, SSPA_ts_anom)
fig.autofmt_xdate()
pl.ylabel('OBP anomaly from Atmosphere (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/SSPressure/SSP_anomaly_global_ocean.png')
pl.close()

#################### GRACE Data ############################

GRACE_file = '/Users/jmh2g09/Documents/PhD/Data/GRACE/GRCTellus.JPL.200208_201603.OCN.RL05_1.DSTvDPC1412.nc'

## Load the GRACE data from the raw file. 
nc = Dataset(GRACE_file, 'r')
GRACE_lat = nc.variables['lat'][:]
GRACE_lon = nc.variables['lon'][:]
GRACE_time = nc.variables['time_bounds'][:]
# Convert GRACE lwe thickness to m (from cm)
GRACE_data_pre_2 = nc.variables['lwe_thickness'][:] / 100# (time, lat, lon)
nc.close()

GRACE_lon = np.append([-.5], GRACE_lon)
GRACE_lon = np.append(GRACE_lon, [360.5])

GRACE_data_pre = np.full((149, 180, 362), fill_value=np.NaN)
for ilon in range(362):
    if ilon == 0:
        GRACE_data_pre[:, :, ilon] = GRACE_data_pre_2[:, :, -1]
    elif 0 < ilon <= 360:
        GRACE_data_pre[:, :, ilon] = GRACE_data_pre_2[:, :, ilon - 1]
    elif ilon == 361:
        GRACE_data_pre[:, :, ilon] = GRACE_data_pre_2[:, :, 0]

# Apply the GMT land mask
nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/Masks/mask.nc', 'r')
# Load the mask (ocean == 1)
ocean_mask = nc.variables['z'][:]
ocean_mask[ocean_mask != 1] = np.NaN
nc.close()

GRACE_data_pre[GRACE_data_pre<-900] = np.NaN

## Resample the GRACE grid
GRACE_data = np.full((len(GRACE_time), 59, 361), fill_value=np.NaN)
for it in range(len(GRACE_time)):
    nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/GRACE/INPUT.nc', 'w')
    nc.createDimension('lat', np.size(GRACE_lat))
    nc.createDimension('lon', np.size(GRACE_lon))
    latitude = nc.createVariable('lat', float, ('lat',))
    longitude = nc.createVariable('lon', float, ('lon',))
    GRACE_save = nc.createVariable('GRACE', float, ('lat','lon',))
    latitude[:] = GRACE_lat
    longitude[:] = GRACE_lon
    GRACE_save[:] = GRACE_data_pre[it, :, :]
    nc.close()
    
    os.system('gmt grdsample /Users/jmh2g09/Documents/PhD/Data/GRACE/INPUT.nc \
        -G/Users/jmh2g09/Documents/PhD/Data/GRACE/OUTPUT.nc -I1.0/0.5 -R0/360/-79/-50')
    os.system('rm /Users/jmh2g09/Documents/PhD/Data/GRACE/INPUT.nc')
    
    nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/GRACE/OUTPUT.nc', 'r')
    lat = nc.variables['y'][:]
    lon = nc.variables['x'][:]
    GRACE_data[it, :, :] = nc.variables['z'][:] * ocean_mask
    nc.close()
    os.system('rm /Users/jmh2g09/Documents/PhD/Data/GRACE/OUTPUT.nc')

# Load the GRACE data within the time bounds that I am working with
GRACE_year = []
GRACE_month = []
GRACE_day = []
GRACE_data_post = []
test = False
t_correction = date.toordinal(date(2002, 1, 1))
for i in range(len(GRACE_time)):
    GRACE_year_start = date.timetuple(date.fromordinal(int(GRACE_time[i, 0]) + t_correction))[0]
    GRACE_month_start = date.timetuple(date.fromordinal(int(GRACE_time[i, 0]) + t_correction))[1]
    GRACE_day_start = date.timetuple(date.fromordinal(int(GRACE_time[i, 0]) + t_correction))[2]
    GRACE_year_end = date.timetuple(date.fromordinal(int(GRACE_time[i, 1]) + t_correction))[0]
    GRACE_month_end = date.timetuple(date.fromordinal(int(GRACE_time[i, 1]) + t_correction))[1]
    GRACE_day_end = date.timetuple(date.fromordinal(int(GRACE_time[i, 1]) + t_correction))[2]
    
    if GRACE_year_start == 2010 and GRACE_month_start == 11:
        test = True
    if GRACE_year_start == 2016 and GRACE_month_start == 3:
        test = False
    
    if test == True:
        GRACE_year.append([GRACE_year_start, GRACE_year_end])
        GRACE_month.append([GRACE_month_start, GRACE_month_end])
        GRACE_day.append([GRACE_day_start, GRACE_day_end])
        GRACE_data_post.append(GRACE_data[i, :, :])

GRACE_data_post = np.array(GRACE_data_post)

GRACE_data_post[np.isnan(GRACE_data_post)] = 999

GRACE_data_post_2 = signal.detrend(GRACE_data_post, axis=0)
GRACE_data_post_2[GRACE_data_post > 100] = np.NaN

# Load the missing month data
missing_years = []
missing_months = []
missing_file = open('/Users/jmh2g09/Documents/PhD/Data/GRACE/GRACE_missing_months.csv', 'r')
for line in missing_file:
    row = line.strip()
    columns = row.split(',')
    missing_years.append(columns[0])
    missing_months.append(columns[1])

# Find the indices where the missing months are located
grace_ts = []
month_ts = []
trigger = False
month_step = 0
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    # Cycle through the months
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if year == '2010' and month == '11':
            trigger = True
        if year == '2016' and month == '03':
            trigger = False
        
        if trigger == True:
            MISSING = False
            for imiss in range(len(missing_years)):
                if int(year) == int(missing_years[imiss]) and int(month) == int(missing_months[imiss]):
                    MISSING = True

            if MISSING == False:
                grace_ts.append(month_step)
            
            month_ts.append(month_step)
            month_step += 1

# Interpolate the missing months
GRACE = np.full((64, 59, 361), fill_value=np.NaN)
for ilat in range(59):
    for ilon in range(361):
        interp_ts = interpolate.interp1d(grace_ts, GRACE_data_post_2[:, ilat, ilon], kind='linear')
        # Interpolate
        GRACE[:, ilat, ilon] = interp_ts(month_ts)

for it in range(64):
    GRACE[it, :, :] = GRACE[it, :, :] * ocean_mask

## Save the interpolated data in seperate .nc files

GRACE_months = np.full((12, len(lat), len(lon), 7), fill_value=np.NaN)

SSP_time_mean = []
trigger = False
A = 0
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    # Cycle through the months
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if year == '2010' and month == '11':
            trigger = True
        if year == '2016' and month == '03':
            trigger = False
        
        if trigger == True:
            # Subtract the global atmospheric pressure anomaly
            grace_corrected = GRACE[A, :, :] - SSPA_ts_anom[A]
            
            GRACE_months[int(month)-1, :, :, int(year)-2010] = grace_corrected

            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/GRACE/' + year + month + '_GRACE.nc', 'w')
            nc.createDimension('lat', np.size(lat))
            nc.createDimension('lon', np.size(lon))
            latitude = nc.createVariable('latitude', float, ('lat',))
            longitude = nc.createVariable('longitude', float, ('lon',))
            GRACE_save = nc.createVariable('GRACE', float, ('lat','lon',))

            latitude.long_name = 'latitude'
            latitude.standard_name = 'latitude'
            latitude.units = 'degrees_north'
            longitude.long_name = 'longitude'
            longitude.standard_name = 'longitude'
            longitude.units = 'degrees_east'
            GRACE_save.long_name = 'GRACE measurements'
            GRACE_save.standard_name = 'ocean_gravity_changes_equivalent_water_thickness'
            GRACE_save.units = 'm'

            latitude[:] = lat
            longitude[:] = lon
            GRACE_save[:] = grace_corrected

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

            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(grace_corrected)), cmap='RdBu_r')
            m.colorbar()
            pl.clim([-0.1, 0.1])
            pl.title('GRACE gravity anomalies Units: seawater thickness (m)')
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/GRACE/Figures/' + year + month + '_GRACE.png', transparent=True, dpi=300, bbox_inches='tight')
            pl.close()
            
            A += 1

GRACE_seasonal = np.nanmean(GRACE_months, axis=3)


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

    m.contourf(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(GRACE_seasonal[imnth, :, :])), 15)
    m.colorbar()
    pl.clim([-0.08, 0.09])
    pl.title('GRACE gravity anomalies Units: seawater thickness (m)')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/GRACE/Figures/' + str(imnth+1) + '_GRACE_months.png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()