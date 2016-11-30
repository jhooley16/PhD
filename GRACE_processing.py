from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
from datetime import date
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as pl
import os
import functions as funct

GRACE_file = '/Users/jmh2g09/Documents/PhD/Data/GRACE/GRCTellus.JPL.200208_201603.OCN.RL05_1.DSTvDPC1412.nc'

## Load the GRACE data from the raw file. 
nc = Dataset(GRACE_file, 'r')
GRACE_lat = nc.variables['lat'][:]
GRACE_lon = nc.variables['lon'][:]
GRACE_time = nc.variables['time_bounds'][:]
# Convert GRACE lwe thickness to m (from cm) 
GRACE_data_pre = nc.variables['lwe_thickness'][:] / 100# (time, lat, lon)
nc.close()

GRACE_lon = np.append(GRACE_lon, 360.5)
GRACE_data_pre = np.dstack((GRACE_data_pre, GRACE_data_pre[:, :, 0]))

SSP_file = '/Users/jmh2g09/Documents/PhD/Data/SSPressure/ERA_Interim_sea_surface_pressure.nc'

## Load SSP data
nc = Dataset(SSP_file, 'r')
SSP_lat = nc.variables['latitude'][:]
SSP_lon = nc.variables['longitude'][:]
SSP_time = nc.variables['time'][:] / 24 # Convert to days
# Convert SSP (Pa) to m of water equivalent
SSP_data_pre = nc.variables['msl'][:] / (-1025*9.81)# (time, lat, lon)
nc.close()

SSP_lon = np.append(SSP_lon, 360)
SSP_data_pre = np.dstack((SSP_data_pre, SSP_data_pre[:, :, 0]))

# Mask out missing values (land, etc) 
GRACE_data_pre[GRACE_data_pre<-999] = np.NaN

## Load the GRACE data within the time bounds that I am working with
GRACE_year = []
GRACE_month = []
GRACE_day = []
GRACE_data = []
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
        GRACE_data.append(GRACE_data_pre[i, :, :])

GRACE_data = np.array(GRACE_data)

## Load the SSP data within the time bounds that I am working with
SSP_year = []
SSP_month = []
SSP_day = []
SSP_data = []
SSP_dates = []
test = False
t_correction = date.toordinal(date(1900, 1, 1))
for i in range(len(SSP_time)):
    iSSP_year = date.timetuple(date.fromordinal(int(SSP_time[i]) + t_correction))[0]
    iSSP_month = date.timetuple(date.fromordinal(int(SSP_time[i]) + t_correction))[1]
    iSSP_day = date.timetuple(date.fromordinal(int(SSP_time[i]) + t_correction))[2]
    
    if iSSP_year == 2010 and iSSP_month == 11:
        test = True
    if iSSP_year == 2016 and iSSP_month == 3:
        test = False
    
    if test == True:
        SSP_year.append(iSSP_year)
        SSP_month.append(iSSP_month)
        SSP_day.append(iSSP_day)
        SSP_dates.append(date(int(iSSP_year), int(iSSP_month), 15))
        SSP_data.append(SSP_data_pre[i, :, :])

SSP_data = np.array(SSP_data)

SSP_mean = np.nanmean(SSP_data, axis=0)
for it in range(64):
    SSP_data[it, :, :] = SSP_data[it, :, :] - SSP_mean
    
## Load the missing month data
missing_years = []
missing_months = []
missing_file = open('/Users/jmh2g09/Documents/PhD/Data/GRACE/GRACE_missing_months.csv', 'r')
for line in missing_file:
    row = line.strip()
    columns = row.split(',')
    missing_years.append(columns[0])
    missing_months.append(columns[1])

## Find the indices where the missing months are located
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

            if MISSING == True:
                pass
            else:
                grace_ts.append(month_step)
            month_ts.append(month_step)
            month_step += 1

## Interpolate the missing months
GRACE = np.full((64, 180, 361), fill_value=np.NaN)
for ilat in range(180):
    for ilon in range(361):
        interp_ts = interpolate.interp1d(grace_ts, GRACE_data[:, ilat, ilon], kind='cubic')
        GRACE[:, ilat, ilon] = interp_ts(month_ts)

GRACE_mean = np.nanmean(GRACE, axis=0)
for it in range(64):
    GRACE[it, :, :] = GRACE[it, :, :] - GRACE_mean

## Save the interpolated data in seperate .nc files
## Also resample onto a 0.5 lat, 1.0 lon grid

SSP_time_mean = []
trigger == False
A = 0
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    # Cycle through the months
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if year == '2010' and month == '11':
            trigger = True
        if year == '2016' and month == '03':
            trigger = False
        
        if trigger == True:
            ## GRACE
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/GRACE/INPUT.nc', 'w')
            nc.createDimension('lat', np.size(GRACE_lat))
            nc.createDimension('lon', np.size(GRACE_lon))
            latitude = nc.createVariable('lat', float, ('lat',))
            longitude = nc.createVariable('lon', float, ('lon',))
            GRACE_save = nc.createVariable('GRACE', float, ('lat','lon',))
            latitude[:] = GRACE_lat
            longitude[:] = GRACE_lon
            GRACE_save[:] = GRACE[A, :, :]
            nc.close()
            
            ## ERA Interim Pressure
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/SSPressure/INPUT.nc', 'w')
            nc.createDimension('lat', np.size(SSP_lat))
            nc.createDimension('lon', np.size(SSP_lon))
            latitude = nc.createVariable('lat', float, ('lat',))
            longitude = nc.createVariable('lon', float, ('lon',))
            SSP_save = nc.createVariable('SSP', float, ('lat','lon',))
            latitude[:] = SSP_lat
            longitude[:] = SSP_lon
            SSP_save[:] = SSP_data[A, :, :]
            nc.close()
            
            os.system('gmt grdsample /Users/jmh2g09/Documents/PhD/Data/GRACE/INPUT.nc \
                -G/Users/jmh2g09/Documents/PhD/Data/GRACE/OUTPUT.nc -I1.0/0.5 -R0/360/-79/-50')
            os.system('rm /Users/jmh2g09/Documents/PhD/Data/GRACE/INPUT.nc')
            
            os.system('gmt grdsample /Users/jmh2g09/Documents/PhD/Data/SSPressure/INPUT.nc \
                -G/Users/jmh2g09/Documents/PhD/Data/SSPressure/OUTPUT.nc -I1.0/0.5 -R0/360/-79/-50')
            os.system('rm /Users/jmh2g09/Documents/PhD/Data/SSPressure/INPUT.nc')
            
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/GRACE/OUTPUT.nc', 'r')
            lat = nc.variables['y'][:]
            lon = nc.variables['x'][:]
            grace_z = nc.variables['z'][:]
            nc.close
            ## Because GRACE lon is 0.5:359.5, no data is binned on 0,
            ## Make 0 equal to 360
            grace_z[:, 0] = grace_z[:, -1]
            
            os.system('rm /Users/jmh2g09/Documents/PhD/Data/GRACE/OUTPUT.nc')

            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/SSPressure/OUTPUT.nc', 'r')
            lat = nc.variables['y'][:]
            lon = nc.variables['x'][:]
            ssp_z = nc.variables['z'][:]
            nc.close()
            os.system('rm /Users/jmh2g09/Documents/PhD/Data/SSPressure/OUTPUT.nc')

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
            GRACE_save[:] = grace_z - ssp_z

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
        
            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(grace_z - ssp_z)), cmap='RdBu_r')
            m.colorbar()
            pl.title('GRACE gravity anomalies Units: seawater thickness (m)')
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/GRACE/Figures/' + year + month + '_GRACE.png', transparent=True, dpi=300)
            pl.close()
            
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
        
            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ssp_z)), cmap='RdBu_r')
            m.colorbar()
            pl.title('Sea Surface Pressure Anomaly Unit: seawater thickness (m)')
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/SSPressure/Figures/' + year + month + '_MSLP.png', transparent=True, dpi=300)
            pl.close()
            
            ## Calculate the surface area of each cell
            # Mesh the lat and lon together calculate the surface area for each cell
            grid_lon, grid_lat = np.meshgrid(lon, lat)
            # Calculate the surface area of each cell
            S = funct.surface_area(grid_lat, grid_lon, 0.5, 1.0)
            # Get total area for the ocean region
            total_area = np.nansum(np.nansum(~np.isnan(grace_z) * S))
            
            SSP_time_mean.append(np.nanmean(ssp_z * S) / total_area)
            
            A += 1

fig = pl.figure()
pl.plot(SSP_dates, SSP_time_mean)
pl.title('Sea Surface Pressure timeseries')
pl.ylabel('Sea Surface Pressure (Sea water equivalent) (m)')
pl.xticks(rotation='vertical')
fig.autofmt_xdate()
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/SSPressure/Figures/SSP_timeseries.png', transparent=True, dpi=300)
pl.close()