import numpy as np
from netCDF4 import Dataset
import os
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
import functions as funct
from datetime import date

plotting = False

# Open the dynamic topography data
dot = np.full((59, 361, 64), fill_value=np.NaN)
dot_2 = np.full((59, 361, 64), fill_value=np.NaN)
ice = np.full((59, 361, 64), fill_value=np.NaN)
it = 0
months = []
years = []
dates = []
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year)
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        file = year + month + '_DOT_filt.nc'
        if os.path.isfile(file):
            months.append(str(month))
            years.append(str(year))
            dates.append(date(int(year), int(month), 15))
            
            nc = Dataset(file, 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            dot[:, :, it] = nc.variables['dynamic_ocean_topography_seasonal_offset'][:]
            dot_2[:, :, it] = nc.variables['dynamic_ocean_topography_constant_offset'][:]
            ice[:, :, it] = nc.variables['sea_ice_concentration'][:]
            nc.close()
            it += 1

# Mask out data not in the Weddell Sea
lon_mask_min = np.where(lon == 50.)[0][0]
lon_mask_max = np.where(lon == 295.)[0][0]

dot[:, lon_mask_min:lon_mask_max, :] = np.NaN
dot_2[:, lon_mask_min:lon_mask_max, :] = np.NaN
ice[:, lon_mask_min:lon_mask_max, :] = np.NaN

# Sort the longitude so it goes from -180 to 180
lon[lon > 180] -= 361
dot = dot[:, np.argsort(lon), :]
dot_2 = dot_2[:, np.argsort(lon), :]
ice = ice[:, np.argsort(lon), :]
lon = lon[np.argsort(lon)]

# Take the mean
mean_dot = np.nanmean(dot, axis=2)
mean_dot_2 = np.nanmean(dot_2, axis=2)

# Define where the gyre is
gyre_bounds = -2.1

# plot the DOT (mean and individuals) in the Weddell Sea
if plotting == True:
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-90, urcrnrlon=55, urcrnrlat=-40, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot[:, :])), cmap='Blues_r')
    c = m.colorbar()
    c.set_label('DOT (m)')
    pl.clim(-2.3, -1)
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/mean_gyre.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

################ Seasonality of the weddell gyre ################

# Mesh the lat and lon together calculate the surface area for each cell
grid_lon, grid_lat = np.meshgrid(lon, lat)
# Calculate the surface area of each cell
S = funct.surface_area(grid_lat, grid_lon, 0.5, 1.0)
weddell_area = np.nansum(S * np.isfinite(mean_dot))

gyre_mean = []
gyre_mean_2 = []
gyre_seasons = np.full((12, 7), fill_value=np.NaN)
gyre_seasons_2 = np.full((12, 7), fill_value=np.NaN)
for it in range(64):
    # Take the mean away to produce anomalies
    dot_dummy = dot[:, :, it] - mean_dot
    dot_dummy_2 = dot_2[:, :, it] - mean_dot_2
    # Blank out everything outside of the "Weddell Gyre"
    dot_dummy[mean_dot > gyre_bounds] = np.NaN
    dot_dummy_2[mean_dot_2 > gyre_bounds] = np.NaN
    # Calculate the seasonal variability of the gyre
    # Calculate the surface area of the gyre
    gyre_area = np.nansum(S * np.isfinite(dot_dummy))
    gyre_area_2 = np.nansum(S * np.isfinite(dot_dummy_2))
    # Calculate the seasonal variability of the gyre
    gyre_mean.append(np.nansum(dot_dummy * S) / gyre_area)
    gyre_mean_2.append(np.nansum(dot_dummy_2 * S) / gyre_area_2)
    gyre_seasons[int(months[it])-1, int(years[it])-2010] = np.nansum(dot_dummy * S) / gyre_area
    gyre_seasons_2[int(months[it])-1, int(years[it])-2010] = np.nansum(dot_dummy_2 * S) / gyre_area_2

# Average all the years to produce a seasonal cycle
gyre_seasonal = np.nanmean(gyre_seasons, axis=1)
gyre_seasonal_2 = np.nanmean(gyre_seasons_2, axis=1)

# Plot the timeseries and seasonal cycle of the DOT anomaly in the Weddell Gyre
if plotting == True:
    fig = pl.figure()
    pl.plot(dates, np.array(gyre_mean) * 10, label='Seasonal Offset')
    pl.plot(dates, np.array(gyre_mean_2) * 10, label='Constant Offset')
    fig.autofmt_xdate()
    pl.legend(loc='best', prop={'size':9})
    pl.grid()
    pl.ylabel('DOT Anomaly (cm)')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/gyre_ts.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.plot(range(1, 13), gyre_seasonal * 10)
    #pl.plot(range(1, 13), gyre_seasonal_2 * 10, label='Constant Offset')
    my_xticks = ['Jan','Feb','Mar','Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']
    pl.xticks(range(1, 13), my_xticks)
    pl.xlim(1, 12)
    #pl.legend(loc='best', prop={'size':9})
    pl.grid()
    pl.ylabel('DOT Anomaly (cm)')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/gyre_seasonal.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

################ Weddell Gyre geostrophic velocities ################

# Define the geostrophic equation elements
a = 6378137. # radius of earth (m)
dy = 0.5 * 60 * 1862 # distance between two latitude points (0.5 grid spacing) (m)
dx = (np.pi / 180) * a * np.cos(lat * np.pi / 180)
dx_3d = np.transpose(np.tile(np.tile(dx, (361, 1)), (64, 1, 1)))
g = 9.81 # gravitational acceleration (meters per square second)
omega = 2 * np.pi / (24 * 60 * 60) # rotation rate of earth (radians per second)

# Calculate the coriolis paramater
f = 2 * omega * np.sin(lat * np.pi / 180)
f_3d = np.transpose(np.tile(np.tile(f, (361, 1)), (64, 1, 1)))

# Take the horizontal gradients
ddot_dy = np.gradient(dot, axis=0) / dy
ddot_dx = np.gradient(dot, axis=1) / dx_3d

# Calculate the geostrophic velocities
u =  - (g / f_3d) * ddot_dy
v = (g / f_3d) * ddot_dx

# Take the time mean
u_bar = np.nanmean(u, axis=2)
v_bar = np.nanmean(v, axis=2)

# Calculate absolute velocity magnitude
velocity_bar = (u_bar ** 2 + v_bar ** 2) ** 0.5

# Plot maps of u, v and absolute velocities in Weddell gyre
if plotting == True:
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(u_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-0.05, 0.05])
    c.set_label('u-velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Velocities/velocity_u_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(v_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-0.06, 0.06])
    c.set_label('v-velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Velocities/velocity_v_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(velocity_bar)), cmap='Reds')
    c = m.colorbar()
    pl.clim([0, 0.2])
    c.set_label('velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    u_blank = u_bar
    v_blank = v_bar
    u_blank[mean_dot > gyre_bounds] = np.NaN
    v_blank[mean_dot > gyre_bounds] = np.NaN

    uproj,vproj,xx,yy = m.transform_vector(u_blank,v_blank,lon,lat, 25, 25,returnxy=True,masked=True)
    m.quiver(xx, yy, uproj, vproj, scale=1)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Velocities/velocity_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

################ 'summer' and 'winter' velocities ################

# Calculate the number of months in each 'summer' and 'winter' period
summer_months = ['01', '02', '08', '09', '10', '11', '12']
winter_months = ['03', '04', '05', '06', '07']
summer_count = 0
winter_count = 0
for month in months:
    if month in summer_months:
        summer_count += 1
    elif month in winter_months:
        winter_count += 1

# Split u and v velocities into summer and winter slices
summer_u = np.full((59, 361, summer_count), fill_value=np.NaN)
summer_v = np.full((59, 361, summer_count), fill_value=np.NaN)
winter_u = np.full((59, 361, winter_count), fill_value=np.NaN)
winter_v = np.full((59, 361, winter_count), fill_value=np.NaN)
summer_step = 0
winter_step = 0
it = 0
# Calculate the summer/winter averages
for mnth in months:
    if mnth in summer_months:
        summer_u[:, :, summer_step] = u[:, :, it]
        summer_v[:, :, summer_step] = v[:, :, it]
        summer_step += 1
    if mnth in winter_months:
        winter_u[:, :, winter_step] = u[:, :, it]
        winter_v[:, :, winter_step] = v[:, :, it]
        winter_step += 1
    it += 1

# Take the average for each summer month
summer_u_bar = np.nanmean(summer_u, axis=2)
summer_v_bar = np.nanmean(summer_v, axis=2)

# Calculate a summer absolute average
summer_velocity_bar = (summer_u_bar ** 2 + summer_v_bar ** 2) ** 0.5

# Take the average for each winter month
winter_u_bar = np.nanmean(winter_u, axis=2)
winter_v_bar = np.nanmean(winter_v, axis=2)

# Calculate a winter absolute average
winter_velocity_bar = (winter_u_bar ** 2 + winter_v_bar ** 2) ** 0.5

# Plot winter and summer u, v and absolute velocities
# Also take (winter - summer) absolute velocity anomaly
if plotting == True:
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(summer_u_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-0.05, 0.05])
    c.set_label('u-velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Velocities/velocity_summer_u_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(summer_v_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-0.06, 0.06])
    c.set_label('v-velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Velocities/velocity_summer_v_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(summer_velocity_bar)), cmap='Reds')
    c = m.colorbar()
    pl.clim([0, 0.2])
    c.set_label('velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    u_blank = summer_u_bar
    v_blank = summer_v_bar
    u_blank[mean_dot > gyre_bounds] = np.NaN
    v_blank[mean_dot > gyre_bounds] = np.NaN

    uproj,vproj,xx,yy = m.transform_vector(u_blank,v_blank,lon,lat, 25, 25,returnxy=True,masked=True)
    m.quiver(xx, yy, uproj, vproj, scale=1)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Velocities/velocity_summer_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_u_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-0.05, 0.05])
    c.set_label('u-velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Velocities/velocity_winter_u_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_v_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-0.06, 0.06])
    c.set_label('v-velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Velocities/velocity_winter_v_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_velocity_bar)), cmap='Reds')
    c = m.colorbar()
    pl.clim([0, 0.2])
    c.set_label('Velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    u_blank = winter_u_bar
    v_blank = winter_v_bar
    u_blank[mean_dot > gyre_bounds] = np.NaN
    v_blank[mean_dot > gyre_bounds] = np.NaN

    uproj,vproj,xx,yy = m.transform_vector(u_blank,v_blank,lon,lat, 25, 25,returnxy=True,masked=True)
    m.quiver(xx, yy, uproj, vproj, scale=1)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Velocities/velocity_winter_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_velocity_bar - summer_velocity_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-0.02, 0.02])
    c.set_label('Velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Velocities/velocity_winter-summer_anomaly.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

################ wind velocities and curl ################

u_wind = np.full((59, 361, 64), fill_value=np.NaN)
v_wind = np.full((59, 361, 64), fill_value=np.NaN)
wind_curl = np.full((59, 361, 64), fill_value=np.NaN)
# Calculate dx
dx_2d_curl = np.transpose(np.tile(dx, (361, 1)))
it = 0
# Load wind data and calculate the curl
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Wind/')
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        file = year + month + '_wind.nc'
        if os.path.isfile(file):
            nc = Dataset(file, 'r')
            lon_wind =  nc.variables['longitude'][:]
            lat_wind =  nc.variables['latitude'][:]
            u_wind[:, :, it] = nc.variables['u_wind'][:]
            v_wind[:, :, it] = nc.variables['v_wind'][:]
            # Calculate horizontal gradients for curl
            dv_dx = np.gradient(nc.variables['v_wind'][:], axis=1) / dx_2d_curl
            du_dy = np.gradient(nc.variables['u_wind'][:], axis=0) / dy
            nc.close()
            # Calculate the wind curl
            wind_curl[:, :, it] = dv_dx - du_dy
            it += 1

# Mask out data not in the Weddell Sea
lon_mask_max = np.where(lon_wind == 50.)[0][0]
lon_mask_min = np.where(lon_wind == -65.)[0][0]

u_wind[:, :lon_mask_min, :] = np.NaN
u_wind[:, lon_mask_max:, :] = np.NaN
v_wind[:, :lon_mask_min, :] = np.NaN
v_wind[:, lon_mask_max:, :] = np.NaN

wind_curl[:, :lon_mask_min, :] = np.NaN
wind_curl[:, lon_mask_max:, :] = np.NaN

# Mask out land
u_wind[np.isnan(mean_dot)] = np.NaN
v_wind[np.isnan(mean_dot)] = np.NaN
wind_curl[np.isnan(mean_dot)] = np.NaN

# Take the time mean
mean_u_wind = np.nanmean(u_wind, axis=2)
mean_v_wind = np.nanmean(v_wind, axis=2)

mean_wind_curl = np.nanmean(wind_curl, axis=2)

# Calculate the absolute wind velocity
mean_wind = (mean_u_wind ** 2 + mean_v_wind ** 2) ** 0.5

# Plot the u, v and absolute wind velocity
# Also plot the wind curl
if plotting == True:
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-90, urcrnrlon=55, urcrnrlat=-40, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    #m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat_wind, lon_wind)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_u_wind)), cmap='RdBu_r')
    c = m.colorbar()
    c.set_label('u-wind (m s$^{-1}$)')
    pl.clim(-8, 8)
    #m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_u_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
    
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-90, urcrnrlon=55, urcrnrlat=-40, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    #m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat_wind, lon_wind)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_wind_curl)), cmap='RdBu_r')
    c = m.colorbar()
    c.set_label('Wind Curl (s$^{-1}$)')
    uproj,vproj,xx,yy = m.transform_vector(mean_u_wind,mean_v_wind,lon,lat, 25, 25,returnxy=True,masked=True)
    m.quiver(xx, yy, uproj, vproj, scale=100)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_curl_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
    
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-90, urcrnrlon=55, urcrnrlat=-40, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    #m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat_wind, lon_wind)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_v_wind)), cmap='RdBu_r')
    c = m.colorbar()
    c.set_label('v-wind (m s$^{-1}$)')
    pl.clim(-3, 3)
    #m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_v_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
    
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_wind)))
    c = m.colorbar()
    pl.clim([0, 10])
    c.set_label('wind velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    u_blank = mean_u_wind
    v_blank = mean_v_wind
    u_blank[mean_dot > gyre_bounds] = np.NaN
    v_blank[mean_dot > gyre_bounds] = np.NaN
    uproj,vproj,xx,yy = m.transform_vector(u_blank,v_blank,lon,lat, 25, 25,returnxy=True,masked=True)
    m.quiver(xx, yy, uproj, vproj, scale=100)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

################ seasonality of weddell wind  ################

# Calculate the seasonality of the wind over the Weddell Gyre
gyre_seasons_u_wind = np.full((12, 7), fill_value=np.NaN)
gyre_seasons_v_wind = np.full((12, 7), fill_value=np.NaN)
for it in range(64):
    # Blank out everything outside of the "Weddell Gyre"
    wind_u_dummy = u_wind[:, :, it] - mean_u_wind
    wind_v_dummy = v_wind[:, :, it] - mean_v_wind
    
    wind_u_dummy[mean_dot > gyre_bounds] = np.NaN
    #wind_u_dummy[wind_u_dummy > 0] = np.NaN
    #wind_u_dummy[49:, :] = np.NaN
    wind_v_dummy[mean_dot > gyre_bounds] = np.NaN
    #wind_v_dummy[wind_u_dummy > 0] = np.NaN
    #wind_v_dummy[49:, :] = np.NaN
    
    gyre_mean_u_wind = np.nansum(wind_u_dummy * S) / weddell_area
    gyre_mean_v_wind = np.nansum(wind_v_dummy * S) / weddell_area
    
    gyre_seasons_u_wind[int(months[it])-1, int(years[it])-2010] = gyre_mean_u_wind
    gyre_seasons_v_wind[int(months[it])-1, int(years[it])-2010] = gyre_mean_v_wind

# Average the years to calculate seasonal wind variation
u_wind_seasonal = np.nanmean(gyre_seasons_u_wind, axis=1)
v_wind_seasonal = np.nanmean(gyre_seasons_v_wind, axis=1)

# Plot the seasonality of the u and v wind velocities
# Also plot the seasonality in each individual year
if plotting == True:
    pl.figure()
    pl.plot(range(1, 13), u_wind_seasonal, label='u-wind')
    pl.plot(range(1, 13), v_wind_seasonal, label='v-wind')
    my_xticks = ['Jan','Feb','Mar','Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']
    pl.xticks(range(1, 13), my_xticks)
    pl.xlim(1, 12)
    pl.legend(loc='best', prop={'size':9})
    pl.grid()
    pl.ylabel('10 m wind velocity over gyre (m s$^{-1}$)')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_seasonal.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
    
    pl.figure()
    years_ind = ['2010', '2011', '2012', '2013', '2014', '2015', '2016']
    for it in range(7):
        pl.plot(range(1, 13),gyre_seasons_u_wind[:, it], label=years_ind[it])
    pl.xticks(range(1, 13), my_xticks)
    pl.xlim(1, 12)
    pl.legend(loc='best', prop={'size':9})
    pl.grid()
    pl.ylabel('10 m u-wind velocity over gyre (m s$^{-1}$)')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_years.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

################ 'summer' and 'winter' weddell wind  ################

# Calculate the difference between winter and summer wind states
summer_u_wind = np.full((59, 361, summer_count), fill_value=np.NaN)
summer_v_wind = np.full((59, 361, summer_count), fill_value=np.NaN)
winter_u_wind = np.full((59, 361, winter_count), fill_value=np.NaN)
winter_v_wind = np.full((59, 361, winter_count), fill_value=np.NaN)
summer_curl = np.full((59, 361, summer_count), fill_value=np.NaN)
winter_curl = np.full((59, 361, winter_count), fill_value=np.NaN)
summer_step = 0
winter_step = 0
it = 0
# Calculate the summer/winter averages
for month in months:
    if month in summer_months:
        summer_u_wind[:, :, summer_step] = u_wind[:, :, it]
        summer_v_wind[:, :, summer_step] = v_wind[:, :, it]
        summer_curl[:, :, summer_step] = wind_curl[:, :, it]
        summer_step += 1
    if month in winter_months:
        winter_u_wind[:, :, winter_step] = u_wind[:, :, it]
        winter_v_wind[:, :, winter_step] = v_wind[:, :, it]
        winter_curl[:, :, winter_step] = wind_curl[:, :, it]
        winter_step += 1
    it += 1

# Calculate summer mean
summer_u_wind_bar = np.nanmean(summer_u_wind, axis=2)
summer_v_wind_bar = np.nanmean(summer_v_wind, axis=2)

# Calculate mean summer wind velocity
summer_velocity_wind_bar = (summer_u_wind_bar ** 2 + summer_v_wind_bar ** 2) ** 0.5

# Calculate mean summer wind curl
summer_curl_bar = np.nanmean(summer_curl, axis=2)

# Calculate winter mean
winter_u_wind_bar = np.nanmean(winter_u_wind, axis=2)
winter_v_wind_bar = np.nanmean(winter_v_wind, axis=2)

# Calculate mean winter wind velocity
winter_velocity_wind_bar = (winter_u_wind_bar ** 2 + winter_v_wind_bar ** 2) ** 0.5

# Calculate mean winter wind curl
winter_curl_bar = np.nanmean(winter_curl, axis=2)

# Plot 'summer' and 'winter' wind curl
# Also plot 'summer' and 'winter' u, v and absolute wind velocities
# Also plot summer - winter anomalies for curl, u and v wind
if plotting == True:
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_curl_bar - summer_curl_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-.000013, .000013])
    c.set_label('wind curl (s$^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_curl_winter-summer_anomaly.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
    
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_u_wind_bar - summer_u_wind_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-3, 3])
    c.set_label('wind velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_winter-summer_u_anomaly.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
    
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_v_wind_bar - summer_v_wind_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-1, 1])
    c.set_label('wind velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_winter-summer_v_anomaly.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
    
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(summer_u_wind_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-8, 8])
    c.set_label('u-wind velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_summer_u_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(summer_v_wind_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-3, 3])
    c.set_label('v-wind velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_summer_v_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_u_wind_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-8, 8])
    c.set_label('u-wind velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_winter_u_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_v_wind_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-3, 3])
    c.set_label('v-wind velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_winter_v_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
    
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(summer_velocity_wind_bar)))
    c = m.colorbar()
    pl.clim([0, 10])
    c.set_label('wind velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    u_blank = summer_u_wind_bar
    v_blank = summer_v_wind_bar
    u_blank[mean_dot > gyre_bounds] = np.NaN
    v_blank[mean_dot > gyre_bounds] = np.NaN

    uproj,vproj,xx,yy = m.transform_vector(u_blank,v_blank,lon,lat, 25, 25,returnxy=True,masked=True)
    m.quiver(xx, yy, uproj, vproj, scale=100)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_summer_velocity_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
    
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(summer_curl_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-.00003, .00003])
    c.set_label('wind curl (s$^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')

    uproj,vproj,xx,yy = m.transform_vector(u_blank,v_blank,lon,lat, 25, 25,returnxy=True,masked=True)
    m.quiver(xx, yy, uproj, vproj, scale=100)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_curl_summer_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_velocity_wind_bar)))
    c = m.colorbar()
    pl.clim([0, 10])
    c.set_label('wind velocity (m s $^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')
    u_blank = winter_u_wind_bar
    v_blank = winter_v_wind_bar
    u_blank[mean_dot > gyre_bounds] = np.NaN
    v_blank[mean_dot > gyre_bounds] = np.NaN

    uproj,vproj,xx,yy = m.transform_vector(u_blank,v_blank,lon,lat, 25, 25,returnxy=True,masked=True)
    m.quiver(xx, yy, uproj, vproj, scale=100)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_winter_velocity_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
    
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_curl_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-.00003, .00003])
    c.set_label('wind curl (s$^{-1}$)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='w')

    uproj,vproj,xx,yy = m.transform_vector(u_blank,v_blank,lon,lat, 25, 25,returnxy=True,masked=True)
    m.quiver(xx, yy, uproj, vproj, scale=100)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/wind_curl_winter_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

################ 'summer' and 'winter' Ekman Transport  ################

rho_air = 1.3 # kg/m^3
U10 = (u_wind**2 + v_wind**2)**0.5
#drag = np.full(np.shape(U10), fill_value=np.NaN)
#drag[U10 <= 6.] = (0.29 + 3.1 / U10[U10 <= 6.] + 7.7 / (U10[U10 <= 6.]**2)) / 1000
drag = (0.6 + 0.071 * U10) / 1000

wind_stress_x = rho_air * drag * u_wind * abs(u_wind)
wind_stress_y = rho_air * drag * v_wind * abs(v_wind)

ekman_y = - (1/f_3d) * wind_stress_x  * dx_3d / (1025 * 1e6)
ekman_x = (1/f_3d) * wind_stress_y * dy / (1025 * 1e6)

#ekman_y[mean_dot > gyre_bounds] = np.NaN
#ekman_x[mean_dot > gyre_bounds] = np.NaN

ekman_y_bar = np.nanmean(ekman_y, axis=2)
ekman_x_bar = np.nanmean(ekman_x, axis=2)

ekman = (ekman_x ** 2 + ekman_y ** 2) ** 0.5
ekman_bar = (ekman_x_bar ** 2 + ekman_y_bar ** 2) ** 0.5

# Calculate the difference between winter and summer ekman transport
summer_x_ekman = np.full((59, 361, summer_count), fill_value=np.NaN)
summer_y_ekman = np.full((59, 361, summer_count), fill_value=np.NaN)
winter_x_ekman = np.full((59, 361, winter_count), fill_value=np.NaN)
winter_y_ekman = np.full((59, 361, winter_count), fill_value=np.NaN)
summer_step = 0
winter_step = 0
it = 0
# Calculate the summer/winter averages
for mnth in months:
    if mnth in summer_months:
        summer_x_ekman[:, :, summer_step] = ekman_x[:, :, it]
        summer_y_ekman[:, :, summer_step] = ekman_y[:, :, it]
        summer_step += 1
    if mnth in winter_months:
        winter_x_ekman[:, :, winter_step] = ekman_x[:, :, it]
        winter_y_ekman[:, :, winter_step] = ekman_y[:, :, it]
        winter_step += 1
    it += 1

# Calculate summer mean
summer_x_ekman_bar = np.nanmean(summer_x_ekman, axis=2)
summer_y_ekman_bar = np.nanmean(summer_y_ekman, axis=2)

# Calculate mean summer wind velocity
summer_ekman_bar = (summer_x_ekman_bar ** 2 + summer_y_ekman_bar ** 2) ** 0.5

# Calculate winter mean
winter_x_ekman_bar = np.nanmean(winter_x_ekman, axis=2)
winter_y_ekman_bar = np.nanmean(winter_y_ekman, axis=2)

# Calculate mean winter wind velocity
winter_ekman_bar = (winter_x_ekman_bar ** 2 + winter_y_ekman_bar ** 2) ** 0.5

if plotting==True:
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ekman_y_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-.05, .05])
    c.set_label('Ekman Transport (Sv)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/ekman_y_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ekman_x_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-.005, .005])
    c.set_label('Ekman Transport (Sv)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/ekman_x_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ekman_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([0, 0.05])
    c.set_label('Ekman Transport (Sv)')
    u_blank = ekman_x_bar
    v_blank = ekman_y_bar
    u_blank[mean_dot > gyre_bounds] = np.NaN
    u_blank[v_blank > 0.015] = np.NaN
    v_blank[mean_dot > gyre_bounds] = np.NaN
    u_blank[v_blank > 0.015] = np.NaN
    uproj,vproj,xx,yy = m.transform_vector(u_blank,v_blank,lon,lat, 25, 25,returnxy=True,masked=True)
    m.quiver(xx, yy, uproj, vproj, scale=0.1)
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/ekman_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(summer_y_ekman_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-.05, .05])
    c.set_label('Ekman Transport (Sv)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/ekman_summer_y_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(summer_x_ekman_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-.005, .005])
    c.set_label('Ekman Transport (Sv)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/ekman_summer_x_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_y_ekman_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-.05, .05])
    c.set_label('Ekman Transport (Sv)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/ekman_winter_y_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_x_ekman_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-.005, .005])
    c.set_label('Ekman Transport (Sv)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/ekman_winter_x_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_ekman_bar - summer_ekman_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-0.01, 0.01])
    c.set_label('Ekman Transport (Sv)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/ekman_winter-summer_anomaly.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
    
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_y_ekman_bar - summer_y_ekman_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-0.01, 0.01])
    c.set_label('Ekman Transport (Sv)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/ekman_winter-summer_y_anomaly.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
    
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_x_ekman_bar - summer_x_ekman_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([-0.01, 0.01])
    c.set_label('Ekman Transport (Sv)')
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/ekman_winter-summer_x_anomaly.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
    
    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(summer_ekman_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([0, 0.05])
    c.set_label('Ekman Transport (Sv)')
    u_blank = summer_x_ekman_bar
    v_blank = summer_y_ekman_bar
    u_blank[mean_dot > gyre_bounds] = np.NaN
    u_blank[v_blank > 0.015] = np.NaN
    v_blank[mean_dot > gyre_bounds] = np.NaN
    u_blank[v_blank > 0.015] = np.NaN
    uproj,vproj,xx,yy = m.transform_vector(u_blank,v_blank,lon,lat, 25, 25,returnxy=True,masked=True)
    m.quiver(xx, yy, uproj, vproj, scale=0.1)
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/ekman_summer_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=-80, urcrnrlon=55, urcrnrlat=-45, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(lat, lon)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(winter_ekman_bar)), cmap='RdBu_r')
    c = m.colorbar()
    pl.clim([0, 0.05])
    c.set_label('Ekman Transport (Sv)')
    u_blank = winter_x_ekman_bar
    v_blank = winter_y_ekman_bar
    u_blank[mean_dot > gyre_bounds] = np.NaN
    u_blank[v_blank > 0.015] = np.NaN
    v_blank[mean_dot > gyre_bounds] = np.NaN
    u_blank[v_blank > 0.015] = np.NaN
    uproj,vproj,xx,yy = m.transform_vector(u_blank,v_blank,lon,lat, 25, 25,returnxy=True,masked=True)
    m.quiver(xx, yy, uproj, vproj, scale=0.1)
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(mean_dot)), [gyre_bounds, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/WeddellGyre/Wind/ekman_winter_bar.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()