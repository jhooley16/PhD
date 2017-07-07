import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, Point, LinearRing, LineString
from geopy.distance import vincenty
import functions as funct

dot_seasonal = np.full((59, 361, 7, 12), fill_value=np.NaN)
ice = np.full((59, 361, 7, 12), fill_value=np.NaN)
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year)
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if os.path.isfile(year + month + '_DOT_filt.nc'):
            nc = Dataset(year + month + '_DOT_filt.nc', 'r')
            dot_seasonal[:, :, int(year)-2010, int(month)-1] = nc.variables['dynamic_ocean_topography_seasonal_offset'][:]
            ice[:, :, int(year)-2010, int(month)-1] = nc.variables['sea_ice_concentration'][:]
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            nc.close()

# Extract nodes from 'crude' Antarctica polygon in Basemap
m = Basemap(projection='spstere', boundinglat=-60, lon_0=180, resolution='c')
coordinates = np.vstack(m.drawcoastlines().get_segments())
# Convert these to longitude/latitude coords from projection coords
lons,lats = m(coordinates[:,0],coordinates[:,1],inverse=True)

# Only take coastline data from South of 60S
AA_lats = []
AA_lons = []
for it in range(len(lats)):
    if lats[it] < -60:
        AA_lats.append(lats[it])
        AA_lons.append(lons[it])

# Convert these back to projection coordinates
AA_stereo_x, AA_stereo_y = m(AA_lons, AA_lats)

# Pair them up for use as polygon coordinates
AA_paired = []
for it in range(len(AA_stereo_x)):
    AA_paired.append((AA_stereo_x[it], AA_stereo_y[it]))

# Define the AA polygon as a LineString for measuring distance along, etc
linepoly = LineString(AA_paired)
# Generate a polygon object based on the new AA coordinates
AA_poly = Polygon(AA_paired)

# measure the total distance around Antarctica
AA_distance = 0.
for ipair in range(len(AA_paired) - 1):
    AA_distance += Point(AA_paired[ipair]).distance(Point(AA_paired[ipair + 1]))

# Interpolate around Polygon with equal distance between points
testpoint = []
for ix in range(1, int(AA_distance), 100000):
    testpoint.append(linepoly.interpolate(ix).xy)


# Get the projection (x, y) coordinates of the interpolated polygon points
testpoints_x = []
testpoints_y = []
for it in range(len(testpoint)):
    # Subtract South Pole coordinate from each point
    testpoints_x.append(tuple(testpoint[it])[0][0])
    testpoints_y.append(tuple(testpoint[it])[1][0])

pl.figure()
pl.clf()
m.drawcoastlines()
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])

m.scatter(AA_stereo_x, AA_stereo_y)
m.scatter(np.array(testpoints_x), np.array(testpoints_y), color='g')
m.scatter(testpoints_x[0], testpoints_y[0], color='r')
m.scatter(testpoints_x[-1], testpoints_y[-1], color='k')

pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/experiment.png', format='png', transparent=True, dpi=300)
pl.close()

######## Test making an orthogonal line from a two-point line #####
for ix in range(len(testpoints_x)):
    print(ix)
    if 2 <= ix <= len(testpoints_x)-3:
        min = ix-2
        max = ix+2
    
    elif ix <= 1:
        max = ix+2
        if ix == 0:
            min = len(testpoints_x)-2
        elif ix == 1:
            min = len(testpoints_x)-1

        
    elif ix >= len(testpoints_x)-2:
        min = ix-2
        if ix == len(testpoints_x)-2:
            max = 0
        elif ix == len(testpoints_x)-1:
            max = 1
    
    A_x = testpoints_x[min]
    A_y = testpoints_y[min]
    B_x = testpoints_x[ix]
    B_y = testpoints_y[ix]
    C_x = testpoints_x[max]
    C_y = testpoints_y[max]

    # Get equation of line AC
    AC_m = (C_y-A_y) / (C_x-A_x)
    AC_c = C_y - AC_m*C_x

    # Get equation of perpendicular line DB
    DB_m = -1/AC_m
    DB_c = B_y - DB_m*B_x

    # Find where two lines intersect
#     D_x = (AC_c-DB_c) / (DB_m-AC_m)
#     D_y = DB_m*D_x + DB_c

    
    B_lon, B_lat = m(B_x, B_y, inverse=True)
    
    test_dist = False
    A = 0
    inside = []
    while test_dist == False:
        E_x_test = B_x + A
        E_y_test = DB_m*E_x_test + DB_c
        E_lon_test, E_lat_test = m(E_x_test, E_y_test, inverse=True)
        if vincenty((B_lat, B_lon), (E_lat_test, E_lon_test)).km > 50:
            test_dist = True
            if Point(E_x_test, E_y_test).within(AA_poly) == True:
                print('Point within Coastline (WRONG)')
                off_coast = False
            elif Point(E_x_test, E_y_test).within(AA_poly) == False:
                print('Point off Coastline (CORRECT)')
                off_coast = True
        A += 50

    good_dist = False
    A = 0
    if off_coast == True:
        print('1')
        while good_dist == False:
            E_x = B_x + A
            E_y = DB_m*E_x + DB_c
            E_lon, E_lat = m(E_x, E_y, inverse=True)
            if vincenty((B_lat, B_lon), (E_lat, E_lon)).km > 1000:
                good_dist = True
            A += 50
        if Point(E_x, E_y).within(AA_poly) == True:
            print('1_2')
            good_dist = False
            A = 0
            while good_dist == False:
                E_x = B_x + A
                E_y = DB_m*E_x + DB_c
                E_lon, E_lat = m(E_x, E_y, inverse=True)
                if vincenty((B_lat, B_lon), (E_lat, E_lon)).km > 1000:
                    good_dist = True
                A -= 50
            
    elif off_coast == False:
        print('2')
        while good_dist == False:
            E_x = B_x + A
            E_y = DB_m*E_x + DB_c
            E_lon, E_lat = m(E_x, E_y, inverse=True)
            if vincenty((B_lat, B_lon), (E_lat, E_lon)).km > 1000:
                good_dist = True
            A -= 50
        if Point(E_x, E_y).within(AA_poly) == True:
            print('2_2')
            good_dist = False
            A = 0
            while good_dist == False:
                E_x = B_x + A
                E_y = DB_m*E_x + DB_c
                E_lon, E_lat = m(E_x, E_y, inverse=True)
                if vincenty((B_lat, B_lon), (E_lat, E_lon)).km > 1000:
                    good_dist = True
                A += 50

    if ix < 10:
        STRING = '00' + str(ix)
    elif 9 < ix < 100:
        STRING = '0' + str(ix)
    elif 99 < ix:
        STRING = str(ix)
    
    pl.figure()
    pl.clf()
    m.drawcoastlines()
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    m.scatter(E_x, E_y, color='r')
    #m.scatter(AA_stereo_x, AA_stereo_y)
    m.scatter(np.array(testpoints_x), np.array(testpoints_y), color='g')

    m.plot([A_x, C_x], [A_y, C_y], color='r')
    m.plot([B_x, E_x], [B_y, E_y], color='k')
    pl.title(ix)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/gif/experiment_2_'+STRING+'.png', format='png', transparent=True)
    pl.close()

funct.gifmaker('test', '/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/gif/')
