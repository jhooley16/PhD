import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, Point, LineString
from geopy.distance import vincenty
import functions as funct

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

# Pair them up as tuples for use as polygon coordinates
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
points_x = []
points_y = []
for it in range(len(testpoint)):
    points_x.append(tuple(testpoint[it])[0][0])
    points_y.append(tuple(testpoint[it])[1][0])

# 'rotate' the points so that points_x and points_y start at the tip
# of the Antarctic Peninsula
points_x = points_x[192:] + points_x[:192]
points_y = points_y[192:] + points_y[:192]

# pl.figure()
# pl.clf()
# m.drawcoastlines()
# m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
# m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
# m.scatter(np.array(points_x), np.array(points_y), color='g')
# m.scatter(points_x[0], points_y[0], color='r')
# m.scatter(points_x[-1], points_y[-1], color='k')
# pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/experiment.png', format='png', transparent=True, dpi=300)
# pl.close()

# Make an orthogonal line from a two-point line
for ix in range(len(points_x)):
    # Take the points two steps ahead and two steps behind
    if 2 <= ix <= len(points_x) - 3:
        min = ix - 2
        max = ix + 2
    # Define these points at the start/end of the list of points
    # If we're at the start...
    elif ix <= 1:
        max = ix + 2
        if ix == 0:
            min = len(points_x) - 2
        elif ix == 1:
            min = len(points_x) - 1
    # If we're at the end...
    elif ix >= len(points_x) - 2:
        min = ix - 2
        if ix == len(points_x) - 2:
            max = 0
        elif ix == len(points_x) - 1:
            max = 1

    # Point A is the before point
    A_x = points_x[min]
    A_y = points_y[min]
    # Point B is the current point
    B_x = points_x[ix]
    B_y = points_y[ix]
    # Point C is the after point
    C_x = points_x[max]
    C_y = points_y[max]

    # Get equation of line AC
    AC_m = (C_y - A_y) / (C_x - A_x)
    AC_c = C_y - AC_m * C_x

    # Get equation of perpendicular line DB
    # Point D is the intersection point on line AC and 
    # its normal, which intersects point B
    DB_m = -1 / AC_m
    DB_c = B_y - DB_m * B_x
    
    # Convert the coordinates of B into lon/lat coordinates
    B_lon, B_lat = m(B_x, B_y, inverse=True)
    
    # Calculate points along DB, increasing the distance slightly each loop
    # until the distance is 50 km from point B
    # If the new point (E) is inside/outside the Antarctic polygon (and 
    # therefore on/off land (Wrong), return True/False.
    test_dist = False
    A = 0
    inside = []
    while test_dist == False:
        E_x_test = B_x + A
        E_y_test = DB_m * E_x_test + DB_c
        E_lon_test, E_lat_test = m(E_x_test, E_y_test, inverse=True)
        if vincenty((B_lat, B_lon), (E_lat_test, E_lon_test)).km > 50:
            test_dist = True
            if Point(E_x_test, E_y_test).within(AA_poly) == True:
                off_coast = False
            elif Point(E_x_test, E_y_test).within(AA_poly) == False:
                off_coast = True
        A += 50

    # If point E is off the coast (and therefore in the ocean) then 
    # run the loop again but increase distance to 1000 km from coast
    good_dist = False
    A = 0
    if off_coast == True:
        while good_dist == False:
            E_x = B_x + A
            E_y = DB_m * E_x + DB_c
            E_lon, E_lat = m(E_x, E_y, inverse=True)
            if vincenty((B_lat, B_lon), (E_lat, E_lon)).km > 1000:
                good_dist = True
            A += 50
        # Just in case: Check if the new point E is inside the coastline
        # If it is, run loop in reverse and that should fix it
        if Point(E_x, E_y).within(AA_poly) == True:
            good_dist = False
            A = 0
            while good_dist == False:
                E_x = B_x + A
                E_y = DB_m * E_x + DB_c
                E_lon, E_lat = m(E_x, E_y, inverse=True)
                if vincenty((B_lat, B_lon), (E_lat, E_lon)).km > 1000:
                    good_dist = True
                A -= 50
    # If point E is on the coast, run the loop in reverse
    elif off_coast == False:
        while good_dist == False:
            E_x = B_x + A
            E_y = DB_m * E_x + DB_c
            E_lon, E_lat = m(E_x, E_y, inverse=True)
            if vincenty((B_lat, B_lon), (E_lat, E_lon)).km > 1000:
                good_dist = True
            A -= 50
        # Just in case, Check if E is inside coastline, if it is, run 
        # loop again
        if Point(E_x, E_y).within(AA_poly) == True:
            good_dist = False
            A = 0
            while good_dist == False:
                E_x = B_x + A
                E_y = DB_m * E_x + DB_c
                E_lon, E_lat = m(E_x, E_y, inverse=True)
                if vincenty((B_lat, B_lon), (E_lat, E_lon)).km > 1000:
                    good_dist = True
                A += 50
    
    # Construct figure string number
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
    m.scatter(B_x, B_y, color='k')
    m.scatter(np.array(points_x), np.array(points_y), color='g')
    m.plot([A_x, C_x], [A_y, C_y], color='r')
    m.plot([B_x, E_x], [B_y, E_y], color='k')
    pl.title(STRING)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/dancing_matchstick/experiment_2_'+STRING+'.png', format='png', transparent=True)
    pl.close()

funct.gifmaker('dancing_matchstick', '/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/dancing_matchstick/')
