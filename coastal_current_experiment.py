import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, Point, LinearRing, LineString

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

AA_lats = []
AA_lons = []
for it in range(len(lats)):
    if lats[it] < -60:
        AA_lats.append(lats[it])
        AA_lons.append(lons[it])

AA_stereo_x, AA_stereo_y = m(AA_lons, AA_lats)

AA_paired = []
for it in range(len(AA_stereo_x)):
    AA_paired.append((AA_stereo_x[it], AA_stereo_y[it]))

# Define an arbitrary point
point_x = AA_stereo_x[0] + 3000000
point_y = AA_stereo_y[0] - 200000

# Define the AA polygon as a LineString for measuring distance, etc
linepoly = LineString(AA_paired)

# measure the total distance around Antarctica
distance = 0.
for ipair in range(len(AA_paired) - 1):
    distance += Point(AA_paired[ipair]).distance(Point(AA_paired[ipair + 1]))

# Interpolate around Polygon with equal distance between points
testpoint = []
for ix in range(1, int(distance), 100000):
    testpoint.append(linepoly.interpolate(ix).xy)

testpoints_x = []
testpoints_y = []
for it in range(len(testpoint)):
    testpoints_x.append(tuple(testpoint[it])[0][0])
    testpoints_y.append(tuple(testpoint[it])[1][0])

#poly_2 = Polygon(testpoints)
point = Point(point_x, point_y)

# pol_ext = LinearRing(poly.exterior.coords)
# d = pol_ext.project(point)
# p = pol_ext.interpolate(d)
# closest_point_coords = list(p.coords)[0]
# 
# closest_point_x = closest_point_coords[0]
# closest_point_y = closest_point_coords[1]

pl.figure()
pl.clf()
m.drawcoastlines()
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])

m.scatter(AA_stereo_x, AA_stereo_y)
m.scatter(testpoints_x, testpoints_y, color='g')
m.scatter(testpoints_x[0], testpoints_y[0], color='r')
m.scatter(testpoints_x[-1], testpoints_y[-1], color='k')

# m.scatter(point_x, point_y, color='r')
# m.scatter(closest_point_x, closest_point_y, color='g')


pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/experiment.png', format='png', transparent=True, dpi=300)
pl.close()

######## Test making an orthogonal line from a two-point line #####

testing_points_x = testpoints_x[:3]
testing_points_y = testpoints_y[:3]

A_x = testing_points_x[0]
A_y = testing_points_y[0]

B_x = testing_points_x[1]
B_y = testing_points_y[1]

C_x = testing_points_x[2]
C_y = testing_points_y[2]

t = - (C_x-A_x)*(B_x-A_x) + (C_y-A_y)*(B_y-A_y)
t /= (C_x - A_x)**2 + (C_y - A_y)**2

D_x = A_x + t*(C_x - A_x)
D_y = A_y + t*(C_y - A_y)

pl.figure()
pl.clf()

pl.scatter(A_x, A_y, label='Point A', color='b')
pl.scatter(B_x, B_y, label='Point B', color='r')
pl.scatter(C_x, C_y, label='Point C', color='g')
pl.scatter(D_x, D_y, label='Point D', color='k')

pl.plot([A_x, C_x], [A_y, C_y], label='Line AC')

pl.legend()

pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/experiment_line.png', format='png', transparent=True, dpi=300)
pl.close()