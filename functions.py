import os
import numpy as np
from mpl_toolkits.basemap import Basemap
import decimal

def month_data(directory):
    """Extract latitude, longitude, sea surface height, and surface ice concentration
    data from all files in 'directory'.
    eg: directory = '/Users/jmh2g09/Desktop/201203_MERGE'"""
    files = os.listdir(directory)
    lat = []
    lon = []
    ssh = []
    type = []
    ice_conc = []
    for file in files:
        f = open(file, 'r')
        for line in f:
            line = line.strip()
            columns = line.split()
            # If data point is from open ocean (1) or from a lead (2)
            if columns[0] == '1' or columns[0] == '2':
                # If data point is listed as 'valid'
                if columns[1] == '1':
                    # If the ssh is less than 3m from the mean ssh
                    #if float(columns[7]) - float(columns[8]) <= 3.:
                    lat.append(float(columns[5]))
                    lon.append(float(columns[6]))
                    ssh.append(float(columns[7]))
                    ice_conc.append(float(columns[11]))
                    type.append(float(columns[0]))
        f.close()
    return {'lat': lat, 'lon': lon, 'ssh': ssh, 'ice_conc': ice_conc, 'type': type}

def day_data(day, directory):
    """Extract latitude, longitude, sea surface height, and surface ice concentration
    data for a certain day in 'directory'. The data is taken from a certain area 
    defined by the latitude and longitude limits. See the code for details.
    eg: directory = '/Users/jmh2g09/Desktop/201203_MERGE'
    day corresponds to the day in a particular month."""
    files = os.listdir(directory)
    if 10 <= day <= 31:
        dayfiles = [file for file in files if file[13:15] == str(day)]
    if day < 10:
        dayfiles = [file for file in files if file[14] == str(day)]
    lat = []
    lon = []
    ssh = []
    ice_conc = []
    for file in dayfiles:
        f = open(file, 'r')
        for line in f:
            line = line.strip()
            columns = line.split()
            # If data point is from open ocean [1] or from a lead [2]
            if columns[0] == '1' or columns[0] == '2':
                # If data point is listed as 'valid'
                if columns[1] == '1':
                    # If the ssh is less than 3m from the mean ssh
                    if float(columns[7]) - float(columns[8]) <= 3.:
                        # If the latitude lies between
                        if -60. > float(columns[5]) > -62.:
                            # If the longitude lies between
                            if -54. > float(columns[6]) > -56.:
                                lat.append(float(columns[5]))
                                lon.append(float(columns[6]))
                                ssh.append(float(columns[7]))
                                ice_conc.append(float(columns[11]))
        f.close()
    return {'lat': lat, 'lon': lon, 'ssh': ssh, 'ice_conc': ice_conc}


def stereo(lat, lon):
    """Convert Lat and Lon (Polar) coordinates to stereographic (x, y) coordinates (km).
        Latitude in degrees.
        Longitude in degrees.
        sgn = 1 for Northern Hempsphere.
        sgn = -1 for Southern Hemisphere."""
    def general(angle):
        epsilon = 0.081816153  # Eccentricity of ellipsoid (Earth)
        T = (np.tan((np.pi/4) - (angle/2)) /
                ((1 - epsilon * np.sin(angle)) /
                    (1 + epsilon * np.sin(angle))) ** epsilon / 2)
        return T
    PHI = 70. * np.pi / 180  # Projection plane in radians (70 degrees)
    theta = np.absolute(lat) * np.pi / 180  # Latitude in radians
    phi = lon * (np.pi / 180)  # Longitude in radians
    epsilon = 0.081816153  # Eccentricity of ellipsoid (Earth)
    R = 6378.273  # Radius of ellipsoid in km (Earth)
    gamma = np.cos(PHI) / np.sqrt(1 - (epsilon**2) * np.sin(PHI)**2)
    rho = R * gamma * general(theta) / general(PHI)
    x = rho * np.sin(phi)  # stereographic x-coordinate in km
    y = rho * np. cos(phi)  # stereographic y-coordinate in km
    return {'X': x, 'Y': y}


def grid05(data, lon_data, lat_data, lat_res, lon_res):
    """A function to grid lon, lat data and produce a masked array.
    
    The data should be of the form (x, y, z), where x is the longitude, y is the latitude
    and z is the value at that position. data, lon_data, and lat_data must all have the  
    same length, and be one dimensional. res is the desired horizontal resolution, for example:
    for a 1-degree grid, use res = 1, for 0.5-degree grid, use res = 0.5.
    
    The gridded data is drawn on an equidistant cylindrical projection. Use your favourite
    projection conversion tool to convert the result to your desired projection.
    
    This simple function draws a square grid on top of the data distribution. Any data
    points that lie within a grid square are averaged and the number of data points
    averaged in that box is returned in xy_count."""
    
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=-30, llcrnrlon=-180,
                urcrnrlon=180, resolution='c')
    x_cyl, y_cyl = m(lon_data, lat_data)
    
    x_range = np.arange(np.round(np.min(x_cyl)) - 0.5, np.round(np.max(x_cyl)) + 0.5 + lon_res, lon_res)
    y_range = np.arange(np.round(np.min(y_cyl)) - 0.5, np.round(np.max(y_cyl)) + 0.5 + lat_res, lat_res)
    xy_grid = np.full([np.size(x_range), np.size(y_range)], fill_value=np.nan)
    xy_count = np.full([np.size(x_range), np.size(y_range)], fill_value=np.nan)

    for i in range(np.size(data)):
        x_coord = float(decimal.Decimal(float(x_cyl[i])).quantize(decimal.Decimal(str(lon_res)), rounding='ROUND_DOWN')) + 0.5
        y_coord = float(decimal.Decimal(float(y_cyl[i])).quantize(decimal.Decimal(str(lat_res)), rounding='ROUND_DOWN')) + 0.5

        ix_range = np.where(np.logical_and(x_range >= x_coord - lon_res/2, x_range <= x_coord + lon_res/2))
        iy_range = np.where(np.logical_and(y_range >= y_coord - lat_res/2, y_range <= y_coord + lat_res/2))

        if np.isnan(xy_grid[ix_range, iy_range]):
            xy_grid[ix_range, iy_range] = data[i]
            xy_count[ix_range, iy_range] = 1

        if not np.isnan(xy_grid[ix_range, iy_range]):
            xy_grid[ix_range, iy_range] = xy_grid[ix_range, iy_range] + data[i]
            xy_count[ix_range, iy_range] = xy_count[ix_range, iy_range] + 1

    xy_grid = xy_grid / xy_count

    xy_grid = np.ma.masked_where(np.isnan(xy_grid), xy_grid)
    xy_count = np.ma.masked_where(np.isnan(xy_count), xy_count)

    return {'Grid':xy_grid, 'Count':xy_count, 'Lon':x_range, 'Lat':y_range}

def grid(data, lon_data, lat_data, res):
    """A function to grid lon, lat data and produce a masked array.
    
    The data should be of the form (x, y, z), where x is the longitude, y is the latitude
    and z is the value at that position. data, lon_data, and lat_data must all have the  
    same length, and be one dimensional. res is the desired horizontal resolution, for example:
    for a 1-degree grid, use res = 1, for 0.5-degree grid, use res = 0.5.
    
    The gridded data is drawn on an equidistant cylindrical projection. Use your favourite
    projection conversion tool to convert the result to your desired projection.
    
    This simple function draws a square grid on top of the data distribution. Any data
    points that lie within a grid square are averaged and the number of data points
    averaged in that box is returned in xy_count."""

    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=-30, llcrnrlon=-180,
                urcrnrlon=180, resolution='c')
    x_cyl, y_cyl = m(lon_data, lat_data)

    x_range = np.arange(np.round(np.min(x_cyl)), np.round(np.max(x_cyl)) + res, res)
    y_range = np.arange(np.round(np.min(y_cyl)), np.round(np.max(y_cyl)) + res, res)
    xy_grid = np.full([np.size(x_range), np.size(y_range)], fill_value=np.nan)
    xy_count = np.full([np.size(x_range), np.size(y_range)], fill_value=np.nan)

    for i in range(np.size(data)):
        x_coord = float(decimal.Decimal(float(x_cyl[i])).quantize(decimal.Decimal(str(res))))
        y_coord = float(decimal.Decimal(float(y_cyl[i])).quantize(decimal.Decimal(str(res))))

        ix_range = np.where(np.logical_and(x_range >= x_coord - res/2, x_range <= x_coord + res/2))
        iy_range = np.where(np.logical_and(y_range >= y_coord - res/2, y_range <= y_coord + res/2))

        if np.isnan(xy_grid[ix_range, iy_range]):
            xy_grid[ix_range, iy_range] = data[i]
            xy_count[ix_range, iy_range] = 1

        if not np.isnan(xy_grid[ix_range, iy_range]):
            xy_grid[ix_range, iy_range] = xy_grid[ix_range, iy_range] + data[i]
            xy_count[ix_range, iy_range] = xy_count[ix_range, iy_range] + 1

    xy_grid = xy_grid / xy_count

    xy_grid = np.ma.masked_where(np.isnan(xy_grid), xy_grid)
    xy_count = np.ma.masked_where(np.isnan(xy_count), xy_count)

    return {'Grid':xy_grid, 'Count':xy_count, 'Lon':x_range, 'Lat':y_range}

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2-1):-(window_len/2)]