import os
import numpy as np
from mpl_toolkits.basemap import Basemap
import decimal
import scipy.ndimage

def month_data(directory, lead_offset=0.):
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
        lat_sub = []
        lon_sub = []
        ssh_sub = []
        type_sub = []
        ice_conc_sub = []
        lat_sub_asc = []
        lon_sub_asc = []
        ssh_sub_asc = []
        type_sub_asc = []
        ice_conc_sub_asc = []
        lat_sub_desc = []
        lon_sub_desc = []
        ssh_sub_desc = []
        type_sub_desc = []
        ice_conc_sub_desc = []
        f = open(file, 'r')
        for line in f:
            line = line.strip()
            columns = line.split()
            # If data point is listed as 'valid'
            if columns[1] == '1':
                # If data point is from open ocean (1) or from a lead (2)
                if columns[0] == '1' or columns[0] == '2':
                    # If the ssh is less than 3 m from the mean ssh
                    if float(columns[7]) - float(columns[8]) <= 3.:
                        lat_sub.append(float(columns[5]))
                        lon_sub.append(float(columns[6]))
                        # If the ssh is from a lead, apply the offset
                        if columns[0] == '2':
                            ssh_sub.append(float(columns[7]) + lead_offset)
                        # If the ssh is from the open ocean, don't apply offset
                        if columns[0] == '1':
                            ssh_sub.append(float(columns[7]))
                        ice_conc_sub.append(float(columns[11]))
                        type_sub.append(float(columns[0]))
        f.close()
        if len(lat_sub) > 3:
            # Do the DESCENDING tracks        
            descending = np.where(np.gradient(lat_sub) < 0.)[0]
            # If there are any descending tracks
            if len(descending) > 0.:
                ssh_sub_desc = ssh_sub[descending[0]:descending[-1]]
                lat_sub_desc = lat_sub[descending[0]:descending[-1]]
                lon_sub_desc = lon_sub[descending[0]:descending[-1]]
                ice_conc_sub_desc = ice_conc_sub[descending[0]:descending[-1]]
                type_sub_desc = type_sub[descending[0]:descending[-1]]
            
                bad_elements = []
                for issh in range(len(ssh_sub_desc)):
                    # If the value is greater than 3 std from the mean
                    if np.mean(ssh_sub_desc) - 3*np.std(ssh_sub_desc) > ssh_sub_desc[issh] > np.mean(ssh_sub_desc) + 3*np.std(ssh_sub_desc):
                        bad_elements.append(issh)
                for issh in range(len(ssh_sub_desc)-1):
                    # If the gradient between this point and the next point is greater than .5 m
                    if abs(ssh_sub_desc[issh] - ssh_sub_desc[issh + 1]) > .5:
                        bad_elements.append(issh + 1)
            
                # remove the points that meet the above criteria
                # In reverse order to avoid index problems
                for bad in sorted(np.unique(bad_elements), reverse=True):
                    del ssh_sub_desc[bad]
                    del lat_sub_desc[bad]
                    del lon_sub_desc[bad]
                    del ice_conc_sub_desc[bad]
                    del type_sub_desc[bad]
                lat += lat_sub_desc
                lon += lon_sub_desc
                type += type_sub_desc
                ice_conc += ice_conc_sub_desc
                # Apply a gaussian filter to the ssh data, with a 40 point (10 km) diameter
                #ssh += list(scipy.ndimage.filters.gaussian_filter1d(ssh_sub_desc, 40.))
                ssh += ssh_sub_desc
        
            ascending = np.where(np.gradient(lat_sub) > 0.)[0]
            # If there are any ascending tracks
            if len(ascending) > 0.:
                ssh_sub_asc = ssh_sub[ascending[1]:ascending[-1]]
                lat_sub_asc = lat_sub[ascending[1]:ascending[-1]]
                lon_sub_asc = lon_sub[ascending[1]:ascending[-1]]
                ice_conc_sub_asc = ice_conc_sub[ascending[1]:ascending[-1]]
                type_sub_asc = type_sub[ascending[1]:ascending[-1]]
        
                bad_elements = []
                # Do the Ascending tracks
                for issh in range(len(ssh_sub_asc)):
                    # If the value is greater than 3 std from the mean
                    if np.mean(ssh_sub_asc) - 3*np.std(ssh_sub_asc) > ssh_sub_asc[issh] > np.mean(ssh_sub_asc) + 3*np.std(ssh_sub_asc):
                        bad_elements.append(issh)
                for issh in range(len(ssh_sub_asc)-1):
                    # If the gradient between this point and the next point is greater than .5 m
                    if abs(ssh_sub_asc[issh] - ssh_sub_asc[issh + 1]) > .5:
                        bad_elements.append(issh + 1)
                
                for bad in sorted(np.unique(bad_elements), reverse=True):
                    del ssh_sub_asc[bad]
                    del lat_sub_asc[bad]
                    del lon_sub_asc[bad]
                    del ice_conc_sub_asc[bad]
                    del type_sub_asc[bad]
                lat += lat_sub_asc
                lon += lon_sub_asc
                type += type_sub_asc
                ice_conc += ice_conc_sub_asc
                # Apply a gaussian filter to the ssh data, with a 40 point (10 km) diameter
                #ssh += list(scipy.ndimage.filters.gaussian_filter1d(ssh_sub_asc, 40.))
                ssh += ssh_sub_asc

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
    
    x_range = np.arange(np.round(np.min(x_cyl)), np.round(np.max(x_cyl))+ lon_res, lon_res)
    y_range = np.arange(np.round(np.min(y_cyl)), np.round(np.max(y_cyl))+ lat_res, lat_res)
    xy_grid = np.full([np.size(x_range), np.size(y_range)], fill_value=np.nan)
    xy_count = np.full([np.size(x_range), np.size(y_range)], fill_value=np.nan)

    for i in range(np.size(data)):
        x_coord = float(decimal.Decimal(float(x_cyl[i])).quantize(decimal.Decimal(str(lon_res)), rounding='ROUND_DOWN'))
        y_coord = float(decimal.Decimal(float(y_cyl[i])).quantize(decimal.Decimal(str(lat_res)), rounding='ROUND_DOWN'))

        ix_range = np.where(np.logical_and(x_range > x_coord - lon_res/2, x_range < x_coord + lon_res/2))
        iy_range = np.where(np.logical_and(y_range > y_coord - lat_res/2, y_range < y_coord + lat_res/2))

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

        ix_range = np.where(np.logical_and(x_range > x_coord - res/2, x_range < x_coord + res/2))
        iy_range = np.where(np.logical_and(y_range > y_coord - res/2, y_range < y_coord + res/2))

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

def ocean_lead_offset(month):
    if month == '01':
        return 0.053
    if month == '02':
        return 0.042
    if month == '03':
        return 0.044
    if month == '04':
        return 0.059
    if month == '05':
        return 0.062
    if month == '06':
        return 0.089
    if month == '07':
        return 0.108
    if month == '08':
        return 0.105
    if month == '09':
        return 0.092
    if month == '10':
        return 0.096
    if month == '11':
        return 0.059
    if month == '12':
        return 0.047
#January offset:  0.0532795957319
#Febuary offset:  0.0424895753339
#March offset:  0.044196774024
#April offset:  0.058855642473
#May offset:  0.0623327136638
#June offset:  0.0885369699514
#July offset:  0.108297693989
#August offset:  0.104956289149
#September offset:  0.0924955874413
#October offset:  0.0964708361887
#November offset:  0.0591056776193
#December offset:  0.0469972667101
