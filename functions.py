import os
import numpy as np
from mpl_toolkits.basemap import Basemap
import decimal
import scipy.ndimage
from netCDF4 import Dataset
import matplotlib.pyplot as pl
import os
import imageio

def month_data(directory, month):
    surface = []
    time = []
    lat = []
    lon = []
    ssh = []
    ice_conc = []
    files = os.listdir(directory)
    for file in files:
        print(file)
        surface_sub = []
        time_sub = []
        lat_sub = []
        lon_sub = []
        ssh_sub = []
        ice_conc_sub = []
    
        surface_sub_asc = []
        time_sub_asc = []
        lat_sub_asc = []
        lon_sub_asc = []
        ssh_sub_asc = []
        ice_conc_sub_asc = []
    
        surface_sub_desc = []
        time_sub_desc = []
        lat_sub_desc = []
        lon_sub_desc = []
        ssh_sub_desc = []
        ice_conc_sub_desc = []
        f = open(file, 'r')
        for line in f:
            line = line.strip()
            columns = line.split()
            # If data point is listed as 'valid'
            if columns[1] == '1':
                # If data point is from open ocean (1) or from a lead (2)
                if columns[0] == '1' or columns[0] == '2':
                    # If the ssh is less than 0.3 m from the mean ssh
                    if abs(float(columns[7]) - float(columns[8])) <= .3:
                        surface_sub.append(float(columns[0]))
                        time_sub.append(float(columns[4]))
                        lat_sub.append(float(columns[5]))
                        lon_sub.append(float(columns[6]))
                        ssh_sub.append(float(columns[7]))
                        ice_conc_sub.append(float(columns[11]))
        f.close()

        if len(lat_sub) > 3:
            ###### DESCENDING tracks ########
            descending = np.where(np.gradient(lat_sub) < 0.)[0]
            # If there are any descending tracks
            if len(descending) > 0.:
                surface_sub_desc = surface_sub[descending[0]:descending[-1]]
                time_sub_desc = time_sub[descending[0]:descending[-1]]
                lat_sub_desc = lat_sub[descending[0]:descending[-1]]
                lon_sub_desc = lon_sub[descending[0]:descending[-1]]
                ssh_sub_desc = ssh_sub[descending[0]:descending[-1]]
                ice_conc_sub_desc = ice_conc_sub[descending[0]:descending[-1]]

                bad_elements = []
                for issh in range(len(ssh_sub_desc)):
                    # If the value is greater than 3 std from the mean
                    if np.mean(ssh_sub_desc) - 3*np.std(ssh_sub_desc) > ssh_sub_desc[issh] > np.mean(ssh_sub_desc) + 3*np.std(ssh_sub_desc):
                        bad_elements.append(issh)
                for issh in range(len(ssh_sub_desc) - 1):
                    # If the gradient between this point and the next point is greater than .5 m
                    if abs(ssh_sub_desc[issh] - ssh_sub_desc[issh + 1]) > .25:
                        bad_elements.append(issh + 1)
        
                # remove the points that meet the above criteria
                # In reverse order to avoid index problems
                for bad in sorted(np.unique(bad_elements), reverse=True):
                    del surface_sub_desc[bad]
                    del time_sub_desc[bad]
                    del lat_sub_desc[bad]
                    del lon_sub_desc[bad]
                    del ssh_sub_desc[bad]
                    del ice_conc_sub_desc[bad]

                ## Filter the track
                input_ssh = open('../INPUT_ssh.dat', 'w')
                input_ice = open('../INPUT_ice.dat', 'w')
                for ilat in range(len(lat_sub_desc)):
                    print(-lat_sub_desc[ilat], ssh_sub_desc[ilat], file=input_ssh)
                    print(-lat_sub_desc[ilat], ice_conc_sub_desc[ilat], file=input_ice)
                input_ssh.close()
                input_ice.close()
                
#                 pl.figure()
#                 
#                 pl.plot(lat_sub_desc, ssh_sub_desc, label='Unfiltered Data')
                
                os.system('gmt filter1d ../INPUT_ssh.dat -Fg0.2 -D0.001 -fi0y -E > ../OUTPUT_ssh.dat')
                os.system('gmt filter1d ../INPUT_ice.dat -Fg0.2 -D0.001 -fi0y -E > ../OUTPUT_ice.dat')
            
                os.system('rm ../INPUT_ssh.dat ../INPUT_ice.dat')
            
                output_ssh = open('../OUTPUT_ssh.dat', 'r')
                lat_sub_desc = []
                ssh_sub_desc = []
                for line in output_ssh:
                    line.strip()
                    columns = line.split()
                    lat_sub_desc.append(-float(columns[0]))
                    ssh_sub_desc.append(float(columns[1]))
                output_ssh.close()
            
                output_ice = open('../OUTPUT_ice.dat', 'r')
                ice_conc_sub_desc = []
                for line in output_ice:
                    line.strip()
                    columns = line.split()
                    ice_conc_sub_desc.append(float(columns[1]))
                output_ice.close()
                
#                 pl.plot(lat_sub_desc, ssh_sub_desc, 'r', label='Filtered Data')
#                 pl.legend(loc='best')
#                 pl.ylabel('Sea Surface Height (m)')
#                 pl.xlabel('Latitude')
#                 pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Processed/Figures/filtered_along_track_example.png', doi=300, transparent=True, bbox_inches='tight')
#                 pl.close()
            
                os.system('rm ../OUTPUT_ssh.dat ../OUTPUT_ice.dat')
#                 pause
                if len(lat_sub_desc) == len(lon_sub_desc):
                    surface += surface_sub_desc
                    time += time_sub_desc
                    lat += lat_sub_desc
                    lon += lon_sub_desc
                    ssh += ssh_sub_desc
                    ice_conc += ice_conc_sub_desc
            ###### ASCENDING tracks ########
            ascending = np.where(np.gradient(lat_sub) > 0.)[0]
            # If there are any ascending tracks
            if len(ascending) > 0.:
                surface_sub_asc = surface_sub[ascending[1]:ascending[-1]]
                time_sub_asc = time_sub[ascending[1]:ascending[-1]]
                lat_sub_asc = lat_sub[ascending[1]:ascending[-1]]
                lon_sub_asc = lon_sub[ascending[1]:ascending[-1]]
                ssh_sub_asc = ssh_sub[ascending[1]:ascending[-1]]
                ice_conc_sub_asc = ice_conc_sub[ascending[1]:ascending[-1]]

                bad_elements = []
                # Do the Ascending tracks
                for issh in range(len(ssh_sub_asc)):
                    # If the value is greater than 3 std from the mean
                    if np.mean(ssh_sub_asc) - 3*np.std(ssh_sub_asc) > ssh_sub_asc[issh] > np.mean(ssh_sub_asc) + 3*np.std(ssh_sub_asc):
                        bad_elements.append(issh)
                for issh in range(len(ssh_sub_asc) - 1):
                    # If the gradient between this point and the next point is greater than .5 m
                    if abs(ssh_sub_asc[issh] - ssh_sub_asc[issh + 1]) > .5:
                        bad_elements.append(issh + 1)
            
                for bad in sorted(np.unique(bad_elements), reverse=True):
                    del surface_sub_asc[bad]
                    del time_sub_asc[bad]
                    del lat_sub_asc[bad]
                    del lon_sub_asc[bad]
                    del ssh_sub_asc[bad]
                    del ice_conc_sub_asc[bad]

                ## Filter the track
                input_ssh = open('../INPUT_ssh.dat', 'w')
                input_ice = open('../INPUT_ice.dat', 'w')
                for ilat in range(len(lat_sub_asc)):
                    print(lat_sub_asc[ilat], ssh_sub_asc[ilat], file=input_ssh)
                    print(lat_sub_asc[ilat], ice_conc_sub_asc[ilat], file=input_ice)
                input_ssh.close()
                input_ice.close()
            
                os.system('gmt filter1d ../INPUT_ssh.dat -Fg0.2 -D0.001 -fi0y -E > ../OUTPUT_ssh.dat')
                os.system('gmt filter1d ../INPUT_ice.dat -Fg0.2 -D0.001 -fi0y -E > ../OUTPUT_ice.dat')
            
                os.system('rm ../INPUT_ssh.dat ../INPUT_ice.dat')
            
                output_ssh = open('../OUTPUT_ssh.dat', 'r')
                lat_sub_asc = []
                ssh_sub_asc = []
                for line in output_ssh:
                    line.strip()
                    columns = line.split()
                    lat_sub_asc.append(float(columns[0]))
                    ssh_sub_asc.append(float(columns[1]))
                output_ssh.close()
            
                output_ice = open('../OUTPUT_ice.dat', 'r')
                ice_conc_sub_asc = []
                for line in output_ice:
                    line.strip()
                    columns = line.split()
                    ice_conc_sub_asc.append(float(columns[1]))
                output_ice.close()
            
                os.system('rm ../OUTPUT_ssh.dat ../OUTPUT_ice.dat')
            
                if len(lat_sub_asc) == len(lon_sub_asc):
                    surface += surface_sub_asc
                    time += time_sub_asc
                    lat += lat_sub_asc
                    lon += lon_sub_asc
                    ssh += ssh_sub_asc
                    ice_conc += ice_conc_sub_asc
                    
    # Calculate the CS-2 mode with which the data points were measured
    mode = mode_points(lat, lon, month)

    return {'lat': lat, 'lon': lon, 'ssh': ssh, 'ice_conc': ice_conc, 'surface': surface, 'time': time, 'mode': mode}

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

def surface_area(lat, lon, cell_size_lat, cell_size_lon):
    """A function to calculate the surface area of a grid cell on a sphere.
        Grid cell has sides (cell_size_lat, cell_size_lon) in degrees.
        Returns the surface area of the grid cell at location (lat, lon) in km^2."""
    # Radius of the Earth (km)
    R = 6371.
    # Step 1, get limits of the grid cell in radians
    lat_side_min = (lat - (cell_size_lat / 2)) * np.pi / 180
    lat_side_max = (lat + (cell_size_lat / 2)) * np.pi / 180

    lon_side_min = (lon - (cell_size_lon / 2)) * np.pi / 180
    lon_side_max = (lon + (cell_size_lon / 2)) * np.pi / 180
    
    # Calculate the surface area of the cell
    S = (R**2)*(lon_side_max - lon_side_min)*(np.sin(lat_side_max) - np.sin(lat_side_min))
    
    return S

def grid05(data, lon_data, lat_data, lat_res, lon_res):
    """A function to grid lon, lat data.
    
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

    return {'Grid':xy_grid, 'Count':xy_count, 'Lon':x_range, 'Lat':y_range}

def apply_offset(month, boundary):
    if month == 'constant':
        Imnth = 13 - 1
    else:
        Imnth = int(month) - 1

    if boundary == 'SAR_ocean':
        f = open('/Users/jmh2g09/Documents/PhD/Data/SeparateModes/LRM-SAR_offset.dat', 'r')
        LRM_SAR_offset = []
        for line in f:
            line.strip()
            column = line.split()
            LRM_SAR_offset = np.append(LRM_SAR_offset, float(column[1]))
        f.close()
        
        return LRM_SAR_offset[Imnth]
        
    if boundary == 'ice':
        f = open('/Users/jmh2g09/Documents/PhD/Data/SeparateModes/ocean-ice_offset.dat', 'r')
        ocean_ice_offset = []
        for line in f:
            line.strip()
            column = line.split()
            ocean_ice_offset = np.append(ocean_ice_offset, float(column[1]))
        f.close()
        
        return ocean_ice_offset[Imnth]

def apply_offset_original(month, boundary):
    if boundary == 'SAR_ocean':
        if month == '01':
            return 0.005
        elif month == '02':
            return 0.007
        elif month == '03':
            return -0.009
        elif month == '04':
            return -0.018
        elif month == '05':
            return -0.010
        elif month == '06':
            return -0.019
        elif month == '07':
            return -0.030
        elif month == '08':
            return -0.022
        elif month == '09':
            return -0.016
        elif month == '10':
            return -0.017
        elif month == '11':
            return 0.001
        elif month == '12':
            return 0.012
        elif month == 'constant':
            return -0.010
    elif boundary == 'ice':
        if month == '01':
            return 0.052
        elif month == '02':
            return 0.044
        elif month == '03':
            return 0.034
        elif month == '04':
            return 0.043
        elif month == '05':
            return 0.047
        elif month == '06':
            return 0.065
        elif month == '07':
            return 0.075
        elif month == '08':
            return 0.075
        elif month == '09':
            return 0.068
        elif month == '10':
            return 0.072
        elif month == '11':
            return 0.050
        elif month == '12':
            return 0.050
        elif month == 'constant':
            return 0.056
#January LRM-SAR offset:  0.00489492491322
#Febuary LRM-SAR offset:  0.00688199790682
#March LRM-SAR offset:  -0.00897106525092
#April LRM-SAR offset:  -0.0181416179892
#May LRM-SAR offset:  -0.0101094678668
#June LRM-SAR offset:  -0.0187603584627
#July LRM-SAR offset:  -0.0295892239469
#August LRM-SAR offset:  -0.0219179266994
#September LRM-SAR offset:  -0.0162279361396
#October LRM-SAR offset:  -0.0174250239166
#November LRM-SAR offset:  0.00133077159719
#December LRM-SAR offset:  0.0124794197335

#January ocean-ice (SAR) offset:  0.0519413074729
#Febuary ocean-ice (SAR) offset:  0.0440994674621
#March ocean-ice (SAR) offset:  0.03386778963
#April ocean-ice (SAR) offset:  0.0432063759616
#May ocean-ice (SAR) offset:  0.0470435583904
#June ocean-ice (SAR) offset:  0.0645406477081
#July ocean-ice (SAR) offset:  0.0747555318892
#August ocean-ice (SAR) offset:  0.0745253922656
#September ocean-ice (SAR) offset:  0.0676124492288
#October ocean-ice (SAR) offset:  0.0723514436851
#November ocean-ice (SAR) offset:  0.0504626853533
#December ocean-ice (SAR) offset:  0.04950011136

def mode_points(lat, lon, month):
    '''Function that takes a list of lat and lon points and returns a list 
    corresponding to the SIRAL retracker mode with which that point was measured. 
    'month' is a string where '01' corresponds to January, '12' for December.
    
    Returns a list of length(lat) where the elements correspond to 
    Low Resolution Mode (0), SAR mode (1) or SARIn mode (2).'''
    
    def in_me(x, y, poly):
        '''Determine if a point is inside a given polygon or not
        Polygon is a list of (x,y) pairs.'''
    
        n = len(poly)
        inside = False

        p1x, p1y = poly[0]
        for i in range(n + 1):
            p2x,p2y = poly[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y

        return inside

    # Load the polygon data
    nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/SeparateModes/ModeMask/SARIn_polygon.nc', 'r')
    lat_poly_SARIn = nc.variables['Lat_SARIn'][:]
    lon_poly_SARIn = nc.variables['Lon_SARIn'][:]
    nc.close()
    nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/SeparateModes/ModeMask/SAR_polygon.nc', 'r')
    lat_poly_SAR = nc.variables['Lat_SAR_' + month][:]
    lon_poly_SAR = nc.variables['Lon_SAR_' + month][:]
    nc.close()

    # Convert to polar stereographic x and y coordinates
    m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
    # For the SAR data
    stereo_x_SAR, stereo_y_SAR = m(lon_poly_SAR, lat_poly_SAR)
    polygon_SAR = list(zip(stereo_x_SAR, stereo_y_SAR))
    # For the SARIn data
    stereo_x_SARIn, stereo_y_SARIn = m(lon_poly_SARIn, lat_poly_SARIn)
    polygon_SARIn = list(zip(stereo_x_SARIn, stereo_y_SARIn))
    # For the track data
    stereo_x, stereo_y = m(lon, lat)

    xy_pair = list(zip(stereo_x, stereo_y))
    point_type = []
    for point in xy_pair:
        point_x = point[0]
        point_y = point[1]
        in_SAR = in_me(point_x, point_y, polygon_SAR)
        in_SARIn = in_me(point_x, point_y, polygon_SARIn)
        
        if in_SAR == False and in_SARIn == False:
            point_type.append(0)
        elif in_SAR == True and in_SARIn == False:
            point_type.append(1)
        elif in_SARIn == True:
            point_type.append(2)
        else:
            print('Length of lat and mode: ', len(mode) == len(lat))
    
    return point_type

def inpaint_nans(y):
    def nan_helper(y):
        """Helper to handle indices and logical indices of NaNs.

        Input:
            - y, 1d numpy array with possible NaNs
        Output:
            - nans, logical indices of NaNs
            - index, a function, with signature indices= index(logical_indices),
              to convert logical indices of NaNs to 'equivalent' indices
        Example:
            >>> # linear interpolation of NaNs
            >>> nans, x= nan_helper(y)
            >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
        """

        return np.isnan(y), lambda z: z.nonzero()[0]
    nans, x= nan_helper(y)
    y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    
    return y

def gifmaker(giftitle, path):
    pngfiles = path + '*.png'
    giffile = path + giftitle + '.gif'
    os.system('magick -dispose previous ' + pngfiles + ' ' + giffile)
