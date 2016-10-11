from xml.dom import minidom
import numpy as np
from netCDF4 import Dataset

def month_number(mnth):
    if mnth == '01':
        return '01'
    if mnth == '02':
        return '03'
    if mnth == '03':
        return '05'
    if mnth == '04':
        return '07'
    if mnth == '05':
        return '09'
    if mnth == '06':
        return '11'
    if mnth == '07':
        return '13'
    if mnth == '08':
        return '15'
    if mnth == '09':
        return '17'
    if mnth == '10':
        return '19'
    if mnth == '11':
        return '21'
    if mnth == '12':
        return '23'

xmlfile = minidom.parse('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Mode Mask/CS-2_mode_mask.xml')
mapregions = xmlfile.getElementsByTagName('MapRegion')

lat_SARIn = []
lon_SARIn = []

for region in mapregions:
    description = region.getElementsByTagName('Description')
    # SARIn in Antarctica
    if region.attributes['id'].firstChild.data == 'CYFSIN06-00':
        mappoints = region.getElementsByTagName('MapPoints')[0]
        #print(region.attributes['id'].firstChild.data + ': ' + description[0].firstChild.data)
        for pair in mappoints.firstChild.data.split(';'):
            lat_SARIn.append(float(pair.split(',')[0]))
            lon_SARIn.append(float(pair.split(',')[1]))

# Make the last point the same as the first, to complete the polygon
lat_SARIn.append(lat_SARIn[0])
lon_SARIn.append(lon_SARIn[0])

nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Mode Mask/SARIn_polygon.nc', 'w', format='NETCDF3_CLASSIC')
nc.description = 'Points describing the polygon of the SARIn boundary around Antarctica'
nc.createDimension('station', np.size(lat_SARIn))

lat = nc.createVariable('Lat_SARIn', float, ('station',))
lon = nc.createVariable('Lon_SARIn', float, ('station',))
    
lat.long_name = 'latitude'
lat.standard_name = 'latitude'
lat.units = 'degrees_north'
lon.long_name = 'longitude'
lon.standard_name = 'longitude'
lon.units = 'degrees_east_-180to180'

lat[:] = lat_SARIn
lon[:] = lon_SARIn
nc.close()

nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Mode Mask/SAR_polygon.nc', 'w', format='NETCDF3_CLASSIC')

for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
    lat_SAR = []
    lon_SAR = []
    for region in mapregions:
        description = region.getElementsByTagName('Description')
        if region.attributes['id'].firstChild.data == 'CYSSAR02-' + month_number(month):
            mappoints = region.getElementsByTagName('MapPoints')[0]
            print(region.attributes['id'].firstChild.data + ': ' + description[0].firstChild.data)
            for pair in mappoints.firstChild.data.split(';'):
                lat_SAR.append(float(pair.split(',')[0]))
                lon_SAR.append(float(pair.split(',')[1]))
    
    # Make the last point the same as the first, to complete the polygon
    lat_SAR.append(lat_SAR[0])
    lon_SAR.append(lon_SAR[0])
    
    nc.createDimension('station_SAR_' + month, np.size(lat_SAR))

    # SARIn
    lat = nc.createVariable('Lat_SAR_' + month, float, ('station_SAR_'+month,))
    lon = nc.createVariable('Lon_SAR_' + month, float, ('station_SAR_'+month,))
    
    lat.long_name = 'latitude'
    lat.standard_name = 'latitude'
    lat.units = 'degrees_north'
    lon.long_name = 'longitude'
    lon.standard_name = 'longitude'
    lon.units = 'degrees_east'

    lat[:] = lat_SAR
    lon[:] = lon_SAR
nc.close()