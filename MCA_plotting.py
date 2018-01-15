from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
from datetime import date

DOT_dates = []
for year in ['2011', '2012', '2013', '2014', '2015']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        DOT_dates.append(date(int(year), int(month), 15))

# Import file
MCA_file = '/Users/jmh2g09/Documents/MATLAB/MCA_ssh_stresscurl.nc'

nc = Dataset(MCA_file, 'r')

ssh_MCA1 = nc.variables['ssh_MCA1'][:]
ssh_MCA2 = nc.variables['ssh_MCA2'][:]
wind_MCA1 = nc.variables['stress_MCA1'][:]
wind_MCA2 = nc.variables['stress_MCA2'][:]

wind_MCA1[:, 0] = wind_MCA1[:, 1]
wind_MCA2[:, 0] = wind_MCA2[:, 1]

wind_MCA1[:, -1] = wind_MCA1[:, -2]
wind_MCA2[:, -1] = wind_MCA2[:, -2]

ssh_PC1 = nc.variables['ssh_PC1'][:]
ssh_PC2 = nc.variables['ssh_PC2'][:]
wind_PC1 = nc.variables['stress_PC1'][:]
wind_PC2 = nc.variables['stress_PC2'][:]

lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
nc.close()

grid_lats, grid_lons = np.meshgrid(lat, lon)

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(np.transpose(ssh_MCA1)), cmap='RdBu_r')
m.colorbar()
pl.clim(5, -5)
pl.title('MCA(1) SSH (cm)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/ssh_MCA1.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(np.transpose(wind_MCA1 * 1e7)), cmap='RdBu_r')
m.colorbar()
pl.clim(1, -1)
pl.title('MCA(1) ocean surface stress curl')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/stress_MCA1.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(np.transpose(ssh_MCA2)), cmap='RdBu_r')
m.colorbar()
pl.clim(5, -5)
pl.title('MCA(2) SSH (cm)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/ssh_MCA2.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(np.transpose(wind_MCA2 * 1e7)), cmap='RdBu_r')
m.colorbar()
pl.clim(1, -1)
pl.title('MCA(2) ocean surface stress curl')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/stress_MCA2.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

fig = pl.figure()
pl.plot(DOT_dates, ssh_PC1, label='PC(1) SSH', color='k', linestyle='-')
pl.plot(DOT_dates, wind_PC1, label='PC(1) ocean surface stress curl', color='k', linestyle='--')
fig.autofmt_xdate()
pl.legend(loc='upper right')
pl.title('PC(1)')
pl.ylim(-5, 5)
pl.ylabel('Normalised PC')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/PC1.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

fig = pl.figure()
pl.plot(DOT_dates, ssh_PC2, label='PC(2) SSH', color='k', linestyle='-')
pl.plot(DOT_dates, wind_PC2, label='PC(2) ocean surface stress curl', color='k', linestyle='--')
fig.autofmt_xdate()
pl.legend(loc='upper right')
pl.title('PC(2)')
pl.ylim(-5, 5)
pl.ylabel('Normalised PC')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/PC2.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()