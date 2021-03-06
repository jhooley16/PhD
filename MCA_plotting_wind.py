from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
from datetime import date
from scipy import stats

DOT_dates = []
for year in ['2011', '2012', '2013', '2014', '2015']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        DOT_dates.append(date(int(year), int(month), 15))

# Import file
MCA_file = '/Users/jmh2g09/Documents/MATLAB/MCA_ssh_stresscurl_wind.nc'

nc = Dataset(MCA_file, 'r')

ssh_MCA1 = nc.variables['ssh_MCA1'][:]
ssh_MCA2 = -nc.variables['ssh_MCA2'][:]
ssh_MCA3 = -nc.variables['ssh_MCA3'][:]
wind_MCA1 = nc.variables['stress_MCA1'][:]
wind_MCA2 = -nc.variables['stress_MCA2'][:]
wind_MCA3 = -nc.variables['stress_MCA3'][:]

wind_MCA1[:, 0] = wind_MCA1[:, 1]
wind_MCA2[:, 0] = wind_MCA2[:, 1]
wind_MCA3[:, 0] = wind_MCA3[:, 1]

wind_MCA1[:, -1] = wind_MCA1[:, -2]
wind_MCA2[:, -1] = wind_MCA2[:, -2]
wind_MCA3[:, -1] = wind_MCA3[:, -2]

ssh_PC1 = nc.variables['ssh_PC1'][:]
ssh_PC2 = -nc.variables['ssh_PC2'][:]
ssh_PC3 = -nc.variables['ssh_PC3'][:]

wind_PC1 = nc.variables['stress_PC1'][:]
wind_PC2 = -nc.variables['stress_PC2'][:]
wind_PC3 = -nc.variables['stress_PC3'][:]

EXPVAR = nc.variables['EXPVAR'][:]

lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
nc.close()

nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/Means/mean.nc', 'r')
ice_lat = nc.variables['lat'][:]
ice_lon = nc.variables['lon'][:]
sic = nc.variables['monthly_mean_sea_ice_concentration'][:]
nc.close()


##### Correlate PCs with SAM

# Load SAM data 
sam_file = '/Users/jmh2g09/Documents/PhD/Data/Wind/SAM_2011-2016.txt'
sam_index = []
f = open(sam_file, 'r')
for line in f:
    line = line.strip()
    sam_index.append(float(line))
f.close()
sam_index = np.array(sam_index[:60])

# Correlate with PCs
xcorr = stats.pearsonr(ssh_PC1, sam_index)
PC1_ssh_xcorr = xcorr[0]
PC1_ssh_xcorr_pvalue = xcorr[1]
xcorr = stats.pearsonr(ssh_PC2, sam_index)
PC2_ssh_xcorr = xcorr[0]
PC2_ssh_xcorr_pvalue = xcorr[1]
xcorr = stats.pearsonr(ssh_PC3, sam_index)
PC3_ssh_xcorr = xcorr[0]
PC3_ssh_xcorr_pvalue = xcorr[1]

xcorr = stats.pearsonr(wind_PC1, sam_index)
PC1_wind_xcorr = xcorr[0]
PC1_wind_xcorr_pvalue = xcorr[1]
xcorr = stats.pearsonr(wind_PC2, sam_index)
PC2_wind_xcorr = xcorr[0]
PC2_wind_xcorr_pvalue = xcorr[1]
xcorr = stats.pearsonr(wind_PC3, sam_index)
PC3_wind_xcorr = xcorr[0]
PC3_wind_xcorr_pvalue = xcorr[1]

# Correlate PCs with each other
#PC1
xcorr = stats.pearsonr(ssh_PC1, wind_PC1)
PC1_xcorr = xcorr[0]
PC1_xcorr_pvalue = xcorr[1]
#PC2
xcorr = stats.pearsonr(ssh_PC2, wind_PC2)
PC2_xcorr = xcorr[0]
PC2_xcorr_pvalue = xcorr[1]
#PC3
xcorr = stats.pearsonr(ssh_PC3, wind_PC3)
PC3_xcorr = xcorr[0]
PC3_xcorr_pvalue = xcorr[1]

grid_lats, grid_lons = np.meshgrid(lat, lon)
grid_lats_ice, grid_lons_ice = np.meshgrid(ice_lat, ice_lon)


pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
stereo_x, stereo_y = m(grid_lons, grid_lats)
stereo_x_ice, stereo_y_ice = m(grid_lons_ice, grid_lats_ice)
m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(np.transpose(ssh_MCA1)), cmap='RdBu_r')
m.colorbar()
pl.clim(5, -5)
m.contour(stereo_x_ice, stereo_y_ice, np.transpose(sic[:, :, 8]), [10, ], colors='k', linestyles='--')
m.contour(stereo_x_ice, stereo_y_ice, np.transpose(sic[:, :, 2]), [10, ], colors='k', linestyles='-')
pl.title('MCA(1) SSH (cm): ' + str(EXPVAR[0]) + '%')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/ssh_MCA1_wind.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
pl.clf()
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(np.transpose(wind_MCA1 * 1e7)), cmap='RdBu_r')
m.colorbar()
pl.clim(1, -1)
m.contour(stereo_x_ice, stereo_y_ice, np.transpose(sic[:, :, 8]), [10, ], colors='k', linestyles='--')
m.contour(stereo_x_ice, stereo_y_ice, np.transpose(sic[:, :, 2]), [10, ], colors='k', linestyles='-')
pl.title('MCA(1) wind stress curl')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/stress_MCA1_wind.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
pl.clf()
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(np.transpose(ssh_MCA2)), cmap='RdBu_r')
m.colorbar()
pl.clim(5, -5)
m.contour(stereo_x_ice, stereo_y_ice, np.transpose(sic[:, :, 8]), [10, ], colors='k', linestyles='--')
m.contour(stereo_x_ice, stereo_y_ice, np.transpose(sic[:, :, 2]), [10, ], colors='k', linestyles='-')
pl.title('MCA(2) SSH (cm): ' + str(EXPVAR[1]) + '%')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/ssh_MCA2_wind.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
pl.clf()
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(np.transpose(wind_MCA2 * 1e7)), cmap='RdBu_r')
m.colorbar()
pl.clim(1, -1)
m.contour(stereo_x_ice, stereo_y_ice, np.transpose(sic[:, :, 8]), [10, ], colors='k', linestyles='--')
m.contour(stereo_x_ice, stereo_y_ice, np.transpose(sic[:, :, 2]), [10, ], colors='k', linestyles='-')
pl.title('MCA(2) wind stress curl')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/stress_MCA2_wind.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()



pl.figure()
pl.clf()
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(np.transpose(ssh_MCA3)), cmap='RdBu_r')
m.colorbar()
pl.clim(5, -5)
m.contour(stereo_x_ice, stereo_y_ice, np.transpose(sic[:, :, 8]), [10, ], colors='k', linestyles='--')
m.contour(stereo_x_ice, stereo_y_ice, np.transpose(sic[:, :, 2]), [10, ], colors='k', linestyles='-')
pl.title('MCA(3) SSH (cm): ' + str(EXPVAR[2]) + '%')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/ssh_MCA3_wind.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
pl.clf()
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(np.transpose(wind_MCA3 * 1e7)), cmap='RdBu_r')
m.colorbar()
pl.clim(1, -1)
m.contour(stereo_x_ice, stereo_y_ice, np.transpose(sic[:, :, 8]), [10, ], colors='k', linestyles='--')
m.contour(stereo_x_ice, stereo_y_ice, np.transpose(sic[:, :, 2]), [10, ], colors='k', linestyles='-')
pl.title('MCA(3) wind stress curl')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/stress_MCA3_wind.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

fig = pl.figure()
pl.plot(DOT_dates, ssh_PC1, label='PC(1) SSH (SAM: ' + str(round(PC1_ssh_xcorr, 2)) +  ')', color='k', linestyle='-')
pl.plot(DOT_dates, wind_PC1, label='PC(1) WSC (SAM: ' + str(round(PC1_wind_xcorr, 2)) + ')', color='k', linestyle='--')
fig.autofmt_xdate()
pl.legend(loc='upper right')
pl.title('PC(1) (Corr: ' + str(round(PC1_xcorr, 2)) + ')')
pl.ylim(-5, 5)
pl.ylabel('Normalised PC ')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/PC1_wind.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

fig = pl.figure()
pl.plot(DOT_dates, ssh_PC2, label='PC(2) SSH (SAM: ' + str(round(PC2_ssh_xcorr, 2)) + ')', color='k', linestyle='-')
pl.plot(DOT_dates, wind_PC2, label='PC(2) WSC (SAM: ' + str(round(PC2_wind_xcorr, 2)) + ')', color='k', linestyle='--')
fig.autofmt_xdate()
pl.legend(loc='upper right')
pl.title('PC(2) (Corr: ' + str(round(PC2_xcorr, 2)) + ')')
pl.ylim(-5, 5)
pl.ylabel('Normalised PC')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/PC2_wind.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

fig = pl.figure()
pl.plot(DOT_dates, ssh_PC3, label='PC(3) SSH (SAM: ' + str(round(PC3_ssh_xcorr, 2)) + ')', color='k', linestyle='-')
pl.plot(DOT_dates, wind_PC3, label='PC(3) WSC (SAM: ' + str(round(PC3_wind_xcorr, 2)) + ')', color='k', linestyle='--')
fig.autofmt_xdate()
pl.legend(loc='upper right')
pl.title('PC(3) (Corr: ' + str(round(PC3_xcorr, 2)) + ')')
pl.ylim(-5, 5)
pl.ylabel('Normalised PC')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/MCA/PC3_wind.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()