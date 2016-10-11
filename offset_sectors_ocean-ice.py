import os
import functions as funct
import numpy as np
import matplotlib.pyplot as pl

average_circumpolar_offset_ocean_ice = np.zeros(12)

WEDD_timeseries_ocean_ice = []
IND_timeseries_ocean_ice = []
ROSS_timeseries_ocean_ice = []
AMBEL_timeseries_ocean_ice = []
timeseries_ocean_ice = []


for year in ['2011', '2012', '2013', '2014', '2015']:
    print(year)
    
    monthly_offset_WEDD_ocean_ice = []
    monthly_offset_IND_ocean_ice = []
    monthly_offset_ROSS_ocean_ice = []
    monthly_offset_AMBEL_ocean_ice = []
    monthly_offset_ocean_ice = []

    hist_offset_WEDD_ocean_ice = []
    hist_offset_IND_ocean_ice = []
    hist_offset_ROSS_ocean_ice = []
    hist_offset_AMBEL_ocean_ice = []
    hist_offset_circum_ocean_ice = []

    month_number = []

    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        print(month)
        offset_WEDD_ocean_ice = []
        offset_IND_ocean_ice = []
        offset_ROSS_ocean_ice = []
        offset_AMBEL_ocean_ice = []
        offset_circum_ocean_ice = []

        if os.path.isdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + year + month + '_elev'):
            os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + year + month + '_elev')
            month_number.append(int(month))
        
            for file in os.listdir():
                
                # 200 points equals the decorrelation scale (50 km)
                SC = 100
                
                ocean_to_ice = [1] * SC + [2] * SC
                ice_to_ocean = [2] * SC + [1] * SC
                
                len_bound = len(ocean_to_ice)             
                
                ssha = []
                lon = []
                lat = []
                surface = []

                f = open(file, 'r')
                for line in f:
                    line = line.strip()
                    columns = line.split()
                    # If data point is from open ocean (1) or from a lead (2)
                    if columns[0] == '1' or columns[0] == '2':
                        # If data point is listed as 'valid'
                        if columns[1] == '1':
                            lat.append(float(columns[5]))
                            lon.append(float(columns[6]))
                            ssha.append((float(columns[7]) - float(columns[8])))
                            surface.append(int(columns[0]))

                tracker_type = funct.mode_points(lat, lon, month)
                
                for point in range(len(tracker_type)):
                    if tracker_type[point] == 1 and surface[point] == '1':
                        ssha[point] += funct.apply_offset(month, 'SAR')

                # Find the boundaries
                iedge_ocean_ice = []
                for it in range(len(surface)):
                    # Find ice boundaries
                    if surface[it:it + len_bound] == ocean_to_ice:
                        iedge_ocean_ice.append(it + len_bound//2)
                    elif surface[it:it + len_bound] == ice_to_ocean:
                        iedge_ocean_ice.append(it + len_bound//2)

                for step in iedge_ocean_ice:
                    # If it's an ocean to ice step
                    if surface[step - len_bound//2:step + len_bound//2] == ocean_to_ice:
                        if np.max(abs(np.gradient(lat[step - len_bound//2:step + len_bound//2]))) < 1:
                            if abs(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2])) < .3:
                                # Calculate circumpolar offset
                                offset_circum_ocean_ice.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                hist_offset_circum_ocean_ice.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                # Choose what the sector the offset point lies within
                                if -60. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 0.:
                                    #print('Weddell')
                                    offset_WEDD_ocean_ice.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    hist_offset_WEDD_ocean_ice.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                elif 0 <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 160.:
                                    #print('Indian')
                                    offset_IND_ocean_ice.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    hist_offset_IND_ocean_ice.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                elif 160. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 180.:
                                    #print('Ross')
                                    offset_ROSS_ocean_ice.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    hist_offset_ROSS_ocean_ice.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                elif -180. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= -130.:
                                    #print('Ross')
                                    offset_ROSS_ocean_ice.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    hist_offset_ROSS_ocean_ice.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                elif -130. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= -60.:
                                    #print('Amundsen-Bellingshausen')
                                    offset_AMBEL_ocean_ice.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    hist_offset_AMBEL_ocean_ice.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))

                    # If it's an ice to ocean step
                    if surface[step - len_bound//2:step + len_bound//2] == ice_to_ocean:
                        if np.max(abs(np.gradient(lat[step - len_bound//2:step + len_bound//2]))) < 1:
                            if abs(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step])) < .3:
                                # Calculate circumpolar offset
                                offset_circum_ocean_ice.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                hist_offset_circum_ocean_ice.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                # Choose what the sector the offset point lies within
                                if -60. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 0.:
                                    #print('Weddell')
                                    offset_WEDD_ocean_ice.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    hist_offset_WEDD_ocean_ice.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                elif 0 <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 160.:
                                    #print('Indian')
                                    offset_IND_ocean_ice.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    hist_offset_IND_ocean_ice.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                elif 160. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 180.:
                                    #print('Ross')
                                    offset_ROSS_ocean_ice.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    hist_offset_ROSS_ocean_ice.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                elif -180. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= -130.:
                                    #print('Ross')
                                    offset_ROSS_ocean_ice.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    hist_offset_ROSS_ocean_ice.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                elif -130. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= -60.:
                                    #print('Amundsen-Bellingshausen')
                                    offset_AMBEL_ocean_ice.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    hist_offset_AMBEL_ocean_ice.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))

            monthly_offset_WEDD_ocean_ice.append(np.mean(offset_WEDD_ocean_ice))
            monthly_offset_IND_ocean_ice.append(np.mean(offset_IND_ocean_ice))
            monthly_offset_ROSS_ocean_ice.append(np.mean(offset_ROSS_ocean_ice))  
            monthly_offset_AMBEL_ocean_ice.append(np.mean(offset_AMBEL_ocean_ice))
            monthly_offset_ocean_ice.append(np.mean(offset_circum_ocean_ice))

            WEDD_timeseries_ocean_ice.append(np.mean(offset_WEDD_ocean_ice))
            IND_timeseries_ocean_ice.append(np.mean(offset_IND_ocean_ice))
            ROSS_timeseries_ocean_ice.append(np.mean(offset_ROSS_ocean_ice))  
            AMBEL_timeseries_ocean_ice.append(np.mean(offset_AMBEL_ocean_ice))
            timeseries_ocean_ice.append(np.mean(offset_circum_ocean_ice))

    pl.figure()
    pl.hist([hist_offset_WEDD_ocean_ice, hist_offset_IND_ocean_ice, hist_offset_ROSS_ocean_ice, hist_offset_AMBEL_ocean_ice], label=['Weddell', 'Indian', 'Ross', 'Amundsen-Bellingshausen'])
    pl.legend()
    pl.title(year + ' ocean - ice (SAR) histogram')
    pl.xlabel('Offset Bin (m)')
    pl.ylabel('Frequency')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/' + year + '_ocean_ice_offset_hist.png', format='png')
    pl.close()

    A = range(np.min(month_number), np.max(month_number) + 1)
    
    pl.figure()
    pl.plot(A, monthly_offset_ocean_ice, label='Circumpolar', marker='.')
    pl.plot(A, monthly_offset_WEDD_ocean_ice, label='Weddell', marker='.')
    pl.plot(A, monthly_offset_IND_ocean_ice, label='Indian', marker='.')
    pl.plot(A, monthly_offset_ROSS_ocean_ice, label='Ross', marker='.')
    pl.plot(A, monthly_offset_AMBEL_ocean_ice, label='Amundsen-Bellingshausen', marker='.')
    pl.legend(loc='lower right')
    pl.title(year + ' ocean - ice (m) (SAR) offset')
    pl.ylabel('Offset (m)')
    pl.xlabel('Month')
    pl.xlim([1, 12])
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/' + year + '_ocean_ice_offset_sectors.png', format='png')
    pl.close()

    average_circumpolar_offset_ocean_ice += np.array(monthly_offset_ocean_ice)

f=open('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/ocean_ice_timeseries.txt', 'w')
for i in range(len(WEDD_timeseries_ocean_ice)):
    print(WEDD_timeseries_ocean_ice[i], IND_timeseries_ocean_ice[i], ROSS_timeseries_ocean_ice[i], AMBEL_timeseries_ocean_ice[i], timeseries_ocean_ice[i], file=f)
f.close()

pl.figure()
pl.plot(timeseries_ocean_ice, label='Circumpolar', marker='.')
pl.plot(WEDD_timeseries_ocean_ice, label='Weddell', marker='.')
pl.plot(IND_timeseries_ocean_ice, label='Indian', marker='.')
pl.plot(ROSS_timeseries_ocean_ice, label='Ross', marker='.')
pl.plot(AMBEL_timeseries_ocean_ice, label='Amundsen-Bellingshausen', marker='.')
pl.title('ocean - ice (m) (SAR) offset timeseries')
pl.xlabel('Month (from Jan 2011)')
pl.ylabel('ocean - ice (m) (SAR) offset')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/ocean_ice_timeseries.png', format='png')
pl.close()

average_circumpolar_offset_ocean_ice /= 5

print(average_circumpolar_offset_ocean_ice)

print('January ocean-ice (SAR) offset: ', average_circumpolar_offset_ocean_ice[0])
print('Febuary ocean-ice (SAR) offset: ', average_circumpolar_offset_ocean_ice[1])
print('March ocean-ice (SAR) offset: ', average_circumpolar_offset_ocean_ice[2])
print('April ocean-ice (SAR) offset: ', average_circumpolar_offset_ocean_ice[3])
print('May ocean-ice (SAR) offset: ', average_circumpolar_offset_ocean_ice[4])
print('June ocean-ice (SAR) offset: ', average_circumpolar_offset_ocean_ice[5])
print('July ocean-ice (SAR) offset: ', average_circumpolar_offset_ocean_ice[6])
print('August ocean-ice (SAR) offset: ', average_circumpolar_offset_ocean_ice[7])
print('September ocean-ice (SAR) offset: ', average_circumpolar_offset_ocean_ice[8])
print('October ocean-ice (SAR) offset: ', average_circumpolar_offset_ocean_ice[9])
print('November ocean-ice (SAR) offset: ', average_circumpolar_offset_ocean_ice[10])
print('December ocean-ice (SAR) offset: ', average_circumpolar_offset_ocean_ice[11])

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

