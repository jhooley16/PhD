import os
import functions as funct
import numpy as np
import matplotlib.pyplot as pl

average_circumpolar_offset_LRM_SAR = np.zeros(12)
average_circumpolar_offset_ocean_ice = np.zeros(12)

WEDD_timeseries_LRM_SAR = []
IND_timeseries_LRM_SAR = []
ROSS_timeseries_LRM_SAR = []
AMBEL_timeseries_LRM_SAR = []
timeseries_LRM_SAR = []
WEDD_timeseries_ocean_ice = []
IND_timeseries_ocean_ice = []
ROSS_timeseries_ocean_ice = []
AMBEL_timeseries_ocean_ice = []
timeseries_ocean_ice = []


for year in ['2011', '2012', '2013', '2014', '2015']:
    print(year)
    
    monthly_offset_WEDD_LRM_SAR = []
    monthly_offset_IND_LRM_SAR = []
    monthly_offset_ROSS_LRM_SAR = []
    monthly_offset_AMBEL_LRM_SAR = []
    monthly_offset_LRM_SAR = []
    monthly_offset_WEDD_ocean_ice = []
    monthly_offset_IND_ocean_ice = []
    monthly_offset_ROSS_ocean_ice = []
    monthly_offset_AMBEL_ocean_ice = []
    monthly_offset_ocean_ice = []

    hist_offset_WEDD_LRM_SAR = []
    hist_offset_IND_LRM_SAR = []
    hist_offset_ROSS_LRM_SAR = []
    hist_offset_AMBEL_LRM_SAR = []
    hist_offset_circum_LRM_SAR = []
    hist_offset_WEDD_ocean_ice = []
    hist_offset_IND_ocean_ice = []
    hist_offset_ROSS_ocean_ice = []
    hist_offset_AMBEL_ocean_ice = []
    hist_offset_circum_ocean_ice = []

    month_number = []

    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        print(month)
        offset_WEDD_LRM_SAR = []
        offset_IND_LRM_SAR = []
        offset_ROSS_LRM_SAR = []
        offset_AMBEL_LRM_SAR = []
        circum_offset_LRM_SAR = []
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
                LRM_to_SAR = [0] * 200 + [1] * 200
                SAR_to_LRM = [1] * 200 + [0] * 200
                
                ocean_to_ice = [1] * 200 + [2] * 200
                ice_to_ocean = [2] * 200 + [1] * 200
                
                len_bound = len(LRM_to_SAR)             
                
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
                            ssha.append(float(columns[7]) - float(columns[8]))
                            surface.append(int(columns[0]))
                
                # Generate a list of retracker modes for this track
                tracker_type = funct.mode_points(lat, lon, month)

                # Find the boundaries
                iedge_LRM_SAR = []
                iedge_ocean_ice = []
                for it in range(len(tracker_type)):
                    # Find retracker boundaries
                    if tracker_type[it:it + len_bound] == LRM_to_SAR:
                        iedge_LRM_SAR.append(it + len_bound//2)
                    elif tracker_type[it:it + len_bound] == SAR_to_LRM:
                        iedge_LRM_SAR.append(it + len_bound//2)
                    # Find ice boundaries
                    if surface[it:it + len_bound] == ocean_to_ice:
                        iedge_ocean_ice.append(it + len_bound//2)
                    elif tracker_type[it:it + len_bound] == ice_to_ocean:
                        iedge_ocean_ice.append(it + len_bound//2)

                # For the LRM-SAR boundary
                for step in iedge_LRM_SAR:
                    # If the boundary is ALL OCEAN
                    if surface[step - len_bound//2:step + len_bound//2] == [1] * len_bound:
                        # If it's a LRM to SAR step
                        if tracker_type[step - len_bound//2:step + len_bound//2] == LRM_to_SAR:
                            if np.max(abs(np.gradient(lat[step - len_bound//2:step + len_bound//2]))) < 0.05:
                                if abs(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2])) < .3:
                                    # Calculate circumpolar offset
                                    circum_offset_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    hist_offset_circum_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    # Choose what the sector the offset point lies within
                                    if -60. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 0.:
                                        #print('Weddell')
                                        offset_WEDD_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                        hist_offset_WEDD_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    elif 0 <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 160.:
                                        #print('Indian')
                                        offset_IND_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                        hist_offset_IND_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    elif 160. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 180.:
                                        #print('Ross')
                                        offset_ROSS_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                        hist_offset_ROSS_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    elif -180. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= -130.:
                                        #print('Ross')
                                        offset_ROSS_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                        hist_offset_ROSS_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    elif -130. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= -60.:
                                        #print('Amundsen-Bellingshausen')
                                        offset_AMBEL_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                        hist_offset_AMBEL_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))

                        # If it's a SAR to LRM step
                        if tracker_type[step - len_bound//2:step + len_bound//2] == SAR_to_LRM:
                            if np.max(abs(np.gradient(lat[step - len_bound//2:step + len_bound//2]))) < 0.05:
                                if abs(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step])) < .3:
                                    # Calculate circumpolar offset
                                    circum_offset_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    hist_offset_circum_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    # Choose what the sector the offset point lies within
                                    if -60. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 0.:
                                        #print('Weddell')
                                        offset_WEDD_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                        hist_offset_WEDD_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    elif 0 <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 160.:
                                        #print('Indian')
                                        offset_IND_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                        hist_offset_IND_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    elif 160. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 180.:
                                        #print('Ross')
                                        offset_ROSS_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                        hist_offset_ROSS_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    elif -180. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= -130.:
                                        #print('Ross')
                                        offset_ROSS_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                        hist_offset_ROSS_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    elif -130. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= -60.:
                                        #print('Amundsen-Bellingshausen')
                                        offset_AMBEL_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                        hist_offset_AMBEL_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))

                for step in iedge_ocean_ice:
                    # If the boundary is ALL SAR
                    if tracker_type[step - len_bound//2:step + len_bound//2] == [1] * len_bound:
                        # If it's an ocean to ice step
                        if surface[step - len_bound//2:step + len_bound//2] == ocean_to_ice:
                            if np.max(abs(np.gradient(lat[step - len_bound//2:step + len_bound//2]))) < 0.05:
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
                            if np.max(abs(np.gradient(lat[step - len_bound//2:step + len_bound//2]))) < 0.05:
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

            monthly_offset_WEDD_LRM_SAR.append(np.mean(offset_WEDD_LRM_SAR))
            monthly_offset_IND_LRM_SAR.append(np.mean(offset_IND_LRM_SAR))
            monthly_offset_ROSS_LRM_SAR.append(np.mean(offset_ROSS_LRM_SAR))  
            monthly_offset_AMBEL_LRM_SAR.append(np.mean(offset_AMBEL_LRM_SAR))
            monthly_offset_LRM_SAR.append(np.mean(circum_offset_LRM_SAR))

            monthly_offset_WEDD_ocean_ice.append(np.mean(offset_WEDD_ocean_ice))
            monthly_offset_IND_ocean_ice.append(np.mean(offset_IND_ocean_ice))
            monthly_offset_ROSS_ocean_ice.append(np.mean(offset_ROSS_ocean_ice))  
            monthly_offset_AMBEL_ocean_ice.append(np.mean(offset_AMBEL_ocean_ice))
            monthly_offset_ocean_ice.append(np.mean(offset_circum_ocean_ice))

            WEDD_timeseries_LRM_SAR.append(np.mean(offset_WEDD_LRM_SAR))
            IND_timeseries_LRM_SAR.append(np.mean(offset_IND_LRM_SAR))
            ROSS_timeseries_LRM_SAR.append(np.mean(offset_ROSS_LRM_SAR)) 
            AMBEL_timeseries_LRM_SAR.append(np.mean(offset_AMBEL_LRM_SAR))
            timeseries_LRM_SAR.append(np.mean(circum_offset_LRM_SAR))
            
            WEDD_timeseries_ocean_ice.append(np.mean(offset_WEDD_ocean_ice))
            IND_timeseries_ocean_ice.append(np.mean(offset_IND_ocean_ice))
            ROSS_timeseries_ocean_ice.append(np.mean(offset_ROSS_ocean_ice))  
            AMBEL_timeseries_ocean_ice.append(np.mean(offset_AMBEL_ocean_ice))
            timeseries_ocean_ice.append(np.mean(offset_circum_ocean_ice))
    
    pl.figure()
    pl.hist([hist_offset_WEDD_LRM_SAR, hist_offset_IND_LRM_SAR, hist_offset_ROSS_LRM_SAR, hist_offset_AMBEL_LRM_SAR], label=['Weddell', 'Indian', 'Ross', 'Amundsen-Bellingshausen'])
    pl.legend()
    pl.title(year + ' LRM - SAR histogram')
    pl.xlabel('Offset Bin (m)')
    pl.ylabel('Frequency')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/' + year + '_LRM_SAR_offset_hist.png', format='png')
    pl.close()

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
    pl.plot(A, monthly_offset_LRM_SAR, label='Circumpolar', marker='.')
    pl.plot(A, monthly_offset_WEDD_LRM_SAR, label='Weddell', marker='.')
    pl.plot(A, monthly_offset_IND_LRM_SAR, label='Indian', marker='.')
    pl.plot(A, monthly_offset_ROSS_LRM_SAR, label='Ross', marker='.')
    pl.plot(A, monthly_offset_AMBEL_LRM_SAR, label='Amundsen-Bellingshausen', marker='.')
    pl.legend(loc='lower right')
    pl.title(year + ' LRM - SAR (m) offset')
    pl.ylabel('Offset (m)')
    pl.xlabel('Month')
    pl.xlim([1, 12])
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/' + year + '_LRM_SAR_offset_sectors.png', format='png')
    pl.close()
    
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
    
    average_circumpolar_offset_LRM_SAR += np.array(monthly_offset_LRM_SAR)
    average_circumpolar_offset_ocean_ice += np.array(monthly_offset_ocean_ice)

f=open('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/LRM_SAR_timeseries.txt', 'w')
print(WEDD_timeseries_LRM_SAR, IND_timeseries_LRM_SAR, ROSS_timeseries_LRM_SAR, AMBEL_timeseries_LRM_SAR, timeseries_LRM_SAR, file=f)
f.close()
f=open('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/ocean_ice_timeseries.txt', 'w')
print(WEDD_timeseries_ocean_ice, IND_timeseries_ocean_ice, ROSS_timeseries_ocean_ice, AMBEL_timeseries_ocean_ice, timeseries_ocean_ice, file=f)
f.close()

pl.figure()
pl.plot(timeseries_LRM_SAR, label='Circumpolar', marker='.')
pl.plot(WEDD_timeseries_LRM_SAR, label='Weddell', marker='.')
pl.plot(IND_timeseries_LRM_SAR, label='Indian', marker='.')
pl.plot(ROSS_timeseries_LRM_SAR, label='Ross', marker='.')
pl.plot(AMBEL_timeseries_LRM_SAR, label='Amundsen-Bellingshausen', marker='.')
pl.title('LRM - SAR (m) offset timeseries')
pl.xlabel('Month (from Jan 2011)')
pl.ylabel('LRM - SAR (m) offset')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/LRMSAR_timeseries.png', format='png')
pl.close()

pl.figure()
pl.plot(timeseries_ocean_ice, label='Circumpolar', marker='.')
pl.plot(WEDD_timeseries_ocean_ice, label='Weddell', marker='.')
pl.plot(IND_timeseries_ocean_ice, label='Indian', marker='.')
pl.plot(ROSS_timeseries_ocean_ice, label='Ross', marker='.')
pl.plot(AMBEL_timeseries_ocean_ice, label='Amundsen-Bellingshausen', marker='.')
pl.title('ocean - ice (m) (SAR) offset timeseries', marker='.')
pl.xlabel('Month (from Jan 2011)')
pl.ylabel('ocean - ice (m) (SAR) offset')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/SARSARIn_timeseries.png', format='png')
pl.close()

average_circumpolar_offset_LRM_SAR /= 5
average_circumpolar_offset_ocean_ice /= 5

print(average_circumpolar_offset_LRM_SAR)
print(average_circumpolar_offset_ocean_ice)

print('January LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[0])
print('Febuary LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[1])
print('March LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[2])
print('April LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[3])
print('May LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[4])
print('June LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[5])
print('July LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[6])
print('August LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[7])
print('September LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[8])
print('October LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[9])
print('November LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[10])
print('December LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[11])

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

#January LRM-SAR offset:  0.00605720510952
#Febuary LRM-SAR offset:  0.0077402495848
#March LRM-SAR offset:  -0.00528115908118
#April LRM-SAR offset:  -0.0122375772566
#May LRM-SAR offset:  -0.004982935448
#June LRM-SAR offset:  -0.0157234755177
#July LRM-SAR offset:  -0.0179909306721
#August LRM-SAR offset:  -0.0129872770328
#September LRM-SAR offset:  -0.0123001013372
#October LRM-SAR offset:  -0.0102008805335
#November LRM-SAR offset:  0.00215810350424
#December LRM-SAR offset:  0.0109285059187
#January SAR-SARIn offset:  -0.0010534507906
#Febuary SAR-SARIn offset:  -0.00611335649884
#March SAR-SARIn offset:  -0.00755736513066
#April SAR-SARIn offset:  -0.00960064774106
#May SAR-SARIn offset:  -0.00895712247207
#June SAR-SARIn offset:  -0.00773984717228
#July SAR-SARIn offset:  -0.00492843462298
#August SAR-SARIn offset:  -0.00431128796423
#September SAR-SARIn offset:  -0.00282061906754
#October SAR-SARIn offset:  -0.00117488669705
#November SAR-SARIn offset:  -0.00406919339594
#December SAR-SARIn offset:  -0.0011911462509

