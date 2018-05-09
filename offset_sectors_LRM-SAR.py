import os
import functions as funct
import numpy as np
import matplotlib.pyplot as pl
from datetime import date

average_circumpolar_offset_LRM_SAR = np.zeros(12)
average_errors = np.zeros(12)

WEDD_timeseries_LRM_SAR = []
IND_timeseries_LRM_SAR = []
ROSS_timeseries_LRM_SAR = []
AMBEL_timeseries_LRM_SAR = []
timeseries_LRM_SAR = []

dates = []

for year in ['2011']:#, '2012', '2013', '2014', '2015', '2016']:
    print(year)
    
    monthly_offset_WEDD_LRM_SAR = []
    monthly_offset_IND_LRM_SAR = []
    monthly_offset_ROSS_LRM_SAR = []
    monthly_offset_AMBEL_LRM_SAR = []
    monthly_offset_LRM_SAR = []
    
    errors = []

    hist_offset_WEDD_LRM_SAR = []
    hist_offset_IND_LRM_SAR = []
    hist_offset_ROSS_LRM_SAR = []
    hist_offset_AMBEL_LRM_SAR = []
    hist_offset_circum_LRM_SAR = []

    month_number = []

    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        offset_WEDD_LRM_SAR = []
        offset_IND_LRM_SAR = []
        offset_ROSS_LRM_SAR = []
        offset_AMBEL_LRM_SAR = []
        circum_offset_LRM_SAR = []
        
        dates.append(date(int(year), int(month), 15))
        
        if os.path.isdir('/Volumes/My Passport/Data/elev_files/' + year + month + '_MERGE'):
            os.chdir('/Volumes/My Passport/Data/elev_files/' + year + month + '_MERGE')
            month_number.append(int(month))
            print(month)
            for file in os.listdir():
                # 200 points equals the decorrelation scale (50 km)
                SC = 100
                LRM_to_SAR = [0] * SC + [1] * SC
                SAR_to_LRM = [1] * SC + [0] * SC
                
                len_bound = len(LRM_to_SAR)             
                
                ssha_pre = []
                lon_pre = []
                lat_pre = []
                surface_pre = []

                f = open(file, 'r')
                for line in f:
                    line = line.strip()
                    columns = line.split()
                    # If data point is listed as 'valid'
                    if columns[1] == '1':
                        # If data point is from open ocean (1) or from a lead (2)
                        if columns[0] == '1' or columns[0] == '2':
                            # If the ssh point is less than 3 m from the mssh
                            if abs(float(columns[7]) - float(columns[8])) < .3:
                                    lat_pre.append(float(columns[5]))
                                    lon_pre.append(float(columns[6]))
                                    surface_pre.append(int(columns[0]))
                                    ssha_pre.append(float(columns[7]) - float(columns[8]))
                f.close()
                
                descending = np.where(np.gradient(lat_pre) < 0.)[0]
                if len(descending) > 10:
                    inflection = descending[-1]
                    #### Descending ####
                
                    ssha_desc = np.array(ssha_pre[:inflection])
                    lat_desc = np.array(lat_pre[:inflection])
                    lon_desc = np.array(lon_pre[:inflection])
                    surface_desc = np.array(surface_pre[:inflection])
                
                    ssha_desc = ssha_desc[np.argsort(-lat_desc)]
                    lon_desc = lon_desc[np.argsort(-lat_desc)]
                    surface_desc = surface_desc[np.argsort(-lat_desc)]
                    lat_desc = lat_desc[np.argsort(-lat_desc)]
                
                    ## Filter the track
                    input_ssh = open('../INPUT_ssh.dat', 'w')
                    for ilat in range(len(lat_desc)):
                        print(-lat_desc[ilat], ssha_desc[ilat], file=input_ssh)
                    input_ssh.close()

                    os.system('gmt filter1d ../INPUT_ssh.dat -Fg0.2 -D0.001 -fi0y -E > ../OUTPUT_ssh.dat')
                    os.system('rm ../INPUT_ssh.dat')
                
                    output_ssh = open('../OUTPUT_ssh.dat', 'r')
                    lat_desc_filt = []
                    ssha_desc_filt = []
                    for line in output_ssh:
                        line.strip()
                        columns = line.split()
                        lat_desc_filt.append(-float(columns[0]))
                        ssha_desc_filt.append(float(columns[1]))
                    output_ssh.close()
                
                    os.system('rm ../OUTPUT_ssh.dat')
                
                    #### Ascending ######
                
                    ssha_asc = np.array(ssha_pre[inflection:])
                    lat_asc = np.array(lat_pre[inflection:])
                    lon_asc = np.array(lon_pre[inflection:])
                    surface_asc = np.array(surface_pre[inflection:])
                
                    ssha_asc = ssha_asc[np.argsort(lat_asc)]
                    lon_asc = lon_asc[np.argsort(lat_asc)]
                    surface_asc = surface_asc[np.argsort(lat_asc)]
                    lat_asc = lat_asc[np.argsort(lat_asc)]
                
                    ## Filter the track
                    input_ssh = open('../INPUT_ssh.dat', 'w')
                    for ilat in range(len(lat_asc)):
                        print(lat_asc[ilat], ssha_asc[ilat], file=input_ssh)
                    input_ssh.close()

                    os.system('gmt filter1d ../INPUT_ssh.dat -Fg0.2 -D0.001 -fi0y -E > ../OUTPUT_ssh.dat')
                    os.system('rm ../INPUT_ssh.dat')
                
                    output_ssh = open('../OUTPUT_ssh.dat', 'r')
                    lat_asc_filt = []
                    ssha_asc_filt = []
                    for line in output_ssh:
                        line.strip()
                        columns = line.split()
                        lat_asc_filt.append(float(columns[0]))
                        ssha_asc_filt.append(float(columns[1]))
                    output_ssh.close()
                
                    os.system('rm ../OUTPUT_ssh.dat')

                    lat =  list(lat_desc_filt) + list(lat_asc_filt)
                    ssha = list(ssha_desc_filt) + list(ssha_asc_filt)
                    lon =  list(lon_desc) + list(lon_asc)
                    surface = list(surface_desc) + list(surface_asc)
                

                    if len(lat) == len(lon):
                        # Generate a list of retracker modes for this track
                        tracker_type = funct.mode_points(lat, lon, month)

                        # Find the boundaries
                        iedge_LRM_SAR = []
                        for it in range(len(tracker_type)):
                            # Find retracker boundaries
                            if tracker_type[it:it + len_bound] == LRM_to_SAR:
                                iedge_LRM_SAR.append(it + len_bound//2)
                            elif tracker_type[it:it + len_bound] == SAR_to_LRM:
                                iedge_LRM_SAR.append(it + len_bound//2)

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

            monthly_offset_WEDD_LRM_SAR.append(np.mean(offset_WEDD_LRM_SAR))
            monthly_offset_IND_LRM_SAR.append(np.mean(offset_IND_LRM_SAR))
            monthly_offset_ROSS_LRM_SAR.append(np.mean(offset_ROSS_LRM_SAR))  
            monthly_offset_AMBEL_LRM_SAR.append(np.mean(offset_AMBEL_LRM_SAR))
            monthly_offset_LRM_SAR.append(np.mean(circum_offset_LRM_SAR))
            
            errors.append(np.nanstd(circum_offset_LRM_SAR) / np.sqrt(len(circum_offset_LRM_SAR)))

            WEDD_timeseries_LRM_SAR.append(np.mean(offset_WEDD_LRM_SAR))
            IND_timeseries_LRM_SAR.append(np.mean(offset_IND_LRM_SAR))
            ROSS_timeseries_LRM_SAR.append(np.mean(offset_ROSS_LRM_SAR)) 
            AMBEL_timeseries_LRM_SAR.append(np.mean(offset_AMBEL_LRM_SAR))
            timeseries_LRM_SAR.append(np.mean(circum_offset_LRM_SAR))
    
    pl.figure()
    pl.hist([hist_offset_WEDD_LRM_SAR, hist_offset_IND_LRM_SAR, hist_offset_ROSS_LRM_SAR, hist_offset_AMBEL_LRM_SAR], label=['Weddell', 'Indian', 'Ross', 'Amundsen-Bellingshausen'])
    pl.legend(loc='upper left')
    pl.title(year + ' LRM - SAR histogram')
    pl.xlabel('Offset Bin (m)')
    pl.ylabel('Frequency')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Offset/Figures/' + year + '_LRM_SAR_offset_hist.png', format='png', doi=300, transparent=True, bbox_inches='tight')
    pl.close()

    A = range(np.nanmin(month_number), np.nanmax(month_number) + 1)
    pl.figure()
    pl.plot(A, monthly_offset_LRM_SAR, label='Circumpolar', marker='.')
    pl.plot(A, monthly_offset_WEDD_LRM_SAR, label='Weddell', marker='.')
    pl.plot(A, monthly_offset_IND_LRM_SAR, label='Indian', marker='.')
    pl.plot(A, monthly_offset_ROSS_LRM_SAR, label='Ross', marker='.')
    pl.plot(A, monthly_offset_AMBEL_LRM_SAR, label='Amundsen-Bellingshausen', marker='.')
    pl.legend(loc='best')
    pl.title(year + ' LRM - SAR (m) offset')
    pl.ylabel('Offset (m)')
    pl.xlabel('Month')
    pl.xlim([1, 12])
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Offset/Figures/' + year + '_LRM_SAR_offset_sectors.png', format='png', doi=300, transparent=True, bbox_inches='tight')
    pl.close()
    
    average_circumpolar_offset_LRM_SAR += np.array(monthly_offset_LRM_SAR)
    average_errors += np.array(errors)

f=open('/Users/jmh2g09/Documents/PhD/Data/Offset/LRM_SAR_timeseries.txt', 'w')
for i in range(len(WEDD_timeseries_LRM_SAR)):
    print(WEDD_timeseries_LRM_SAR[i], IND_timeseries_LRM_SAR[i], ROSS_timeseries_LRM_SAR[i], AMBEL_timeseries_LRM_SAR[i], timeseries_LRM_SAR[i], file=f)
f.close()

fig = pl.figure()
pl.plot(dates, timeseries_LRM_SAR, label='Circumpolar', marker='.')
pl.plot(dates, WEDD_timeseries_LRM_SAR, label='Weddell', marker='.')
pl.plot(dates, IND_timeseries_LRM_SAR, label='Indian', marker='.')
pl.plot(dates, ROSS_timeseries_LRM_SAR, label='Ross', marker='.')
pl.plot(dates, AMBEL_timeseries_LRM_SAR, label='Amundsen-Bellingshausen', marker='.')
pl.legend(loc='lower right', prop={'size':6})
fig.autofmt_xdate()
pl.ylabel('LRM - SAR$_\mathrm{ocean}$ offset (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Offset/Figures/LRMSAR_timeseries.png', format='png', doi=300, transparent=True, bbox_inches='tight')
pl.close()

average_circumpolar_offset_LRM_SAR /= 6
average_errors /= 6

print('offset')
print(average_circumpolar_offset_LRM_SAR)
print('errors')
print(average_errors)

f = open('/Users/jmh2g09/Documents/PhD/Data/Offset/LRM-SAR_offset.dat', 'w')
for mnth in range(len(average_circumpolar_offset_LRM_SAR)):
    print(mnth + 1, average_circumpolar_offset_LRM_SAR[mnth], average_errors[mnth], file=f)

## Calculate the average offset for use as a constant (time) offset
constant_LRM_SAR_offset = np.nanmean(timeseries_LRM_SAR)
print(mnth + 2, constant_LRM_SAR_offset, file=f)
f.close()

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

print('Constant LRM-SAR offset: ', constant_LRM_SAR_offset)
