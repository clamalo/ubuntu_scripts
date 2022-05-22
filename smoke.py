import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as colors
import cartopy
import os
from metpy.plots import USCOUNTIES
import multiprocessing
from datetime import timedelta
import datetime
import time

def datestr_and_cycle():
    import datetime
    now = datetime.datetime.now()
    datestr = now.strftime("%Y%m%d")
    hour = now.hour

    if 7 < hour < 19:
        cycle = '12'

    elif 19<=hour or 7>hour:
        cycle = '00'
        #if it's early but not in time for 0z, subtract 1 day and create new datestr
        if 14<=hour:
            datestr = datetime.datetime.strptime(datestr, '%Y%m%d')
            days = timedelta(1)
            datestr = str(datestr+days)
            datestr = datestr.split('-')
            year = str(datestr[0])
            month = str(datestr[1])
            day = str((datestr[2].split(' '))[0])
            datestr = str(datestr[2].split(' ')[1])
            datestr = str(year+month+day)
        else:
            datestr = datestr

    return datestr,cycle

#create the colormap
def create_colormap():
    import matplotlib.colors
    #six color listed colormap: green, yellow, orange, red, purple, maroon
    cmap = colors.ListedColormap(['#00FF00','#FFFF00','#FFA500','#FF0000','#800080','#800000'])
    return cmap

#FUNCTION: crop any given dataset to bounding dimensions set by domain_select
def crop(ds):
    # c = [49.41097, -125.70557,46.07323, -118.91602]
    # c = [40.84706, -102.08496,23.68477, -74.48730]
    # c = [55,-125,30,-60]
    c = [50,-125,31,-90]
    # c = [36.95688,-109.03387,31.93721,-100.89521]
    # c = [41.93032,-124.75193,33,-113.99778]
    # c = [41,-109,37,-104]
    max_mins = [(c[0]),(c[2]),(c[3]+360),(c[1]+360)]
    max_lat = float(max_mins[0])
    min_lat = float(max_mins[1])
    max_lon = float(max_mins[2])
    min_lon = float(max_mins[3])

    mask_lon = (ds.longitude >= min_lon) & (ds.longitude <= max_lon)
    mask_lat = (ds.latitude >= min_lat) & (ds.latitude <= max_lat)
    ds = ds.where(mask_lat, drop=True)
    ds = ds.where(mask_lon, drop=True)
    return ds

#function to add 0 to number if less than 10
def add_zero(num):
    if num < 10:
        return '0' + str(num)
    else:
        return str(num)

#read grib idx file
def read_idx(datestr,cycle,frame,idx_file):
    with open(idx_file, 'r') as f:
        lines = f.readlines()

    n = 76

    line1 = lines[n-1]
    line2 = lines[n]
    line1 = line1.split(':')
    line2 = line2.split(':')
    start_bytes = int(line1[1])
    end_bytes = int(line2[1])

    url = 'https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.'+datestr+'/conus/hrrr.t'+cycle+'z.wrfsfcf'+frame+'.grib2'
    file_name = '/root/hrrr'+str(int(frame))+'.grib2'

    os.system("curl "+url+" -r "+str(start_bytes)+"-"+str(end_bytes)+" > "+file_name)

# for n in range(0,48):
def ingest_and_plot(i):
    n = i
    datestr = datestr_and_cycle()[0]
    cycle = datestr_and_cycle()[1]
    #os system wget rtma file from noaa
    # os.system('wget -O /root/hrrr.grib2 "https://nomads.ncep.noaa.gov/cgi-bin/filter_hrrr_2d.pl?file=hrrr.t00z.wrfsfcf'+add_zero(n)+'.grib2&lev_8_m_above_ground=on&var_MASSDEN=on&leftlon=0&rightlon=360&toplat=90&bottomlat=-90&dir=%2Fhrrr.20220515%2Fconus"')

    idx_url = 'https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.'+datestr+'/conus/hrrr.t'+cycle+'z.wrfsfcf'+add_zero(n)+'.grib2.idx'
    os.system("curl "+idx_url+" > /root/hrrr"+str(n)+".idx")
    idx_file = '/root/hrrr'+str(n)+'.idx'

    read_idx(datestr,cycle,add_zero(n),idx_file)

    os.remove(idx_file)

    os.system('gdalwarp -t_srs EPSG:4326 hrrr'+str(n)+'.grib2 hrrr_test'+str(n)+'.grib2')

    ds = xr.load_dataset('/root/hrrr_test'+str(n)+'.grib2',engine='cfgrib')
    #remove rtma test
    # os.system('rm /root/hrrr_test'+str(n)+'.grib2')
    os.remove('/root/hrrr_test'+str(n)+'.grib2')
    os.remove('/root/hrrr'+str(n)+'.grib2')
    os.remove('/root/hrrr_test'+str(n)+'.grib2.923a8.idx')

    ds = crop(ds)
    new_lon = np.linspace(ds.longitude[0], ds.longitude[-1], ds.dims["longitude"] * 3)
    new_lat = np.linspace(ds.latitude[0], ds.latitude[-1], ds.dims["latitude"] * 3)
    ds = ds.interp(latitude=new_lat, longitude=new_lon)
    #0-12 = green, 12-35.4 = yellow, 35.4-55.4 = orange, 55.4-150.4 = red, 150.4-250.4 = purple, 250.4-500+ = maroon

    precip = ds.variables['unknown']*1000000000
    lats = ds.variables['latitude'][:]
    lons = ds.variables['longitude'][:]

    #initialize plot using cartopy plate caree
    fig = plt.figure(figsize=(20, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    #plot precip
    newcmp = create_colormap()
    bounds = [0,12,35.4,55.4,150.4,250.4,500.4]
    #norm = colors.BoundaryNorm(boundaries=bounds, ncolors=41)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=6)
    #norm = colors.BoundaryNorm(bounds, newcmp.N, clip=True)
    cf = ax.pcolormesh(lons, lats, precip, norm=norm, cmap=newcmp)
    # cf = ax.pcolormesh(lons,lats,precip,cmap='gnuplot')
    #add cartopy states
    ax.add_feature(cartopy.feature.STATES, linewidth=0.5, edgecolor='black')
    #add counties using metpy
    # ax.add_feature(USCOUNTIES.with_scale('500k'),linewidth=1,edgecolor='gray')
    #add colorbar
    cbar = plt.colorbar(cf, orientation='horizontal', pad=0.05, shrink=0.4)
    # cbar.set_ticks([12,35.4,55.4,150.4,250.4,500.4])
    cbar.set_ticks([6,23.7,45.4,102.9,200.4,374.4])
    cbar.set_ticklabels(['Good','Moderate','Unhealthy For \n Sensitive Groups','Unhealthy','Very Unhealthy','Hazardous'])
    cbar.ax.tick_params(labelsize=8)

    #create initialization time and valid time for given frame for plot title
    datestr = datestr_and_cycle()[0]
    cycle = datestr_and_cycle()[1]
    init_label = datestr[0:4]+'-'+datestr[4:6]+'-'+datestr[6:8]+' '+cycle+'z'
    datestr = str(datestr)+str(cycle)
    datestr = datetime.datetime.strptime(datestr, '%Y%m%d%H')
    hours_added = timedelta(hours = int(n))
    datestr = str(datestr+hours_added)
    valid_label = datestr[0:4]+'-'+datestr[5:7]+'-'+datestr[8:13]+'z'

    #set plot title
    plt.title("Gridded AQI Thresholds || Forecast Hour "+str(n)+" || Init "+init_label+" || Valid "+valid_label,fontsize=10)

    plt.savefig('/root/script/smoke/'+str(n)+'hrrr_plot.png',dpi=500,bbox_inches='tight')

if __name__ == '__main__':
    datestr_and_cycle()

    p = multiprocessing.Pool(5)
    p.map(ingest_and_plot, range(0,49))

    os.chdir('/root/script/smoke/')
    os.system('git add *')
    os.system('git commit -m "auto-push"')
    os.system('git checkout master')
    os.system('git pull')
    os.system('git config --global core.askpass "git-gui--askpass"')
    os.system('git push')
