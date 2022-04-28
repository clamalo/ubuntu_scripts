import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import math
from datetime import datetime, timedelta
import os
import requests
import matplotlib.colors as colors
from metpy.plots import USCOUNTIES
from osgeo import gdal

#create the colormap
def create_colormap():
    import matplotlib.colors
    cmap = colors.ListedColormap(['#ffffff', '#bdbfbd','#aba5a5', '#838383', '#6e6e6e', '#b5fbab', '#95f58b','#78f572', '#50f150','#1eb51e', '#0ea10e', '#1464d3', '#2883f1', '#50a5f5','#97d3fb', '#b5f1fb','#fffbab', '#ffe978', '#ffc13c', '#ffa100', '#ff6000','#ff3200', '#e11400','#c10000', '#a50000', '#870000', '#643c31', '#8d6558','#b58d83', '#c7a095','#f1ddd3', '#cecbdc'])#, '#aca0c7', '#9b89bd', '#725ca3','#695294', '#770077','#8d008d', '#b200b2', '#c400c4', '#db00db'])
    return cmap

#FUNCTION: make frame as a string in the correct format for NOMADS request
def name_frame(frame):
    if len(str(frame)) == 1:
        frame = '0'+str(frame)
    else:
        frame = str(frame)
    return frame

#FUNCTION: determine the date and time from which to pull data from
def datestr_and_cycle():
    #pull current year, month, day, and hour in UTC time
    datestr = str(datetime.utcnow())
    datestr = datestr.split('-')
    year = str(datestr[0])
    month = str(datestr[1])
    day = str((datestr[2].split(' '))[0])
    datestr = str(datestr[2].split(' ')[1])
    hour = int(datestr[0:2])
    #datestr to pull from NOMADS
    datestr = str(year+month+day)

    #logic checks to make sure it's not pulling data before it's done on NOMADS
    if 3<=hour<15:
        cycle = '00'
    # if 9<=hour<15:
    #     cycle = '06'
    elif 15<=hour or hour<3:
        cycle = '12'
    # elif 21<=hour or 3>hour:
    #     cycle = '18'
        #if it's early but not in time for 0z, subtract 1 day and create new datestr
        if 3>hour:
            datestr = datetime.strptime(datestr, '%Y%m%d')
            days = timedelta(1)
            datestr = str(datestr-days)
            datestr = datestr.split('-')
            year = str(datestr[0])
            month = str(datestr[1])
            day = str((datestr[2].split(' '))[0])
            datestr = str(datestr[2].split(' ')[1])
            datestr = str(year+month+day)
        else:
            datestr = datestr
    return datestr,cycle,hour

#read grib idx file
def read_idx(idx_file,model,frame,cycle,datestr):
    # datestr = datestr_and_cycle()[0]
    # cycle = datestr_and_cycle()[1]

    with open(idx_file, 'r') as f:
        lines = f.readlines()

    if model == 'arw':
        if frame==5:
            n = 16
        elif frame==1 or frame==2 or frame==4 or frame==5 or frame==7 or frame==8 or frame==10 or frame==11 or frame==13 or frame==14 or frame==16 or frame==17 or frame==19 or frame==20 or frame==22 or frame==23 or frame==25 or frame==26 or frame==28 or frame==29 or frame==31 or frame==32 or frame==34 or frame==35:
            n = 20
        else:
            n = 15

    if model == 'fv3':
        if frame%3 != 0:
            n = 17
        else:
            n = 12

    if model == 'nam':
        n = -2
    
    if model == 'hrrr':
        n = 90

    if model == 'arw5k_1' or model == 'arw5k_2':
        n = 35
    
    if model == 'fv35k':
        n = 34
    
    if model == 'arw2p5k':
        if frame%3 == 0:
            n = 15
        else:
            n = 20
    
    if model == 'fv32p5k':
        if frame%3 == 0:
            n = 12
        else:
            n = 17

    line1 = lines[n-1]
    line2 = lines[n]
    line1 = line1.split(':')
    line2 = line2.split(':')
    start_bytes = int(line1[1])
    end_bytes = int(line2[1])
    print(line1)
    print(line2)
    print(model,frame)

    if model == 'arw':
        url = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/hiresw/prod/hiresw.'+datestr+'/hiresw.t'+cycle+'z.arw_2p5km.f'+name_frame(frame)+'.conus.grib2'

    if model == 'nam':
        url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/nam/prod/nam.'+datestr+'/nam.t'+cycle+'z.conusnest.hiresf'+name_frame(frame)+'.tm00.grib2'

    if model == 'fv3':
        url = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/hiresw/prod/hiresw.'+datestr+'/hiresw.t'+cycle+'z.fv3_2p5km.f'+name_frame(frame)+'.conus.grib2'

    if model == 'hrrr':
        url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/hrrr/prod/hrrr.'+datestr+'/conus/hrrr.t'+cycle+'z.wrfsfcf'+name_frame(frame)+'.grib2'

    if model == 'arw5k_1':
        url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.'+datestr+'/hiresw.t'+cycle+'z.arw_5km.f'+name_frame(frame)+'.conus.grib2'

    if model == 'arw5k_2':
        url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.'+datestr+'/hiresw.t'+cycle+'z.arw_5km.f'+name_frame(frame)+'.conusmem2.grib2'

    if model == 'fv35k':
        url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.'+datestr+'/hiresw.t'+cycle+'z.fv3_5km.f'+name_frame(frame)+'.conus.grib2'
    
    if model == 'arw2p5k':
        url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.'+datestr+'/hiresw.t'+cycle+'z.arw_2p5km.f'+name_frame(frame)+'.conus.grib2'
    
    if model == 'fv32p5k':
        url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.'+datestr+'/hiresw.t'+cycle+'z.fv3_2p5km.f'+name_frame(frame)+'.conus.grib2'

    file_name = '/root/current.grib2'
    curl_message = ('curl '+url+' -r '+str(start_bytes)+'-'+str(end_bytes)+' > '+file_name)
    print(curl_message)
    os.system(curl_message)
    # ds = xr.load_dataset(file_name,engine='cfgrib')

def nam3k(chelsa_ds,frame,cycle,datestr,offset):

    #remove files
    file_exists = os.path.exists('/root/current_.grib2')
    if file_exists == True:
        os.remove('/root/current_.grib2')
    file_exists = os.path.exists('/root/minus_one_.grib2')
    if file_exists == True:
        os.remove('/root/minus_one_.grib2')

    frame = name_frame(int(frame)+offset)
    idx_url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/nam/prod/nam.'+datestr+'/nam.t'+cycle+'z.conusnest.hiresf'+frame+'.tm00.grib2.idx'
    os.system('curl "'+idx_url+'" --output "/root/nam.t'+cycle+'z.conusnest.hiresf'+frame+'.tm00.grib2.idx"')
    idx_file = '/root/nam.t'+cycle+'z.conusnest.hiresf'+frame+'.tm00.grib2.idx'
    read_idx(idx_file,'nam',int(frame),cycle,datestr)
    (xr.load_dataset('/root/current.grib2')).to_netcdf('/root/current.nc')
    os.system('/root/anaconda3/envs/blend/bin/gdalwarp -t_srs EPSG:4326 /root/current.nc /root/current_.tif')
    inputfile = '/root/current_.tif'
    outputfile = '/root/current_.nc'
    ds = gdal.Translate(outputfile, inputfile, format='NetCDF')
    os.remove(idx_file)
    dataset = xr.load_dataset('/root/current_.nc')
    if 'crs' in str(dataset):
        dataset = dataset.drop(['crs'])
    
    for n in range(len(dataset.lat)):
        values = dataset.tp[n].values
        output = []
        for value in zip(values):
            if str(value[0]) == 'nan':
                value = 0
            else:
                value = value[0]
            output.append(value)
            # print(value)
        dataset['tp'][n] = output

    if (int(frame)-1)%3 != 0:
        frame = name_frame((int(frame)-1))
        idx_url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/nam/prod/nam.'+datestr+'/nam.t'+cycle+'z.conusnest.hiresf'+frame+'.tm00.grib2.idx'
        os.system('curl "'+idx_url+'" --output "/root/nam.t'+cycle+'z.conusnest.hiresf'+frame+'.tm00.grib2.idx"')
        idx_file = '/root/nam.t'+cycle+'z.conusnest.hiresf'+frame+'.tm00.grib2.idx'
        read_idx(idx_file,'nam',int(frame),cycle,datestr)
        (xr.load_dataset('/root/current.grib2')).to_netcdf('/root/current.nc')
        os.system('/root/anaconda3/envs/blend/bin/gdalwarp -t_srs EPSG:4326 /root/current.nc /root/minus_one_.tif')
        inputfile = '/root/minus_one_.tif'
        outputfile = '/root/minus_one_.nc'
        ds = gdal.Translate(outputfile, inputfile, format='NetCDF')
        prior_dataset = xr.load_dataset('/root/minus_one_.nc')
        if 'crs' in str(prior_dataset):
            prior_dataset = prior_dataset.drop(['crs'])
        os.remove(idx_file)
        for n in range(len(prior_dataset.lat)):
            values = prior_dataset.tp[n].values
            output = []
            for value in zip(values):
                if str(value[0]) == 'nan':
                    value = 0
                else:
                    value = value[0]
                # print(value)
                output.append(value)
            prior_dataset['tp'][n] = output
        for n in range(len(dataset.lat)):
            current_list = dataset.tp[n].values
            prior_list = prior_dataset.tp[n].values
            output = []
            for current,prior in zip(current_list,prior_list):
                current = current-prior
                output.append(current)
            dataset['tp'][n] = output
        # dataset['tp'] = dataset['tp']-prior_dataset['tp']
        # for n in range(len(dataset.lat)):
        #     print(max(dataset.tp[n].values))

    dataset['lon'] = dataset['lon']+360
    dataset = crop_ds(dataset,'180_chelsa')
    for n in range(len(dataset.lat)):
        print(max(dataset.tp[n].values))
    dataset = dataset.interp(lat=chelsa_ds["lat"], lon=chelsa_ds["lon"])
    dataset['tp'] = dataset['tp']*chelsa_ds['precip']

    return dataset

def hrrr3k(chelsa_ds,frame,cycle,datestr,offset):

    #remove files
    file_exists = os.path.exists('/root/current_.grib2')
    if file_exists == True:
        os.remove('/root/current_.grib2')
    file_exists = os.path.exists('/root/minus_one_.grib2')
    if file_exists == True:
        os.remove('/root/minus_one_.grib2')

    frame = name_frame(int(frame)+offset)
    idx_url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/hrrr/prod/hrrr.'+datestr+'/conus/hrrr.t'+cycle+'z.wrfsfcf'+frame+'.grib2.idx'
    os.system('curl "'+idx_url+'" --output "/root/hrrr.t'+cycle+'z.wrfsfcf'+frame+'.grib2.idx"')
    idx_file = '/root/hrrr.t'+cycle+'z.wrfsfcf'+frame+'.grib2.idx'
    read_idx(idx_file,'hrrr',int(frame),cycle,datestr)
    (xr.load_dataset('/root/current.grib2')).to_netcdf('/root/current.nc')
    os.system('/root/anaconda3/envs/blend/bin/gdalwarp -t_srs EPSG:4326 /root/current.nc /root/current_.tif')
    inputfile = '/root/current_.tif'
    outputfile = '/root/current_.nc'
    ds = gdal.Translate(outputfile, inputfile, format='NetCDF')
    os.remove(idx_file)
    dataset = xr.load_dataset('/root/current_.nc')
    if 'crs' in str(dataset):
        dataset = dataset.drop(['crs'])

    dataset['lon'] = dataset['lon']+360
    dataset = crop_ds(dataset,'180_chelsa')
    dataset = dataset.interp(lat=chelsa_ds["lat"], lon=chelsa_ds["lon"])
    dataset['tp'] = dataset['tp']*chelsa_ds['precip']

    return dataset

def arw5k_1(chelsa_ds,frame,cycle,datestr,offset):

    #remove files
    file_exists = os.path.exists('/root/current_.grib2')
    if file_exists == True:
        os.remove('/root/current_.grib2')
    file_exists = os.path.exists('/root/minus_one_.grib2')
    if file_exists == True:
        os.remove('/root/minus_one_.grib2')

    frame = name_frame(int(frame)+offset)
    idx_url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.'+datestr+'/hiresw.t'+cycle+'z.arw_5km.f'+frame+'.conus.grib2.idx'
    os.system('curl "'+idx_url+'" --output "/root/hiresw.t'+cycle+'z.arw_5km.f'+frame+'.conus.grib2.idx"')
    idx_file = '/root/hiresw.t'+cycle+'z.arw_5km.f'+frame+'.conus.grib2.idx'
    read_idx(idx_file,'arw5k_1',int(frame),cycle,datestr)
    (xr.load_dataset('/root/current.grib2')).to_netcdf('/root/current.nc')
    os.system('/root/anaconda3/envs/blend/bin/gdalwarp -t_srs EPSG:4326 /root/current.nc /root/current_.tif')
    inputfile = '/root/current_.tif'
    outputfile = '/root/current_.nc'
    ds = gdal.Translate(outputfile, inputfile, format='NetCDF')
    os.remove(idx_file)
    dataset = xr.load_dataset('/root/current_.nc')
    if 'crs' in str(dataset):
        dataset = dataset.drop(['crs'])

    dataset['lon'] = dataset['lon']+360
    dataset = crop_ds(dataset,'180_chelsa')
    dataset = dataset.interp(lat=chelsa_ds["lat"], lon=chelsa_ds["lon"])
    dataset['tp'] = dataset['tp']*chelsa_ds['precip']

    return dataset

def arw5k_2(chelsa_ds,frame,cycle,datestr,offset):

    #remove files
    file_exists = os.path.exists('/root/current_.grib2')
    if file_exists == True:
        os.remove('/root/current_.grib2')
    file_exists = os.path.exists('/root/minus_one_.grib2')
    if file_exists == True:
        os.remove('/root/minus_one_.grib2')

    frame = name_frame(int(frame)+offset)
    idx_url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.'+datestr+'/hiresw.t'+cycle+'z.arw_5km.f'+frame+'.conusmem2.grib2.idx'
    os.system('curl "'+idx_url+'" --output "/root/hiresw.t'+cycle+'z.arw_5km.f'+frame+'.conusmem2.grib2.idx"')
    idx_file = '/root/hiresw.t'+cycle+'z.arw_5km.f'+frame+'.conusmem2.grib2.idx'
    read_idx(idx_file,'arw5k_2',int(frame),cycle,datestr)
    (xr.load_dataset('/root/current.grib2')).to_netcdf('/root/current.nc')
    os.system('/root/anaconda3/envs/blend/bin/gdalwarp -t_srs EPSG:4326 /root/current.nc /root/current_.tif')
    inputfile = '/root/current_.tif'
    outputfile = '/root/current_.nc'
    ds = gdal.Translate(outputfile, inputfile, format='NetCDF')
    os.remove(idx_file)
    dataset = xr.load_dataset('/root/current_.nc')
    if 'crs' in str(dataset):
        dataset = dataset.drop(['crs'])

    dataset['lon'] = dataset['lon']+360
    dataset = crop_ds(dataset,'180_chelsa')
    dataset = dataset.interp(lat=chelsa_ds["lat"], lon=chelsa_ds["lon"])
    dataset['tp'] = dataset['tp']*chelsa_ds['precip']

    return dataset

def fv35k(chelsa_ds,frame,cycle,datestr,offset):

    #remove files
    file_exists = os.path.exists('/root/current_.grib2')
    if file_exists == True:
        os.remove('/root/current_.grib2')
    file_exists = os.path.exists('/root/minus_one_.grib2')
    if file_exists == True:
        os.remove('/root/minus_one_.grib2')

    frame = name_frame(int(frame)+offset)
    idx_url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.'+datestr+'/hiresw.t'+cycle+'z.fv3_5km.f'+frame+'.conus.grib2.idx'
    os.system('curl "'+idx_url+'" --output "/root/hiresw.t'+cycle+'z.fv3_5km.f'+frame+'.conus.grib2.idx"')
    idx_file = '/root/hiresw.t'+cycle+'z.fv3_5km.f'+frame+'.conus.grib2.idx'
    read_idx(idx_file,'fv35k',int(frame),cycle,datestr)
    (xr.load_dataset('/root/current.grib2')).to_netcdf('/root/current.nc')
    os.system('/root/anaconda3/envs/blend/bin/gdalwarp -t_srs EPSG:4326 /root/current.nc /root/current_.tif')
    inputfile = '/root/current_.tif'
    outputfile = '/root/current_.nc'
    ds = gdal.Translate(outputfile, inputfile, format='NetCDF')
    os.remove(idx_file)
    dataset = xr.load_dataset('/root/current_.nc')
    if 'crs' in str(dataset):
        dataset = dataset.drop(['crs'])

    dataset['lon'] = dataset['lon']+360
    dataset = crop_ds(dataset,'180_chelsa')
    dataset = dataset.interp(lat=chelsa_ds["lat"], lon=chelsa_ds["lon"])
    dataset['tp'] = dataset['tp']*chelsa_ds['precip']

    return dataset

def arw2p5k(chelsa_ds,frame,cycle,datestr,offset):

    #remove files
    file_exists = os.path.exists('/root/current_.grib2')
    if file_exists == True:
        os.remove('/root/current_.grib2')
    file_exists = os.path.exists('/root/minus_one_.grib2')
    if file_exists == True:
        os.remove('/root/minus_one_.grib2')

    frame = name_frame(int(frame)+offset)
    idx_url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.'+datestr+'/hiresw.t'+cycle+'z.arw_2p5km.f'+frame+'.conus.grib2.idx'
    os.system('curl "'+idx_url+'" --output "/root/hiresw.t'+cycle+'z.arw_2p5km.f'+frame+'.conus.grib2.idx"')
    idx_file = '/root/hiresw.t'+cycle+'z.arw_2p5km.f'+frame+'.conus.grib2.idx'
    read_idx(idx_file,'arw2p5k',int(frame),cycle,datestr)
    (xr.load_dataset('/root/current.grib2')).to_netcdf('/root/current.nc')
    os.system('/root/anaconda3/envs/blend/bin/gdalwarp -t_srs EPSG:4326 /root/current.nc /root/current_.tif')
    inputfile = '/root/current_.tif'
    outputfile = '/root/current_.nc'
    ds = gdal.Translate(outputfile, inputfile, format='NetCDF')
    os.remove(idx_file)
    dataset = xr.load_dataset('/root/current_.nc')
    if 'crs' in str(dataset):
        dataset = dataset.drop(['crs'])

    dataset['lon'] = dataset['lon']+360
    dataset = crop_ds(dataset,'180_chelsa')
    dataset = dataset.interp(lat=chelsa_ds["lat"], lon=chelsa_ds["lon"])
    dataset['tp'] = dataset['tp']*chelsa_ds['precip']

    return dataset


def fv32p5k(chelsa_ds,frame,cycle,datestr,offset):

    #remove files
    file_exists = os.path.exists('/root/current_.grib2')
    if file_exists == True:
        os.remove('/root/current_.grib2')
    file_exists = os.path.exists('/root/minus_one_.grib2')
    if file_exists == True:
        os.remove('/root/minus_one_.grib2')

    frame = name_frame(int(frame)+offset)
    idx_url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/hiresw/prod/hiresw.'+datestr+'/hiresw.t'+cycle+'z.fv3_2p5km.f'+frame+'.conus.grib2.idx'
    os.system('curl "'+idx_url+'" --output "/root/hiresw.t'+cycle+'z.fv3_2p5km.f'+frame+'.conus.grib2.idx"')
    idx_file = '/root/hiresw.t'+cycle+'z.fv3_2p5km.f'+frame+'.conus.grib2.idx'
    read_idx(idx_file,'fv32p5k',int(frame),cycle,datestr)
    (xr.load_dataset('/root/current.grib2')).to_netcdf('/root/current.nc')
    os.system('/root/anaconda3/envs/blend/bin/gdalwarp -t_srs EPSG:4326 /root/current.nc /root/current_.tif')
    inputfile = '/root/current_.tif'
    outputfile = '/root/current_.nc'
    ds = gdal.Translate(outputfile, inputfile, format='NetCDF')
    os.remove(idx_file)
    dataset = xr.load_dataset('/root/current_.nc')
    if 'crs' in str(dataset):
        dataset = dataset.drop(['crs'])

    dataset['lon'] = dataset['lon']+360
    dataset = crop_ds(dataset,'180_chelsa')
    dataset = dataset.interp(lat=chelsa_ds["lat"], lon=chelsa_ds["lon"])
    dataset['tp'] = dataset['tp']*chelsa_ds['precip']

    return dataset





def crop_ds(ds,type):
    # if type == '180_chelsa':
        # for n in range(len(ds.lat)):
        #     print(max(ds.tp[n].values))
    # topleft_bottomright = [45,-125,35,-115]
    # topleft_bottomright = [41,-109,37,-102]
    topleft_bottomright = [50,-125,46.5,-120]
    # topleft_bottomright = []
    # topleft_bottomright = [50,-125,25,-60]
    # topleft_bottomright = [50,-125,45,-120]
    # topleft_bottomright = [41,-115,37,-102]

    if type == '360_grib' or type == '180_chelsa':
        min_lon = topleft_bottomright[1]+360
        max_lon = topleft_bottomright[3]+360
    elif type == '360_chelsa':
        min_lon = topleft_bottomright[1]+180
        max_lon = topleft_bottomright[3]+180
    elif type == 'resolutions':
        min_lon = topleft_bottomright[1]+180
        max_lon = topleft_bottomright[3]+180
    else:
        min_lon = topleft_bottomright[1]
        max_lon = topleft_bottomright[3]
    min_lat = topleft_bottomright[2]
    max_lat = topleft_bottomright[0]

    if type == '360_grib':
        mask_lon = (ds.longitude >= min_lon) & (ds.longitude <= max_lon)
        mask_lat = (ds.latitude >= min_lat) & (ds.latitude <= max_lat)
    elif type == '360_chelsa' or type == '180_chelsa' or type == 'resolutions':
        mask_lon = (ds.lon >= min_lon) & (ds.lon <= max_lon)
        mask_lat = (ds.lat >= min_lat) & (ds.lat <= max_lat)
    ds = ds.where(mask_lat, drop=True)
    ds = ds.where(mask_lon, drop=True)
    if type == '360_chelsa':
        ds['lon'] = ds['lon']+180
    if type=='180_chelsa':
        for n in range(len(ds.lat)):
            values = ds.tp[n].values
            output = []
            for value in zip(values):
                if str(value[0]) == 'nan':
                    value = 0
                else:
                    value = value[0]
                output.append(value)
            ds['tp'][n] = output
    return ds

def resolutions():
    resolutions = [3,5,2.5,12]
    for resolution in resolutions:
        ds = xr.load_dataset('/root/chelsa2.nc')
        latitude = 45
        km = math.cos(latitude*0.0174533)*111.321543
        degrees = resolution/km
        print(degrees)

        coarsen_resolution = degrees
        lons = int(43200/(360/coarsen_resolution))
        lats = int(20880/(180/coarsen_resolution))

        ds['lon'] = ds['lon']+180
        ds = crop_ds(ds,'resolutions')

        ds_coarse = ds.coarsen(lon=lons, lat=lats, boundary='pad').mean()

        ds_coarse = ds_coarse.interp(lat=ds["lat"], lon=ds["lon"])

        ds['precip'] = ds['precip']/ds_coarse['precip']

        for x in range(len(ds.lat)):
            ds_list = ds.precip[x].values
            output = []
            for value in zip(ds_list):
                value = value[0]
                if str(value) == 'nan':
                    value = 1
                output.append(value)
            ds.precip[x] = output

        ds.to_netcdf('/root/'+str(resolution)+'chelsa.nc')

def create_master_ds():
    datestr = datestr_and_cycle()[0]
    cycle = datestr_and_cycle()[1]
    sub_models = [['nam',3]]
    for model in sub_models:
        file_exists = os.path.exists('/root/current_.grib2')
        if file_exists == True:
            os.remove('/root/current_.grib2')
        file_exists = os.path.exists('/root/current.grib2')
        if file_exists == True:
            os.remove('/root/current.grib2')
        file_exists = os.path.exists('/root/master.grib2')
        if file_exists == True:
            os.remove('/root/master.grib2')
        if model[0] == 'nam':

            frame = 3

            frame = name_frame(int(frame))
            idx_url = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/nam/prod/nam.'+datestr+'/nam.t'+cycle+'z.conusnest.hiresf'+frame+'.tm00.grib2.idx'
            os.system('curl "'+idx_url+'" --output "/root/nam.t'+cycle+'z.conusnest.hiresf'+frame+'.tm00.grib2.idx"')
            idx_file = '/root/nam.t'+cycle+'z.conusnest.hiresf'+frame+'.tm00.grib2.idx'
            read_idx(idx_file,'nam',int(frame),cycle,datestr)
            (xr.load_dataset('/root/current.grib2')).to_netcdf('/root/current.nc')
            os.system('/root/anaconda3/envs/blend/bin/gdalwarp -t_srs EPSG:4326 /root/current.nc /root/master.tif')
            inputfile = '/root/master.tif'
            outputfile = '/root/master.nc'
            ds = gdal.Translate(outputfile, inputfile, format='NetCDF')
            os.remove(idx_file)
            dataset = xr.load_dataset('/root/master.nc')
            if 'crs' in str(dataset):
                dataset = dataset.drop(['crs'])

        # os.system('/root/anaconda3/envs/blend/bin/gdalwarp -t_srs EPSG:4326 current.grib2 current_.grib2')
        # dataset = xr.load_dataset('/root/current_.grib2',engine='cfgrib')
        dataset['lon'] = dataset['lon']+360
        dataset = crop_ds(dataset,'180_chelsa')
        dataset['tp'] = dataset['tp']*0
        chelsa_ds = xr.load_dataset('/root/3chelsa.nc')
        # chelsa_ds['lon'] = chelsa_ds['lon']-180
        # print(chelsa_ds,dataset)
        chelsa_ds = crop_ds(chelsa_ds,'360_chelsa')

        new_lon = np.linspace(chelsa_ds.lon[0], chelsa_ds.lon[-1], chelsa_ds.dims["lon"] * 2)
        new_lat = np.linspace(chelsa_ds.lat[0], chelsa_ds.lat[-1], chelsa_ds.dims["lat"] * 2)
        chelsa_ds = chelsa_ds.interp(lat=new_lat, lon=new_lon)

        dataset = dataset.interp(lat=chelsa_ds["lat"], lon=chelsa_ds["lon"])
        dataset['nam3k_1'] = dataset['tp']
        dataset['nam3k_2'] = dataset['tp']
        dataset['nam3k_3'] = dataset['tp']
        dataset['nam3k_4'] = dataset['tp']
        dataset['nam3k_5'] = dataset['tp']
        dataset['hrrr3k_1'] = dataset['tp']
        dataset['hrrr3k_2'] = dataset['tp']
        dataset['hrrr3k_3'] = dataset['tp']
        dataset['arw5k_1_1'] = dataset['tp']
        dataset['arw5k_1_2'] = dataset['tp']
        dataset['arw5k_2_1'] = dataset['tp']
        dataset['arw5k_2_2'] = dataset['tp']
        dataset['fv35k_1'] = dataset['tp']
        dataset['fv35k_2'] = dataset['tp']
        dataset['fv35k_3'] = dataset['tp']
        dataset['arw2.5k_1'] = dataset['tp']
        dataset['arw2.5k_2'] = dataset['tp']
        dataset['fv32.5k_1'] = dataset['tp']
        dataset['fv32.5k_2'] = dataset['tp']
        dataset['fv32.5k_3'] = dataset['tp']

    return dataset


def ingest_gribs(frame,master_ds):
    frame = name_frame(int(frame))
    datestr = datestr_and_cycle()[0]
    cycle = datestr_and_cycle()[1]
    # resolutions = [3,5,2.5]
    resolutions = [3]
    for resolution in resolutions:
        if resolution == 3:
            #load downscaling file
            chelsa_ds = xr.load_dataset('/root/'+str(resolution)+'chelsa.nc')
            # chelsa_ds['lon'] = chelsa_ds['lon']-180
            chelsa_ds = crop_ds(chelsa_ds,'360_chelsa')

            new_lon = np.linspace(chelsa_ds.lon[0], chelsa_ds.lon[-1], chelsa_ds.dims["lon"] * 2)
            new_lat = np.linspace(chelsa_ds.lat[0], chelsa_ds.lat[-1], chelsa_ds.dims["lat"] * 2)
            chelsa_ds = chelsa_ds.interp(lat=new_lat, lon=new_lon)

            # sub_models = [['nam3k_1',3],['nam3k_2',3],['hrrr3k',3]]
            # sub_models = [['nam3k_1',3],['nam3k_2',3],['nam3k_3',3],['nam3k_4',3],['nam3k_5',3],['nam3k_6',3],['nam3k_7',3],['nam3k_8',3]]
            sub_models = [['hrrr3k',3],['nam3k',3]]
            for model in sub_models:

                datestr = datestr_and_cycle()[0]
                cycle = datestr_and_cycle()[1]

                #remove files
                file_exists = os.path.exists('/root/current_.grib2')
                if file_exists == True:
                    os.remove('/root/current_.grib2')
                file_exists = os.path.exists('/root/minus_one_.grib2')
                if file_exists == True:
                    os.remove('/root/minus_one_.grib2')

                datasets = []
                #nam ingest
                if model[0] == 'nam3k':
                    if cycle == '00':
                        dataset_one = nam3k(chelsa_ds,frame,'00',datestr,0)
                        datestr = ((datetime.strptime(datestr, '%Y%m%d'))-timedelta(days=1)).strftime('%Y%m%d')
                        dataset_two = nam3k(chelsa_ds,frame,'18',datestr,6)
                        dataset_three = nam3k(chelsa_ds,frame,'12',datestr,12)
                        dataset_four = nam3k(chelsa_ds,frame,'06',datestr,18)
                        dataset_five = nam3k(chelsa_ds,frame,'00',datestr,24)
                    else:
                        dataset_one = nam3k(chelsa_ds,frame,'12',datestr,0)
                        dataset_two = nam3k(chelsa_ds,frame,'06',datestr,6)
                        dataset_three = nam3k(chelsa_ds,frame,'00',datestr,12)
                        datestr = ((datetime.strptime(datestr, '%Y%m%d'))-timedelta(days=1)).strftime('%Y%m%d')
                        dataset_four = nam3k(chelsa_ds,frame,'18',datestr,18)
                        dataset_five = nam3k(chelsa_ds,frame,'12',datestr,24)

                    datasets = [dataset_one,dataset_two,dataset_three,dataset_four,dataset_five]

                #hrrr ingest
                elif model[0] == 'hrrr3k':
                    if cycle == '00':
                        dataset_one = hrrr3k(chelsa_ds,frame,'00',datestr,0)
                        datestr = ((datetime.strptime(datestr, '%Y%m%d'))-timedelta(days=1)).strftime('%Y%m%d')
                        dataset_two = hrrr3k(chelsa_ds,frame,'18',datestr,6)
                        dataset_three = hrrr3k(chelsa_ds,frame,'12',datestr,12)
                    else:
                        dataset_one = hrrr3k(chelsa_ds,frame,'12',datestr,0)
                        dataset_two = hrrr3k(chelsa_ds,frame,'06',datestr,6)
                        dataset_three = hrrr3k(chelsa_ds,frame,'00',datestr,12)

                    datasets = [dataset_one,dataset_two,dataset_three]

                if model[0] == 'nam3k':
                    r = 5
                elif model[0] == 'hrrr3k':
                    r = 3
                for n in range(r):
                    master_ds[model[0]+'_'+str(n+1)] = master_ds[model[0]+'_'+str(n+1)]+datasets[n]['tp']

        elif resolution == 5:
            #load downscaling file
            chelsa_ds = xr.load_dataset('/root/'+str(resolution)+'chelsa.nc')
            # chelsa_ds['lon'] = chelsa_ds['lon']-180
            chelsa_ds = crop_ds(chelsa_ds,'360_chelsa')

            new_lon = np.linspace(chelsa_ds.lon[0], chelsa_ds.lon[-1], chelsa_ds.dims["lon"] * 2)
            new_lat = np.linspace(chelsa_ds.lat[0], chelsa_ds.lat[-1], chelsa_ds.dims["lat"] * 2)
            chelsa_ds = chelsa_ds.interp(lat=new_lat, lon=new_lon)

            sub_models = [['arw5k_1',5],['arw5k_2',5],['fv35k',5]]
            for model in sub_models:

                datestr = datestr_and_cycle()[0]
                cycle = datestr_and_cycle()[1]

                #remove files
                file_exists = os.path.exists('/root/current_.grib2')
                if file_exists == True:
                    os.remove('/root/current_.grib2')
                file_exists = os.path.exists('/root/minus_one_.grib2')
                if file_exists == True:
                    os.remove('/root/minus_one_.grib2')

                datasets = []
                #wrf arw5k 1 ingest
                if model[0] == 'arw5k_1':
                    if cycle == '00':
                        dataset_one = arw5k_1(chelsa_ds,frame,'00',datestr,0)
                        datestr = ((datetime.strptime(datestr, '%Y%m%d'))-timedelta(days=1)).strftime('%Y%m%d')
                        dataset_two = arw5k_1(chelsa_ds,frame,'12',datestr,12)
                    else:
                        dataset_one = arw5k_1(chelsa_ds,frame,'12',datestr,0)
                        dataset_two = arw5k_1(chelsa_ds,frame,'00',datestr,12)

                    datasets = [dataset_one,dataset_two]

                #wrf arw5k 2 ingest
                elif model[0] == 'arw5k_2':
                    if cycle == '00':
                        dataset_one = arw5k_2(chelsa_ds,frame,'00',datestr,0)
                        datestr = ((datetime.strptime(datestr, '%Y%m%d'))-timedelta(days=1)).strftime('%Y%m%d')
                        dataset_two = arw5k_2(chelsa_ds,frame,'12',datestr,12)
                    else:
                        dataset_one = arw5k_2(chelsa_ds,frame,'12',datestr,0)
                        dataset_two = arw5k_2(chelsa_ds,frame,'00',datestr,12)

                    datasets = [dataset_one,dataset_two]

                #wrf fv3 ingest
                elif model[0] == 'fv35k':
                    if cycle == '00':
                        dataset_one = fv35k(chelsa_ds,frame,'00',datestr,0)
                        datestr = ((datetime.strptime(datestr, '%Y%m%d'))-timedelta(days=1)).strftime('%Y%m%d')
                        dataset_two = fv35k(chelsa_ds,frame,'12',datestr,12)
                        dataset_three = fv35k(chelsa_ds,frame,'00',datestr,24)
                    else:
                        dataset_one = fv35k(chelsa_ds,frame,'12',datestr,0)
                        dataset_two = fv35k(chelsa_ds,frame,'00',datestr,12)
                        datestr = ((datetime.strptime(datestr, '%Y%m%d'))-timedelta(days=1)).strftime('%Y%m%d')
                        dataset_three = fv35k(chelsa_ds,frame,'12',datestr,24)

                    datasets = [dataset_one,dataset_two,dataset_three]

                if model[0] == 'arw5k_1' or model[0] == 'arw5k_2':
                    r = 2
                elif model[0] == 'fv35k':
                    r = 3
                for n in range(r):
                    master_ds[model[0]+'_'+str(n+1)] = master_ds[model[0]+'_'+str(n+1)]+datasets[n]['tp']

        elif resolution == 2.5:
            #load downscaling file
            chelsa_ds = xr.load_dataset('/root/'+str(resolution)+'chelsa.nc')
            # chelsa_ds['lon'] = chelsa_ds['lon']-180
            chelsa_ds = crop_ds(chelsa_ds,'360_chelsa')

            new_lon = np.linspace(chelsa_ds.lon[0], chelsa_ds.lon[-1], chelsa_ds.dims["lon"] * 2)
            new_lat = np.linspace(chelsa_ds.lat[0], chelsa_ds.lat[-1], chelsa_ds.dims["lat"] * 2)
            chelsa_ds = chelsa_ds.interp(lat=new_lat, lon=new_lon)

            sub_models = [['arw2.5k',2.5],['fv32.5k',2.5]]
            # sub_models = [['arw2.5k',2.5]]
            for model in sub_models:

                datestr = datestr_and_cycle()[0]
                cycle = datestr_and_cycle()[1]

                #remove files
                file_exists = os.path.exists('/root/current_.grib2')
                if file_exists == True:
                    os.remove('/root/current_.grib2')
                file_exists = os.path.exists('/root/minus_one_.grib2')
                if file_exists == True:
                    os.remove('/root/minus_one_.grib2')

                datasets = []
                #wrf arw2.5k 1 ingest
                if model[0] == 'arw2.5k':
                    if cycle == '00':
                        dataset_one = arw2p5k(chelsa_ds,frame,'00',datestr,0)
                        datestr = ((datetime.strptime(datestr, '%Y%m%d'))-timedelta(days=1)).strftime('%Y%m%d')
                        dataset_two = arw2p5k(chelsa_ds,frame,'12',datestr,12)
                    else:
                        dataset_one = arw2p5k(chelsa_ds,frame,'12',datestr,0)
                        dataset_two = arw2p5k(chelsa_ds,frame,'00',datestr,12)

                    datasets = [dataset_one,dataset_two]


                #wrf fv32.5k ingest
                elif model[0] == 'fv32.5k':
                    if cycle == '00':
                        dataset_one = fv32p5k(chelsa_ds,frame,'00',datestr,0)
                        datestr = ((datetime.strptime(datestr, '%Y%m%d'))-timedelta(days=1)).strftime('%Y%m%d')
                        dataset_two = fv32p5k(chelsa_ds,frame,'12',datestr,12)
                        dataset_three = fv32p5k(chelsa_ds,frame,'00',datestr,24)
                    else:
                        dataset_one = fv32p5k(chelsa_ds,frame,'12',datestr,0)
                        dataset_two = fv32p5k(chelsa_ds,frame,'00',datestr,12)
                        datestr = ((datetime.strptime(datestr, '%Y%m%d'))-timedelta(days=1)).strftime('%Y%m%d')
                        dataset_three = fv32p5k(chelsa_ds,frame,'12',datestr,24)

                    datasets = [dataset_one,dataset_two,dataset_three]


                if model[0] == 'arw2.5k':
                    r = 2
                elif model[0] == 'fv32.5k':
                    r = 3
                for n in range(r):
                    master_ds[model[0]+'_'+str(n+1)] = master_ds[model[0]+'_'+str(n+1)]+datasets[n]['tp']


    return master_ds

# resolutions()

frame = '03'
master_master_ds = create_master_ds()
# master_ds = create_master_ds()
# print(master_ds)
for n in range(8,36):
    master_ds = create_master_ds()
    frame = name_frame(n)
    master_ds = ingest_gribs(frame,master_ds)
    # master_ds['tp'] = (master_ds['nam3k']+master_ds['hrrr3k']+master_ds['arw5k_1']+master_ds['arw5k_2']+master_ds['fv35k']+master_ds['arw2.5k']+master_ds['fv32.5k'])/7
    # master_ds['tp'] = (master_ds['nam3k_1']+master_ds['nam3k_2']+master_ds['nam3k_3']+master_ds['nam3k_4']+master_ds['nam3k_5']+master_ds['hrrr3k_1']+master_ds['hrrr3k_2']+master_ds['hrrr3k_3']+master_ds['arw5k_1_1']+master_ds['arw5k_1_2']+master_ds['arw5k_2_1']+master_ds['arw5k_2_2']+master_ds['fv35k_1']+master_ds['fv35k_2']+master_ds['fv35k_3']+master_ds['arw2.5k_1']+master_ds['arw2.5k_2']+master_ds['fv32.5k_1']+master_ds['fv32.5k_2']+master_ds['fv32.5k_3'])/20
    master_ds['tp'] = (master_ds['nam3k_1']+master_ds['nam3k_2']+master_ds['hrrr3k_1']+master_ds['hrrr3k_2'])/4
    # master_ds['tp'] = (master_ds['nam3k_1']+master_ds['hrrr3k_1']+master_ds['arw5k_1_1']+master_ds['arw5k_2_1']+master_ds['fv35k_1']+master_ds['arw2.5k_1']+master_ds['fv32.5k_1'])/7
    master_ds.to_netcdf('/root/master_ds.nc')
    print(master_ds)

    # master_master_ds = xr.concat([master_master_ds,master_ds], dim="hour")

    # ds = master_master_ds.isel(hour=n-1)

    # new_lon = np.linspace(ds.lon[0], ds.lon[-1], ds.dims["lon"] * 2)
    # new_lat = np.linspace(ds.lat[0], ds.lat[-1], ds.dims["lat"] * 2)
    # ds = ds.interp(lat=new_lat, lon=new_lon)

    ds = master_ds
    # ds =

    lats = ds['lat']
    lons = ds['lon']
    # if int(frame) == 2:
    #     tp = ds['hrrr3k_1']*.0393701
    # elif int(frame) == 3:
    #     tp = ds['nam3k_1']*.0393701
    # elif int(frame) == 4:
    #     tp = ds['arw5k_1_1']*.0393701
    # elif int(frame) == 5:
    #     tp = ds['arw5k_2_1']*.0393701
    # #**********
    # elif int(frame) == 6:
    #     tp = ds['fv35k_1']*.0393701
    # #**********
    # elif int(frame) == 7:
    #     tp = ds['arw2.5k_1']*.0393701
    # elif int(frame) == 8:
    #     tp = ds['fv32.5k_1']*.0393701
    # elif int(frame) == 9:
    #     tp = ds['tp']*.0393701

    tp = ds['nam3k_1']*.0393701

    fig = plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    newcmp = create_colormap()
    bounds = [0,0.01,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,5,6,7,8,9,10]#,11,12,13,14,15,16,17,18,19,20]
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=32)
    #cf = ax.contourf(lons,lats,precip,us,norm=norm,cmap=newcmp)
    cf = ax.pcolormesh(lons, lats, tp, norm=norm, cmap=newcmp)

    # cf = ax.pcolormesh(lons, lats, tp, cmap='jet', vmin=0, vmax=50)
    ax.coastlines()
    ax.add_feature(cartopy.feature.STATES)
    ax.add_feature(USCOUNTIES.with_scale('500k'),linewidth=1)
    cbar = plt.colorbar(cf, shrink=0.7, orientation="horizontal", pad=0.03)
    plt.savefig('/root/script/hrcamef/tp_'+frame+'.png',dpi=500,bbox_inches='tight')
    plt.clf()

    os.chdir('/root/script')
    os.system('git add hrcamef')
    os.system('git commit -m "auto-push"')
    os.system('git checkout master')
    os.system('git pull git@github.com:clamalo/ubuntu_scripts.git master')
    os.system('git config --global core.askpass "git-gui--askpass"')
    os.system('git push git@github.com:clamalo/ubuntu_scripts.git master')
