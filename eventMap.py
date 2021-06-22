"""
Before running the program, enter your specifications into
SETUP_eventMap.txt
"""
#########################################################################
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4
import numpy as np
import obspy
from obspy.clients.fdsn.client import Client
from obspy.core import UTCDateTime
from obspy.core.event import Event
from obspy.core.event import FocalMechanism
from obspy.core.event import Origin
from obspy.imaging.beachball import beach
import os

#### Read program specs from the setup file
with open('SETUP/SETUP_eventMap.txt', 'r') as f:
    setup=f.read().splitlines()
PROJECTION=setup[0][11:]
MIN_LATITUDE=float(setup[1][13:])
MAX_LATITUDE=float(setup[2][13:])
MIN_LONGITUDE=float(setup[3][14:])
MAX_LONGITUDE=float(setup[4][14:])
QUAKEML_FILE=setup[5][13:]
NETCDF4_FILE=setup[6][13:]
IMAGE_FILE=setup[7][11:]

start=UTCDateTime() #time the program

#### Create basemap
def polar_stere(lon_w, lon_e, lat_s, lat_n, **kwargs):
    ## Returns a basemap in polar stereographic of a specific box not centered at a pole
    lon_0 = lon_w + (lon_e - lon_w) / 2.
    ref = lat_s if abs(lat_s) > abs(lat_n) else lat_n
    lat_0 = np.copysign(90., ref)
    prj = Basemap(projection=PROJECTION, lon_0=lon_0, lat_0=lat_0,
                          boundinglat=0, resolution='l')
    lons = [lon_w, lon_e, lon_w, lon_e, lon_0, lon_0]
    lats = [lat_s, lat_s, lat_n, lat_n, lat_s, lat_n]
    x, y = prj(lons, lats)
    ll_lon, ll_lat = prj(min(x), min(y), inverse=True)
    ur_lon, ur_lat = prj(max(x), max(y), inverse=True)
    return Basemap(projection='stere', lat_0=lat_0, lon_0=lon_0,
                           llcrnrlon=ll_lon, llcrnrlat=ll_lat,
                           urcrnrlon=ur_lon, urcrnrlat=ur_lat, **kwargs)
if (PROJECTION == 'spstere' or PROJECTION == 'npstere'):
    my_map = polar_stere(MIN_LONGITUDE, MAX_LONGITUDE, MIN_LATITUDE, MAX_LATITUDE, area_thresh=0.1)
else:
    my_map = Basemap(projection=PROJECTION, lat_0=MIN_LATITUDE, 
            lon_0=MIN_LONGITUDE, resolution='l', area_thresh=0.1,
            llcrnrlon=MIN_LONGITUDE, llcrnrlat=MIN_LATITUDE,
            urcrnrlon=MAX_LONGITUDE, urcrnrlat=MAX_LATITUDE)
my_map.drawcoastlines()
my_map.drawcountries()
my_map.fillcontinents(color=[1,1,1])
my_map.drawmapboundary()
my_map.drawmeridians(range(120,129,8), labels=[0,0,0,1], fontsize=3)

#### Read events from saved file from 'fetchdata.py'
events = obspy.read_events(QUAKEML_FILE, 'QUAKEML')
lons, lats, deps, mags = [], [], [], []
for event in events:
    if (len(event.focal_mechanisms) > 0):
        lons.append(event.preferred_origin().longitude)
        lats.append(event.preferred_origin().latitude)
        deps.append(event.preferred_origin().depth)
        mags.append(event.magnitudes[0].mag)
print('\n___________________\nFirst event: '+str(events[0].resource_id)+
        '\ntime: '+str(events[0].preferred_origin().time)+'\nlat,lon: '+
        str(lats[0])+','+str(lons[0])+'\nmagnitude: '+str(mags[0])+
        '\ndepth: '+str(deps[0])+'\n')

#### Plot the focal mechanisms
x,y = my_map(lons, lats)
ax = plt.gca()
maxDepth = max(deps)
skips=0
for xxx, yyy, dep, mag in zip(x, y, deps, mags):
    try:
        my_map.plot(xxx, yyy, marker='o', markersize=1.5*mag, 
        markerfacecolor='none', markeredgecolor=[1,0,0], 
        markeredgewidth=2.0)
        #plt.text(xxx, yyy, str(dep), fontsize=2)
    except: skips+=1
print('\n...'+str(skips)+' events skipped...\n')
#### Read topobathymetry data from saved file
nc = netCDF4.Dataset(NETCDF4_FILE,'r',format='NETCDF4')
bathLat =  nc.variables['latitude'][::10]
bathLon = nc.variables['longitude'][::10]
bathLon, bathLat = np.meshgrid(bathLon, bathLat)
bathX, bathY = my_map(bathLon, bathLat)

#### Add bathymetry to map
#interval=range(np.amin(nc.variables['ROSE'][::10,::10]),0,25)
interval=range(-6000,25,25)
cs1 = my_map.contourf(bathX,bathY,nc.variables['ROSE'][::10,::10], interval)
cbar1 = my_map.colorbar(cs1,pad="5%")
cbar1.set_label('Ocean Depth (m)')

#### save the plot and display runtime
plt.savefig(IMAGE_FILE, dpi=2400)
print('\n\nRuntime: '+str(UTCDateTime()-start)+' seconds\n\n')
os.system('say "Event mapping complete"')
##########################################################################

"""
cite smith and sandwell (bathymetry)
cite polar_stere
"""