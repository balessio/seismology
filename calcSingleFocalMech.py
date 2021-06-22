"""
Before running the program, enter your specifications into
SETUP_focalMechMap.txt
"""
#########################################################################
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4
import numpy as np
import obspy
from obspy import read
from obspy.clients.fdsn.client import Client
from obspy.core import UTCDateTime
from obspy.core.event import Event
from obspy.core.event import FocalMechanism
from obspy.core.event import Origin
from obspy.core.stream import Stream
from obspy.core.trace import Trace
from obspy.imaging.beachball import beach
from obspy.geodetics.base import gps2dist_azimuth as separation
from obspy.signal.trigger import ar_pick
import os

#### Read program specs from the setup file
with open('SETUP/SETUP_calcSingleFocalMech.txt', 'r') as f:
    setup=f.read().splitlines()
PROJECTION=setup[0][11:]
MIN_LATITUDE=float(setup[1][13:])
MAX_LATITUDE=float(setup[2][13:])
MIN_LONGITUDE=float(setup[3][14:])
MAX_LONGITUDE=float(setup[4][14:])
QUAKEML_FILE=setup[5][13:]
NETCDF4_FILE=setup[6][13:]
IMAGE_PATH=setup[7][11:]
WAVE_PATH=setup[8][10:]

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
lons, lats, deps, mags, fmcs = [], [], [], [], []
for event in events:
    if (len(event.focal_mechanisms) > 0):
        lons.append(event.preferred_origin().longitude)
        lats.append(event.preferred_origin().latitude)
        deps.append(event.preferred_origin().depth)
        mags.append(event.magnitudes[0].mag)
        fmcs.append(event.focal_mechanisms[0])
print('\n___________________\nFirst event: '+str(events[0].resource_id)+
        '\ntime: '+str(events[0].preferred_origin().time)+'\nlat,lon: '+
        str(lats[0])+','+str(lons[0])+'\nmagnitude: '+str(mags[0])+
        '\ndepth: '+str(deps[0])+'\nfocal mec: '+
        str(fmcs[0].resource_id)+'\n___________________\n')

#### Plot the focal mechanisms
x,y = my_map(lons, lats)
ax = plt.gca()
maxDepth = max(deps)
skips=0
for xxx, yyy, dep, mag, fmec in zip(x, y, deps, mags, fmcs):
    try:
        t = fmec.moment_tensor.tensor
        print(t)
        b = beach([t.m_rr,t.m_tt,t.m_pp,t.m_rt,t.m_rp,t.m_tp], xy=(xxx, yyy), 
                width=100000*mag, linewidth=0.1, facecolor=[dep/maxDepth,0,0])
        b.set_zorder(10)
        ax.add_collection(b)
    except: skips+=1
print('\n...'+str(skips)+' events skipped...\n')
#### Read topobathymetry data from saved file
nc = netCDF4.Dataset(NETCDF4_FILE,'r',format='NETCDF4')
bathLat =  nc.variables['latitude'][::10]
bathLon = nc.variables['longitude'][::10]
bathLon,bathLat = np.meshgrid(bathLon,bathLat)
bathX, bathY = my_map(bathLon, bathLat)

#### Add bathymetry to map
#interval=range(np.amin(nc.variables['ROSE'][::10,::10]),0,25)
interval=range(-6000,25,25)
cs1 = my_map.contourf(bathX,bathY,nc.variables['ROSE'][::10,::10], interval)
cbar1 = my_map.colorbar(cs1,pad="5%")
cbar1.set_label('Ocean Depth (m)')

#### save the plot and display runtime
plt.savefig(IMAGE_PATH+'initial.png', dpi=2400)
print('\nfirst plot complete\n')

#### organize the wave data
client=Client('IRIS')
traces=[]
for mseed in os.listdir(WAVE_PATH): # reads all mseed files in path
        if '.mseed' in mseed:
            stream=read(WAVE_PATH+mseed)
            for trace in stream: traces.append(trace)
traceIDdict={}

for trace in traces:
    idList=trace.id.split('.')
    try: traceIDdict[idList[0]+idList[1]+idList[2]].append(trace)
    except: traceIDdict[idList[0]+idList[1]+idList[2]]=[trace]
print('\ndictionary filled\n')

firstArrivals=[] ## (Radial first arrival amplitude with polarity, distance^3, azimuth)
groupedStreams=[] ## (Z,N,E,R,T)
posZandAzm=[] ## (+stZ[0].max(),distance^3,azimuth)
for key in traceIDdict:
    if len(traceIDdict[key])==5:
        for trace in traceIDdict[key]:
            traceInfo=trace.id.split('.')
            if traceInfo[3][2] == 'Z':
                station=client.get_stations(network=traceInfo[0],
                                        station=traceInfo[1],level='response')
                coords=station.get_coordinates(seed_id=trace.id)
                sep=separation(lats[0],lons[0],coords['latitude'],coords['longitude'])
                stZ=Stream(traces=[trace])
            elif traceInfo[3][2] == 'E': stE=Stream(traces=[trace])
            elif traceInfo[3][2] == 'N': stN=Stream(traces=[trace])
            elif traceInfo[3][2] == 'R': stR=Stream(traces=[trace])
            elif traceInfo[3][2] == 'T': stT=Stream(traces=[trace])
        arrivals=ar_pick(a=stZ[0].data,b=stN[0].data,c=stE[0].data,
                samp_rate=stE[0].stats.sampling_rate,
                f1=0.02,f2=1.0,
                lta_p=60.0,sta_p=2.0,
                lta_s=60.0,sta_s=2.0,
                m_p=4,m_s=4,
                l_p=4,l_s=4,
                s_pick=True)
        firstPindex=int(np.round(arrivals[0]*stE[0].stats.sampling_rate))
        firstPamp=stR[0].data[firstPindex]
        firstArrivals.append((firstPamp,sep[0]**3,sep[1]))
        ZvR.append((abs(stZ[0].data[firstPindex]/firstPamp),1,sep[1]))
        #groupedStreams.append((stZ,stN,stE,stR,stT))
        zMax=stZ[0].max()
        if zMax>0: posZandAzm.append((zMax,sep[0]**3,sep[1]))
    else: 
        for trace in traceIDdict[key]:
            traceInfo=trace.id.split('.')
            if traceInfo[3][2] == 'Z':
                station=client.get_stations(network=traceInfo[0],
                                        station=traceInfo[1],level='response')
                coords=station.get_coordinates(seed_id=trace.id)
                sep=separation(lats[0],lons[0],coords['latitude'],coords['longitude'])
                stZ=Stream(traces=[trace])
                zMax=stZ[0].max()
                if zMax>0: posZandAzm.append((zMax,sep[0]**3,sep[1]))
print(firstArrivals[0])

#### write the first arrival amplitudes and locations into a text file
posRad,negRad=[],[] # separate the entries by direction of first motion
with open(IMAGE_PATH+'info.txt','w') as f:
    f.write('Radial first arrival amplitude, distance^3, azimuth\n\n')
    for entry in firstArrivals:
        f.write(str(entry)+'\n')
        if entry[0]>=0: posRad.append(entry)
        else: negRad.append(entry)

#### merge sort functions for tuples
def merge(a=[],b=[],index=0):
    c=[(Trace(),0.0)]*(len(a)+len(b))
    countC,countA,countB=0,0,0
    while (countC<len(c)):
        if (countA<len(a) and countB<len(b)):
            if (a[countA][index]<b[countB][index]):
                c[countC]=a[countA]
                countA+=1
                countC+=1
        if (countA<len(a) and countB<len(b)):
            if (a[countA][index]>b[countB][index]):
                c[countC]=b[countB]
                countB+=1
                countC+=1
        if (countA<len(a) and countB<len(b)):
            if (a[countA][index]==b[countB][index]):
                c[countC]=a[countA]
                countC+=1
                countA+=1
                c[countC]=b[countB]
                countC+=1
                countB+=1
        elif (countA<len(a)):
            while (countA<len(a)):
                c[countC]=a[countA]
                countA+=1
                countC+=1
        elif (countB<len(b)):
            while (countB<len(b)):
                c[countC]=b[countB]
                countB+=1
                countC+=1
    return c
def split_list(to_be_split=[]):
    if (len(to_be_split)==1): return to_be_split
    x=0
    a,b=[],[]
    split=int(len(to_be_split)/2)
    while (x<len(to_be_split)):
        if (x<split): a.append(to_be_split[x])
        else: b.append(to_be_split[x])
        x+=1
    return a,b
def merge_sort(to_be_merged=[],index=0):
    if (len(to_be_merged)<2): return to_be_merged
    left,right=split_list(to_be_merged)
    left=merge_sort(left,index)
    right=merge_sort(right,index)
    return merge(left,right,index)

#### recalculate the focal mechanism from the wave data
incAzmPosRad=merge_sort(posRad,2) ## 2 is the index of azimuth
incAzmNegRad=merge_sort(negRad,2)
incZ=merge_sort(posZandAzm,2)
inZvR=merge_sort(ZvR,2)

## returns the average amplitude of an index of tuples of a list between
## previously set range of azimuths, also dividing by the distance^3
def local_average(values=[],index=0):
    avg=0.0
    for value in values:
        avg+=abs(value[index])/value[1]
    return avg/len(values)
## returns an azimuth specifying orientation of focal mechanism
def search_azm(values=[],index=0):
    a,b=split_list(values)
    avgA=local_average(a,index)
    avgB=local_average(b,index)
    if len(a)>7 and len(b)>7:
        if avgA>avgB: azm=search_azm(a,index)
        else: azm=search_azm(b,index)
    else:
        if avgA>avgB: return a[0][2]+((a[len(a)-1][2]-a[0][2])/2)
        return b[0][2]+((b[len(b)-1][2]-b[0][2])/2)
    return azm

strikeAngle=search_azm(incAzmPosRad,0)
#negRazm=search_azm(incAzmNegRad,0)
slipAngle=search_azm(incZ,0)
if slipAngle>180: zBazm=slipAngle-180
else: zBazm=slipAngle+180
print('\nstrike angle: '+str(strikeAngle)+'\nslip angle: '+str(slipAngle))

## returns boolean for if within width number of degrees of an azimuth
def within(to_fit=0.0,azm=0.0,width=5):
    if azm<=360-width and azm>=width:
        return to_fit<=azm+width and to_fit>=azm-width
    elif azm>360-width:
        return to_fit>=azm-width or to_fit>azm or to_fit<=width+azm-360
    else:
        return to_fit<=azm or to_fit>=360+azm-width

## calculating the dip angle
ZatAzm,ZatBazm=[],[]
try:
    for entry in incZvR:
            if within(entry[2],slipAngle,80): ZatAzm.append(entry)
            elif within(entry[2],zBazm,80): ZatBazm.append(entry)
    dipAngle=90*abs(local_average(ZatBazm,0)/local_average(ZatAzm,0))
except:
    for entry in incZ:
            if within(entry[2],slipAngle,80): ZatAzm.append(entry)
            elif within(entry[2],zBazm,80): ZatBazm.append(entry)
    dipAngle=90*abs(local_average(ZatBazm,0)/local_average(ZatAzm,0))

## calculating the moment tensor
angles=np.radians((strikeAngle,slipAngle,dipAngle))
normal=(-1*np.sin(angles[2])*np.sin(angles[0]),-1*np.sin(angles[2])*np.cos(angles[0]),np.cos(angles[2]))
slip=(np.cos(angles[1])*np.cos(angles[0])+np.sin(angles[1])*np.cos(angles[2])*np.sin(angles[0]),
    -1*np.cos(angles[1])*np.sin(angles[0])+np.sin(angles[1])*np.cos(angles[2])*np.cos(angles[0]),
    np.sin(angles[1])*np.sin(angles[2]))

fmec.moment_tensor.tensor.m_rr=2*normal[0]*slip[0]
fmec.moment_tensor.tensor.m_tt=2*normal[1]*slip[1]
fmec.moment_tensor.tensor.m_pp=2*normal[2]*slip[2]
fmec.moment_tensor.tensor.m_rt=normal[0]*slip[1]+slip[0]*normal[1]
fmec.moment_tensor.tensor.m_rp=normal[0]*slip[2]+slip[0]*normal[2]
fmec.moment_tensor.tensor.m_tp=normal[1]*slip[2]+slip[1]*normal[2]
t = fmec.moment_tensor.tensor
print(t)
b = beach([t.m_rr,t.m_tt,t.m_pp,t.m_rt,t.m_rp,t.m_tp], xy=(xxx, yyy), 
        width=100000*mag, linewidth=0.1, facecolor=[dep/maxDepth,0,0])
b.set_zorder(10)
ax.add_collection(b)
plt.savefig(IMAGE_PATH+'_recalc.png', dpi=2400)

print('\n\nRuntime: '+str(UTCDateTime()-start)+' seconds\n\n')
os.system('say "Focal mechanism mapping complete"')
##########################################################################

"""
cite smith and sandwell (bathymetry)
cite polar_stere
"""