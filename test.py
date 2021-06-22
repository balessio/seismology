########################
"""
import matplotlib.pyplot as plt
from PIL import Image
plt.plot([1,2,3,4,7])
plt.savefig('plot.png', dpi = 1200)
Image.open('plot.png').show()
"""
#######################
"""
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
data = Image.open('s05_w173_1arc_v3.tif')
d = np.array(data)
print(np.where(d==1))
 """
#######################
"""
from obspy.clients.fdsn import Client
client = Client("IRIS")
from obspy import UTCDateTime
t = UTCDateTime("2010-02-27T06:45:00.000")
st = client.get_waveforms("IU", "MBWA", "00", "BHZ", t, t+3600)
#st+=client.get_waveforms("IU", "ADK", "00", "LHZ", t, t+3600)
print(st)
st.plot()
"""
####################### 
""" 
from obspy.core import UTCDateTime as dt
a=dt()
print(dt()-a)
"""
#######################
"""
Before running the program, check that it is using the right file set as
in lines 44 and 45
######################################################################
import obspy
from obspy.clients.fdsn.client import Client
from obspy.core import UTCDateTime
from obspy.core.event import Event
from obspy.core.event import FocalMechanism
from obspy.core.event import Origin

start=UTCDateTime() # Time the program

#### Download event data from ISC
client = Client('ISC')
now = UTCDateTime()
x=0
while (x > -1):
    now-=x
    timeBefore = now - 86400
    try:
        events = client.get_events(starttime=timeBefore, endtime=now, 
        minlatitude=-85, maxlatitude=0, 
        minlongitude=102, maxlongitude=175, 
        includeallmagnitudes=True, includeallorigins=True, minmagnitude=5)
        #### Check for usable data
        print('\n...'+str(len(events))+' events found...\n')
        cleanEvents=client.get_events(eventid='609301')
        cleanEvents.clear()
        numBrokenEvents=0
        for event in events:
            a=0.0
            b=1.0
            try:
                a=float(event.magnitudes[0].mag)
                a=event.focal_mechanisms[0].moment_tensor.tensor.m_rr
                a=event.preferred_origin()
                b/=float(a.depth)
                cleanEvents.append(event)
            except: numBrokenEvents+=1
        print('\n...Removed '+str(numBrokenEvents)+' broken events...\n')
        try: 
            print(cleanEvents[0].resource_id)
            x = -1
        except: x+=86400#raise Exception('\n\nError: There are no useful events!\n\n')
    except: x+=86400
#### Write it into the file
cleanEvents.write("dataRecent.xml", format="QUAKEML") 
with open('outputs/dataRecent.txt', 'w') as textFile:
    textFile.write(cleanEvents.__str__(print_all=True))

print('\n\nRuntime: '+str(UTCDateTime()-start)+' seconds\n\n')
######################################################################
 """
#######################
"""
from obspy import read
from obspy.core import UTCDateTime
import obspy
from obspy.clients.fdsn.client import Client
from obspy.core.event import Event
from obspy.core.event import Origin

start=UTCDateTime()
client=Client('IRIS')
station=client.get_stations(network='II',station='TAU',level='response')
asdf=station.get_coordinates(seed_id='II.TAU.00.BHZ')
print(asdf['elevation'])
print('\nRuntime: '+str(UTCDateTime()-start))
"""
#######################
"""
import numpy as np

a = np.array([0,1,2])
b=np.array([3,1,4])
print(np.linalg.norm(a-b))
print(np.sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+
        (a[2]-b[2])*(a[2]-b[2])))
"""
#######################
"""
from obspy import read, read_inventory
#from obspy.clients.fdsn.client import Client

#client=Client('IRIS')
st = read("inputs/MSEED/square_SEIR_waves_2002to2018/2017L-49,164M6.5/epDist/transRad/G.WUS.00.BHE.mseed")
st2=read("inputs/MSEED/square_SEIR_waves_2002to2018/2017L-49,164M6.5/epDist/transRad/G.WUS.00.BHN.mseed")
for trace in st2: st.append(trace)
#stations=[]
#for trace in st:
 #   traceInfo=trace.id.split('.')
    #stations.append(client.get_stations(network=traceInfo[0],station=traceInfo[1],level='response'))

st.rotate(method='NE->RT',back_azimuth=12.0)
st.plot()
"""
#######################
"""
dic={}
dic['hello']=(1.0,[1,2,3])
print('\n1')
print(dic)
dic['hello'][1][1]*=3
print('\n2')
print(dic)
dic['goodbye']=(2.0,[2,3,4])
dic['goodbye'][1].append(5)
print('\n3')
print(dic)
for key in dic: print(dic[key])
"""
#######################
"""
from obspy import read, read_inventory
from obspy.geodetics.base import gps2dist_azimuth as separation

st = read('inputs/MSEED/square_SEIR_waves_2002to2018/2017L-49,164M6.5/epDist/transRad/G.WUS.00.BHR.mseed')
st[0].stats.distance=100.0
print(separation(-49,164,41.200716,79.216498))
st.filter(type='bandpass',freqmin=0.02,freqmax=0.1)
#st.filter(type='lowpass',freq=0.02)
print(st[0].max())
#st.plot()
"""
######################
"""
from obspy.clients.fdsn.client import Client
from obspy.core.event import Event

client=Client('ISC')
event=client.get_events(eventid=609261308)
event.write('inputs/QUAKEML/event_in_2016.xml', format="QUAKEML")
"""
######################
"""
import numpy as np
import matplotlib.pyplot as plt # don't use pylab
import matplotlib.colors
import matplotlib.cm

maxRad=0.0
r=[]
theta=[]
colors=[]
for value in values: 
    if abs(value[0]) > maxRad: maxRad=abs(value[0])
    theta.append(value[1])
    if value[0]>0:
        colors.append('blue')
    else: colors.append('red')
for value in values: r.append(value[0]/maxRad)
fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8],polar=True)
ax.set_ylim(0,1)  
sc = ax.scatter(r,theta,c =colors, s=49, zorder=2)
#plt.polar(values,values,alpha=0.8)
plt.show()
"""
######################
"""
import numpy as np
import matplotlib.pyplot as plt # don't use pylab
import matplotlib.colors
import matplotlib.cm

dec = []
ra = []
n = []
maxD=0.0
for point in points:
    if abs(point[0])>maxD: maxD=abs(point[0])
    ra.append(point[1])
    if point[0]<=0: n.append(15)
    else: n.append(50)
for point in points: dec.append((90*abs(point[0]))/maxD)
ra = [x/180.0*np.pi for x in ra]
fig = plt.figure()
ax = fig.add_subplot(111,polar=True)
ax.set_ylim(0,90)

# 4. only show 0,30, 60 ticks
ax.set_yticks([0,30,60])
# 4. orient ylabels along horizontal line
ax.set_rlabel_position(0)

# 1. prepare cmap and norm
colors= ["red"] * 31 + ["gold"] * 10 + ["limegreen"] * 20
cmap=matplotlib.colors.ListedColormap(colors)
norm = matplotlib.colors.Normalize(vmin=0, vmax=60)   
# 2. make circles bigger, using `s` argument
# 1. set different colors according to `n`
sc = ax.scatter(ra,dec,c =n, s=49, cmap=cmap, norm=norm, zorder=2,alpha=0.5)

# 1. make colorbar
cax = fig.add_axes([0.8,0.1,0.01,0.2])
fig.colorbar(sc, cax=cax, label="n", ticks=[0,30,40,60])
# 3. move title upwards, then adjust top spacing
ax.set_title("Graph Title here", va='bottom', y=1.1)
plt.subplots_adjust(top=0.8)

plt.show()
"""
######################

from obspy import read, read_inventory
from obspy.geodetics.base import gps2dist_azimuth as separation
from obspy.signal.trigger import ar_pick
import numpy as np

stE = read('inputs/MSEED/square_SEIR_waves_2002to2018/2016L-49,126M6.0/S.AUTAR..BHE.mseed')[0]
stN = read('inputs/MSEED/square_SEIR_waves_2002to2018/2016L-49,126M6.0/S.AUTAR..BHN.mseed')[0]
stZ = read('inputs/MSEED/square_SEIR_waves_2002to2018/2016L-49,126M6.0/S.AUTAR..BHZ.mseed')[0]
stR = read('inputs/MSEED/square_SEIR_waves_2002to2018/2016L-49,126M6.0/S.AUTAR..BHR.mseed')[0]
stT = read('inputs/MSEED/square_SEIR_waves_2002to2018/2016L-49,126M6.0/S.AUTAR..BHT.mseed')[0]
dataE=stE.data
dataN=stN.data
dataZ=stZ.data
dataR=stR.data
dataT=stT.data
arrivals=ar_pick(a=dataZ,b=dataN,c=dataE,
                samp_rate=stE.stats.sampling_rate,
                f1=0.02,f2=1.0,
                lta_p=60.0,sta_p=2.0,
                lta_s=60.0,sta_s=2.0,
                m_p=4,m_s=4,
                l_p=4,l_s=4,
                s_pick=True)
firstPindex=int(np.round(arrivals[0]*stE.stats.sampling_rate))
firstPamp=dataR[firstPindex]
print(firstPamp)

######################


