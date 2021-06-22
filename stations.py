"""
Before running the program, enter your specifications into
SETUP_stations.txt
"""
###############################################################
import matplotlib.pyplot as plt
import numpy as np
from obspy import read
from obspy.core import UTCDateTime
import os
import obspy
from obspy.clients.fdsn.client import Client
from obspy.core.event import Event
from obspy.core.event import Origin
from obspy.core.stream import Stream
from obspy.core.trace import Trace
from obspy.geodetics.base import gps2dist_azimuth as separation

start = UTCDateTime() # time the program

#### Read program specs from the setup file
with open('SETUP/SETUP_stations.txt', 'r') as f:
    setup=f.read().splitlines()
IN_PATH=setup[0][8:]
OUT_PATH=setup[1][9:]
EVENT_ID=setup[2][9:]

client=Client('ISC')
event=client.get_events(eventid=EVENT_ID)[0]
print(str(event.preferred_origin().longitude)+','+str(event.preferred_origin().latitude))

client=Client('IRIS')
failures=0
IiIuG=[]
AuSIS=[]
RT_IiIuG=[]
RT_AuSIS=[]
with open(OUT_PATH+'info.txt', 'w') as textFile:
    for mseed in os.listdir(IN_PATH): # reads all mseed files in path
        try:
            if '.mseed' in mseed:
                stream = read(IN_PATH+mseed)
                textFile.write(stream.__str__(extended=True)+'\n')
                for trace in stream:
                    traceInfo=trace.id.split('.')
                    station=client.get_stations(network=traceInfo[0],
                                                station=traceInfo[1],level='response')
                    #textFile.write(station.get_response(seedid=trace.id))
                    try:
                        station.plot_response(min_freq=0.001,network=traceInfo[0],
                                                station=traceInfo[1],
                                                location=traceInfo[2],channel=traceInfo[3],
                                                outfile=OUT_PATH+'response'+trace.id+'.pdf')
                    except: print('skipping response of '+trace.id+'\n')
                    try:
                        if (traceInfo[0]=='S'):
                            AuSIS.append((trace,station.get_coordinates(seed_id=trace.id)))
                        else: IiIuG.append((trace,station.get_coordinates(seed_id=trace.id)))
                    except: 
                        if (traceInfo[0]=='S'):
                            RT_AuSIS.append((trace,traceInfo))
                        else: RT_IiIuG.append((trace,traceInfo))
        except: failures+=1
print('\nFailed '+str(failures)+' time(s)\n')

#### fill in the epicentral distance for any rotated traces which are missing it
fixes=[]
for trace in RT_IiIuG:
    for entry in IiIuG:
        matchID=entry[0].id.split('.')
        if trace[1][0]==matchID[0] and trace[1][1]==matchID[1] and trace[1][2]==matchID[2]:
            fixes.append((trace[0],entry[1]))
for entry in fixes: IiIuG.append(entry)
fixes=[]
for trace in RT_AuSIS:
    for entry in AuSIS:
        matchID=entry[0].id.split('.')
        if trace[1][0]==matchID[0] and trace[1][1]==matchID[1] and trace[1][2]==matchID[2]:
            fixes.append((trace[0],entry[1]))
for entry in fixes: AuSIS.append(entry)

#event=obspy.read_events('inputs/QUAKEML/event_in_2017.xml', 'QUAKEML')[0]
lonE=event.preferred_origin().longitude
latE=event.preferred_origin().latitude
altE=event.preferred_origin().depth * -1
for x in range(0,len(IiIuG)):
    IiIuG[x]=(IiIuG[x][0],separation(IiIuG[x][1]['latitude'],
                IiIuG[x][1]['longitude'],latE,lonE)[0])
for x in range(0,len(AuSIS)):
    AuSIS[x]=(AuSIS[x][0],separation(AuSIS[x][1]['latitude'],
                AuSIS[x][1]['longitude'],latE,lonE)[0])

## merge sort functions
def merge(a=[],b=[]):
    c=[(Trace(),0.0)]*(len(a)+len(b))
    countC,countA,countB=0,0,0
    while (countC<len(c)):
        if (countA<len(a) and countB<len(b)):
            if (a[countA][1]<b[countB][1]):
                c[countC]=a[countA]
                countA+=1
                countC+=1
        if (countA<len(a) and countB<len(b)):
            if (a[countA][1]>b[countB][1]):
                c[countC]=b[countB]
                countB+=1
                countC+=1
        if (countA<len(a) and countB<len(b)):
            if (a[countA][1]==b[countB][1]):
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
def merge_sort(to_be_merged=[]):
    if (len(to_be_merged)<2): return to_be_merged
    left,right=split_list(to_be_merged)
    left=merge_sort(left)
    right=merge_sort(right)
    return merge(left,right)
IiIuG=merge_sort(IiIuG)
AuSIS=merge_sort(AuSIS)

traces=[]
for entry in IiIuG: traces.append(entry[0])
for entry in AuSIS: traces.append(entry[0])

stream=Stream(traces=traces)
stream.plot(outfile=OUT_PATH+'broadband'+EVENT_ID+'.pdf',equal_scale=False, automerge=False)

streamLow=Stream(traces=traces)
streamLow.filter('lowpass',freq=0.05)
stream.plot(outfile=OUT_PATH+'lowpass_0.05_'+EVENT_ID+'.pdf',equal_scale=False, automerge=False)

streamHigh=Stream(traces=traces)
streamHigh.filter('highpass',freq=4.0)
stream.plot(outfile=OUT_PATH+'highpass_4.0_'+EVENT_ID+'.pdf',equal_scale=False, automerge=False)

streamMid=Stream(traces=traces)
streamMid.filter('bandpass',freqmin=0.05,freqmax=0.15)
stream.plot(outfile=OUT_PATH+'bandpass_0.05,0.15_'+EVENT_ID+'.pdf',equal_scale=False, automerge=False)

print('\n\nRuntime: '+str(UTCDateTime()-start)+' seconds\n\n')
os.system('say "Waveform plotting complete"')
#################################################################