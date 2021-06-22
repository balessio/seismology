"""
Before running the program, enter your specifications into
SETUP_NEtoTR.txt
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
with open('SETUP/SETUP_NEtoTR.txt', 'r') as f:
    setup=f.read().splitlines()
PATH=setup[0][5:]
EVENT_ID=setup[1][9:]

client=Client('ISC')
event=client.get_events(eventid=EVENT_ID)[0]
lonE=event.preferred_origin().longitude
latE=event.preferred_origin().latitude
print(str(lonE)+','+str(latE))

client=Client('IRIS')
failures=0
traceAndBazm=[]
for mseed in os.listdir(PATH): # reads all mseed files in path
    try:
        if '.mseed' in mseed:
            st = read(PATH+mseed)
            for trace in st:
                traceInfo=trace.id.split('.')
                if traceInfo[3][2] != 'Z':
                    station=client.get_stations(network=traceInfo[0],
                                                station=traceInfo[1],level='response')
                    coords=station.get_coordinates(seed_id=trace.id)
                    print(separation(latE,lonE,coords['latitude'],coords['longitude']))
                    traceAndBazm.append((trace,separation(latE,lonE,coords['latitude'],
                                        coords['longitude'])[2]))
    except: failures+=1
print('\nFailed '+str(failures)+' time(s)\n')
print(traceAndBazm)

traceIDdict={} ## example: {'IU.KIP.00':(stream,back azimuth)}
for entry in traceAndBazm: 
    traceInfo=entry[0].id.split('.')
    try: traceIDdict[traceInfo[0]+traceInfo[1]+traceInfo[2]][0].append(entry[0])
    except: traceIDdict[traceInfo[0]+traceInfo[1]+traceInfo[2]]=(Stream(entry[0]),entry[1])
print(traceIDdict)
for key in traceIDdict:
    try:
        stream=traceIDdict[key][0].rotate(method='NE->RT',back_azimuth=traceIDdict[key][1])
        for trace in stream: trace.write(PATH+trace.id+'.mseed', format='MSEED')
    except: 
        brokenTraces=''
        for trace in traceIDdict[key][0]: brokenTraces+=trace.id+','
        print('\nThe data in '+brokenTraces+' was too cruddy to perform the rotation\n')

print('\n\nRuntime: '+str(UTCDateTime()-start)+' seconds\n\n')
os.system('say "Coordinate system rotation complete"')
#################################################################