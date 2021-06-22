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
from obspy.core.stream import Stream
from obspy.core.trace import Trace

start = UTCDateTime() # time the program

#### Read program specs from the setup file
with open('SETUP/SETUP_indStnPlt.txt', 'r') as f:
    setup=f.read().splitlines()
IN_PATH=setup[0][8:]
OUT_PATH=setup[1][9:]

client=Client('IRIS')
failures=0
traces=[]
with open(OUT_PATH+'info.txt', 'w') as textFile:
    for mseed in os.listdir(IN_PATH): # reads all mseed files in path
        try:
            if '.mseed' in mseed:
                stream = read(IN_PATH+mseed)
                for trace in stream: traces.append(trace)
        except: failures+=1
print('\nFailed '+str(failures)+' time(s)\n')

traceIDdict={}
for trace in traces:
    idList=trace.id.split('.')
    try: traceIDdict[idList[0]+idList[1]+idList[2]].append(trace)
    except: traceIDdict[idList[0]+idList[1]+idList[2]]=[trace]

for key in traceIDdict:
    traces=traceIDdict[key]
    idList=traces[0].id.split('.')
    stationName=idList[0]+'.'+idList[1]+'.'+idList[2]

    stream=Stream(traces=traces)
    stream.plot(outfile=OUT_PATH+stationName+'broadband'+'.pdf')

    streamLow=Stream(traces=traces)
    streamLow.filter('lowpass',freq=0.05)
    stream.plot(outfile=OUT_PATH+stationName+'lowpass_0.05_'+'.pdf')

    streamHigh=Stream(traces=traces)
    streamHigh.filter('highpass',freq=4.0)
    stream.plot(outfile=OUT_PATH+stationName+'highpass_4.0_'+'.pdf')

    streamMid=Stream(traces=traces)
    streamMid.filter('bandpass',freqmin=0.05,freqmax=0.15)
    stream.plot(outfile=OUT_PATH+stationName+'bandpass_0.05,0.15_'+'.pdf')

print('\n\nRuntime: '+str(UTCDateTime()-start)+' seconds\n\n')
os.system('say "Waveform plotting complete"')
#################################################################