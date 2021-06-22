"""
Before running the program, enter your specifications into
.csv files within SETUP/fetchWaveforms
with the following format:
blank,directory (with optional filename),
skip row
Network,Station,Location,Channel,Start time,End time
as many extra rows as wanted

end on '/' before 'xxx.mseed' if you want separate files for each trace, just include 
the folder within inputs/MSEED, to be used for separateWaveformPlot.py instead of 
bigWaveformPlot.py
"""
###############################################################
import csv
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
import os

start = UTCDateTime() # time the program

for SETUP_FILE in os.listdir('SETUP/fetchWaveforms/'):
    #### Read the program specs from the setup file
    bulk=[]
    client = Client("IRIS")
    try:
        with open('SETUP/fetchWaveforms/'+SETUP_FILE, 'r') as f:
            setup = csv.reader(f, delimiter=',')
            rowNavigator=0
            for row in setup:
                if (rowNavigator==0):
                    MSEED_FILE=row[1]
                elif (rowNavigator>1):
                    bulk.append(((row[0],row[1],row[2],row[3],UTCDateTime(row[4]),UTCDateTime(row[5]))))
                rowNavigator+=1
        if (len(bulk)<1): raise Exception('\n\nYou must specify what you want!\n\n')

        #### Download the stream from IRIS
        stream = client.get_waveforms_bulk(bulk)

        #### If the MSEED filename is left blank in the setup file, then create separate files for each 
        #### trace. Otherwise create one file to hold the entire stream.
        if '.mseed' not in MSEED_FILE:
            for trace in stream:
                trace.write('inputs/MSEED/'+MSEED_FILE+trace.id+'.mseed', format='MSEED')
        else: stream.write('inputs/MSEED/'+MSEED_FILE, format='MSEED')
    except: print('\nSkipped '+SETUP_FILE)

print('\n\nRuntime: '+str(UTCDateTime()-start)+' seconds\n\n')
os.system('say "Waveforms collected"')
#################################################################


## attach and remove response
## spectral density