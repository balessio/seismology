"""
Before running the program, enter your specifications into
SETUP_waveforms.txt
"""
###############################################################
import matplotlib.pyplot as plt
from obspy import read
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
import os

start = UTCDateTime() # time the program

#### Read program specs from the setup file
with open('SETUP/SETUP_separateWaveformPlot.txt', 'r') as f:
    setup=f.read().splitlines()
IN_PATH=setup[0][8:]
OUT_PATH=setup[1][9:]

failures=0
with open(OUT_PATH+'info.txt', 'w') as textFile:
    for mseed in os.listdir(IN_PATH): # reads all mseed files in path
        try:
            stream = read(IN_PATH+mseed)
            textFile.write(stream.__str__(extended=True)+'\n')
            for trace in stream: stream.plot(outfile=OUT_PATH+trace.id+'.pdf')
        except: failures+=1
print ('\nFailed '+str(failures)+' time(s)\n')
print('\n\nRuntime: '+str(UTCDateTime()-start)+' seconds\n\n')
os.system('say "Waveform plotting complete"')
#################################################################