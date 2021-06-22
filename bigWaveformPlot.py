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
with open('SETUP/SETUP_bigWaveformPlot.txt', 'r') as f:
    setup=f.read().splitlines()
MSEED_FILE=setup[0][11:]
TEXT_FILE=setup[1][10:]
IMAGE_FILE=setup[2][11:]

#### read in stream from file
stream = read(MSEED_FILE)
print(stream[0].stats.location)
with open(TEXT_FILE, 'w') as textFile:
    textFile.write(stream.__str__(extended=True))
stream.plot(outfile=IMAGE_FILE)
#stream.plot(type='section',outfile=IMAGE_FILE)

print('\n\nRuntime: '+str(UTCDateTime()-start)+' seconds\n\n')
os.system('say "Waveform plotting complete"')
#################################################################