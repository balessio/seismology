"""
Before running the program, enter your specifications into
SETUP_fetchEventData.txt
"""
######################################################################
import obspy
from obspy.clients.fdsn.client import Client
from obspy.core import UTCDateTime
from obspy.core.event import Event
from obspy.core.event import FocalMechanism
from obspy.core.event import Origin
import os

start=UTCDateTime() # Time the program

#### Read program specs from the setup file
with open('SETUP/SETUP_fetchEventData.txt', 'r') as f:
    setup=f.read().splitlines()
START_TIME=UTCDateTime(setup[0][29:])
END_TIME=UTCDateTime(setup[1][27:])
MIN_LATITUDE=float(setup[2][13:])
MAX_LATITUDE=float(setup[3][13:])
MIN_LONGITUDE=float(setup[4][14:])
MAX_LONGITUDE=float(setup[5][14:])
MIN_MAGNITUDE=float(setup[6][14:])
MAX_MAGNITUDE=float(setup[7][14:])
QUAKEML_FILE=setup[8][13:]
TEXT_FILE=setup[9][10:]

#### Download event data from ISC
client = Client('ISC')
events = client.get_events(starttime=START_TIME, endtime=END_TIME, 
    minlatitude=MIN_LATITUDE, maxlatitude=MAX_LATITUDE, 
    minlongitude=MIN_LONGITUDE, maxlongitude=MAX_LONGITUDE, 
    includeallmagnitudes=True, includeallorigins=True, 
    minmagnitude=MIN_MAGNITUDE, maxmagnitude=MAX_MAGNITUDE)
#### Check for usable data
print('\n...'+str(len(events))+' events found...\n')
cleanEvents=client.get_events(eventid='609301') # just to create an
cleanEvents.clear() # empty catalog...
numBrokenEvents=0
for event in events:
    a=0.0
    b=1.0
    try:
        a=float(event.preferred_origin().longitude)
        a=float(event.preferred_origin().latitude)
        a=float(event.magnitudes[0].mag)
        a=event.preferred_origin()
        b/=float(a.depth)
        cleanEvents.append(event)
    except: numBrokenEvents+=1
print('\n...Removed '+str(numBrokenEvents)+' broken events...\n')
try: print('\nFirst event: '+str(cleanEvents[0].resource_id))
except: raise Exception('\n\nThere are no useful events!\n\n')

#### Write it into the file
cleanEvents.write(QUAKEML_FILE, format="QUAKEML") 
with open(TEXT_FILE, 'w') as textFile:
    textFile.write(cleanEvents.__str__(print_all=True))

print('\n\nRuntime: '+str(UTCDateTime()-start)+' seconds\n\n')
os.system('say "Events collected"')
######################################################################
