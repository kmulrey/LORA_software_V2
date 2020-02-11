import numpy as np
import re
from datetime import datetime


eventlist=np.genfromtxt(open('eventlist_xmaxfit_done_701.txt','r'))
LORA_directory='/vol/astro7/lofar/lora/LORA_V2/'
event_id_offset=1262304000


for e in np.arange(len(eventlist)):
    event=eventlist[e]
    try:
        timestamp=event+event_id_offset # event timestamp in UTC- used to find LORA file
        dt_object = datetime.utcfromtimestamp(timestamp)
    
        year=dt_object.year
        month=dt_object.month
        day=dt_object.day
        hour=dt_object.hour
        minute=dt_object.minute
        sec=dt_object.second
        LORA_file=LORA_directory+'LORAdata-'+str(dt_object.year)+str(dt_object.month).zfill(2)+str(dt_object.day).zfill(2)+'T'+str(dt_object.hour).zfill(2)+str(dt_object.minute).zfill(2)+str(dt_object.second).zfill(2)+'.dat'



        f=open(LORA_file,'r')
        lines=f.readlines()
        header_info=lines[1].strip().split()
        f.close()
        #print 'reading LORA file: ',LORA_file

    except:
        #print 'can\'t read event {0}'.format(int(event))
        print int(event), LORA_file
