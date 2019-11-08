import numpy as np
from datetime import datetime
import os
import read_functions as read
import ROOT
import plot_functions as plot

data_dir='/Users/kmulrey/LOFAR/LORA/LORAraw/'
LORA4times=np.genfromtxt(data_dir+'LORAtime4')

print len(LORA4times)

for t in np.arange(50,51):
    timestamp = int(LORA4times[t][0])
    ns_timestamp = int(LORA4times[t][1])
    tag=read.find_tag(timestamp,ns_timestamp,data_dir)
    print '{0}   {1}:  {2}'.format(timestamp,ns_timestamp,tag)
    test=read.return_root(tag,timestamp,ns_timestamp,data_dir)



#print test[0].keys()
#plot.plot_traces(test)



    

