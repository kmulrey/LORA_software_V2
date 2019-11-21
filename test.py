import numpy as np
from datetime import datetime
import os
import read_functions as read
import ROOT
import plot_functions as plot
import LORAparameters as LORA
import process_functions as process
import detector as det
import event as event

data_dir='/Users/kmulrey/LOFAR/LORA/LORAraw/'
LORA4times=np.genfromtxt(data_dir+'LORAtime4')


detectors=[]
lasas=[]

for i in np.arange(LORA.nDetA):
    detectors.append(det.Detector('Det'+(str(i+1))))

lasas=[]
for i in np.arange(LORA.nLasaA):
    lasas.append(det.Lasa('Lasa'+(str(i+1))))
                                  
det.load_positions(detectors)
det.load_signal(detectors)
det.load_gain(detectors)



for t in np.arange(30,31):
    timestamp = int(LORA4times[t][0])
    ns_timestamp = int(LORA4times[t][1])
    
    
    LOFAR_id=str(int(timestamp-LORA.event_id_offset))
    ev=event.Event(LOFAR_id,'V1')
    
    tag=read.find_tag(timestamp,ns_timestamp,data_dir)
    print '_____________________'

    print '{3}:    {0}   {1}:  {2}'.format(timestamp,ns_timestamp,tag,t)
        #try:
    info=read.return_root(tag,timestamp,ns_timestamp,data_dir)
    sec_info0,sec_info1,sec_info2=read.return_second_data(tag,timestamp,ns_timestamp,data_dir)
    
    '''
    print '................'
    print int(sec_info0[0]['GPS_time_stamp']),sec_info0[0]['sync'],sec_info0[0]['CTP'],sec_info0[0]['quant']
    print int(sec_info1[0]['GPS_time_stamp']),sec_info1[0]['sync'],sec_info1[0]['CTP'],sec_info1[0]['quant']
    print int(sec_info2[0]['GPS_time_stamp']),sec_info2[0]['sync'],sec_info2[0]['CTP'],sec_info2[0]['quant']
    print '_____________________'
    '''

    

    if (sec_info0[0]['GPS_time_stamp']-info[0]['gps'])>1:
        print 'Can\'t find corresponding 1 Sec message: {0}  {1}'.format(int(info[0]['gps']),int(sec_info0[0]['GPS_time_stamp']))
        for l in np.arange(5):
           lasas[l].sec_flag=1

        continue
    
    det.load_event_information(info,detectors)
    det.load_sec_information(sec_info0,sec_info1,sec_info2,lasas)


    for d in np.arange(LORA.nDetA):
        event.find_counts(detectors[d])
        event.get_arrival_time(detectors[d])
        print 'getting timestamp for----> {0}  {1}'.format(d, int((detectors[d].number-1)/4))
        event.get_event_timestamp(detectors[d],lasas[int((detectors[d].number-1)/4)])

    event.cal_event_timestamp(detectors,lasas[0])
    event.cal_event_timestamp(detectors,lasas[1])
    event.cal_event_timestamp(detectors,lasas[2])
    event.cal_event_timestamp(detectors,lasas[3])
    event.cal_event_timestamp(detectors,lasas[4])

    #plot.plot_traces(detectors)

#first_hit=process.find_event_time(info)
#first_hit=process.find_event_time(info)
#plot.plot_traces(info)
#except:
#print 'no data'

    

