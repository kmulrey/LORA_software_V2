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
from pytz import timezone
import ROOT


data_dir='/Users/kmulrey/LOFAR/LORA/LORAraw/'
LORA4times=np.genfromtxt(data_dir+'LORAtime4')
outputdir='LORAnew/'

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
        #print 'getting timestamp for----> {0}  {1}'.format(d, int((detectors[d].number-1)/4))
        event.get_event_timestamp(detectors[d],lasas[int((detectors[d].number-1)/4)])

    for l in np.arange(LORA.nLasaA):
        event.cal_event_timestamp(detectors,lasas[l])
   

    event.do_arrival_time_diff(detectors)
    event.do_arrival_direction(detectors,ev)

    event.do_COM_core(detectors,ev)
    event.find_density(detectors,ev)

    event.fit_arrival_direction(detectors,ev)
#final event time
# arrival time diff
#  arrival direction

    #make output file

    '''
    dt_object = datetime.utcfromtimestamp(timestamp)

    year=dt_object.year
    month=dt_object.month
    day=dt_object.day
    hour=dt_object.hour
    minute=dt_object.minute
    sec=dt_object.second



    output_file=open(outputdir+'LORAdata-'+str(dt_object.year)+str(dt_object.month).zfill(2)+str(dt_object.day).zfill(2)+'T'+str(dt_object.hour).zfill(2)+str(dt_object.minute).zfill(2)+str(dt_object.second).zfill(2)+'.dat','w')
    output_file.write('//UTC_Time(secs)\tnsecs\t\tCore(X)\t\tCore(Y)\t\tElevation\tAzimuth\t\tEnergy(eV)\tCoreE(X)\tCoreE(Y)\tMoliere_rad(m)\t\tElevaErr\tAziErr\tEnergyErr(eV)\tNe\t\tNeErr\t\tCorCoef_XY\t\tNe_RefA\t\tNeErr_RefA\t\tEnergy_RefA\t\tEnergyErr_RefA\n')


    output_file.close()
    '''

    #plot.plot_traces(detectors)

#first_hit=process.find_event_time(info)
#first_hit=process.find_event_time(info)
#plot.plot_traces(info)
#except:
#print 'no data'

    

