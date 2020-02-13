import numpy as np
from datetime import datetime
import os
import read_functions as read
import ROOT
#import plot_functions as plot
import LORAparameters as LORA
import process_functions as process
import detector as det
import event as event
from pytz import timezone
import ROOT
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-e', '--event',type='int',help='event numbert',default=136950320)
parser.add_option('-i', '--ind',type='int',help='line number of event',default=-1)

(options, args) = parser.parse_args()
e=int(options.event)
ind=int(options.ind)





#data_dir='/Users/kmulrey/LOFAR/LORA/LORAraw/'
data_dir='/vol/astro3/lofar/vhecr/lora_triggered/LORAraw/'
LORA4times=np.genfromtxt('/vol/astro3/lofar/vhecr/lora_triggered/LORA/'+'LORAtime4')
#outputdir='LORAnew/'
outputdir='/vol/astro7/lofar/lora/LORA_V2/'

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

#event.read_attenuation()

eventlist=np.genfromtxt(open('eventlist_xmaxfit_done_701.txt','r'))

if ind>-1:
    LOFAR_id=str(eventlist[ind])

else:   
    LOFAR_id=str(e)

timestamp=int(float(LOFAR_id))+LORA.event_id_offset
#LOFAR_id=str(int(timestamp-LORA.event_id_offset))

#print LORA4times.T[0]
index=-1

for i in np.arange(len(LORA4times)):
    #print int(LORA4times[i][0]), timestamp
    if int(LORA4times[i][0])== timestamp:
        index=i
        continue
print 'index: ',index

if i>-1:
    run=1
else:
    run=0
    print 'can\'t  find timestamp'
for t in np.arange(run):


    ns_timestamp = int(LORA4times[index][1])

    #for i in np.arange(1):
    #for t in np.arange(30,31):
    #for t in np.arange(30,30):

    #timestamp = int(LORA4times[t][0])
    #ns_timestamp = int(LORA4times[t][1])
    
    
    #LOFAR_id=str(int(timestamp-LORA.event_id_offset))
    ev=event.Event(LOFAR_id,'V1')
    
    tag=read.find_tag(timestamp,ns_timestamp,data_dir)
    print '_____________________'

    print '{2}:    {0}   {1}'.format(timestamp,ns_timestamp,tag)
    
    log_file_info=read.log_file(tag,data_dir)
    info=read.return_root(tag,timestamp,ns_timestamp,data_dir)
    sec_info0,sec_info1,sec_info2=read.return_second_data(tag,timestamp,ns_timestamp,data_dir)
    log_info=read.return_log_data(tag,timestamp,ns_timestamp,data_dir)
    noise_info=read.return_noise_data(tag,timestamp,ns_timestamp,data_dir)

    LOFAR_trig=log_file_info['LOFAR_trig']
    lasa1_status=log_file_info['lasa1_status']
    lasa2_status=log_file_info['lasa2_status']
    lasa3_status=log_file_info['lasa3_status']
    lasa4_status=log_file_info['lasa4_status']
    lasa5_status=log_file_info['lasa5_status']

    print LOFAR_trig,lasa1_status,lasa2_status,lasa3_status,lasa4_status,lasa5_status


    if (sec_info0[0]['GPS_time_stamp']-info[0]['gps'])>1:
        print 'Can\'t find corresponding 1 Sec message: {0}  {1}'.format(int(info[0]['gps']),int(sec_info0[0]['GPS_time_stamp']))
        for l in np.arange(5):
           lasas[l].sec_flag=1
        continue
    
    det.load_event_information(info,detectors)
    det.load_sec_information(sec_info0,sec_info1,sec_info2,lasas)
    det.load_log_information(log_info,detectors)
    det.load_noise_information(noise_info,detectors)


    for d in np.arange(LORA.nDetA):
        event.find_counts(detectors[d])
        event.get_arrival_time(detectors[d])
        #print 'getting timestamp for----> {0}  {1}'.format(d, int((detectors[d].number-1)/4))
        event.get_event_timestamp(detectors[d],lasas[int((detectors[d].number-1)/4)])
        
        
        print 'threshold: {0:.2f},  corr peak: {1:.2f}, baseline: {2:.2f}, mean: {6:.2f}, sigma: {7:.2f}, peak: {3:.2f},  trigger: {4}, {5}'.format(detectors[d].threshold,detectors[d].corrected_peak,detectors[d].trace_mean, detectors[d].peak,detectors[d].trig,detectors[d].trigg_condition,detectors[d].sec_mean,detectors[d].sec_sigma)

    for l in np.arange(LORA.nLasaA):
        event.cal_event_timestamp(detectors,lasas[l])
   
   
    event.do_arrival_time_diff(detectors)
    event.do_arrival_direction(detectors,ev)

    event.do_COM_core(detectors,ev)
    event.find_density(detectors,ev)

    event.fit_arrival_direction(detectors,ev)
    event.fit_NKG(detectors,ev)


    

    dt_object = datetime.utcfromtimestamp(timestamp)

    year=dt_object.year
    month=dt_object.month
    day=dt_object.day
    hour=dt_object.hour
    minute=dt_object.minute
    sec=dt_object.second



    output_file=open(outputdir+'LORAdata-'+str(dt_object.year)+str(dt_object.month).zfill(2)+str(dt_object.day).zfill(2)+'T'+str(dt_object.hour).zfill(2)+str(dt_object.minute).zfill(2)+str(dt_object.second).zfill(2)+'.dat','w')
    output_file.write('//UTC_Time(secs)\tnsecs\t\tCore(X)\t\tCore(Y)\t\tElevation\tAzimuth\t\tEnergy(eV)\tCoreE(X)\tCoreE(Y)\tMoliere_rad(m)\t\tElevaErr\tAziErr\tEnergyErr(eV)\tNe\t\tNeErr\t\tCorCoef_XY\t\tNe_RefA\t\tNeErr_RefA\t\tEnergy_RefA\t\tEnergyErr_RefA\t\tLOFAR_Trigger\t\tLASA1_trig\t\tLASA2_trig\t\tLASA3_trig\t\tLASA4_trig\t\tLASA5_trig\t\tLASA6_trig\t\tLASA7_trig\t\tLASA8_trig\t\tLASA9_trig\t\tLASA10_trig\t\tLASA1_status\t\tLASA2_status\t\tLASA3_status\t\tLASA4_status\t\tLASA5_status\t\tLASA6_status\t\tLASA7_status\t\tLASA8_status\t\tLASA9_status\t\tLASA10_status\n')
    output_file.write('{0}\t{1}\t{2:.2f}\t{3:.2f}\t{4:.2f}\t{5:.2f}\t{6}\t{7:.2f}\t{8:.2f}\t{9:.2f}\t{10:.2f}\t{11:.2f}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\t{27},\t{28}\t{29}\t{30}\t{31}\t{32}\t{33}\t{34}\t{35}\t{36}\t{37}\t{38}\t{39}\t{40}\n'.format(int(ev.UTC_min), int(ev.nsec_min/10), ev.x_core, ev.y_core, ev.fit_elevation, ev.fit_phi, ev.energy*np.power(10,15), ev.x_core_err, ev.y_core_err, ev.Rm, ev.fit_elevation_err, ev.fit_phi_err, ev.energy_err*np.power(10,15), ev.Ne, ev.Ne_err, ev.CorCoef_xy, ev.Ne_RefA, ev.NeErr_RefA, ev.Energy_RefA*np.power(10,15), ev.EnergyErr_RefA*np.power(10,15), LOFAR_trig, int(detectors[0].trigg_condition),int(detectors[4].trigg_condition),int(detectors[8].trigg_condition),int(detectors[12].trigg_condition),int(detectors[16].trigg_condition),0,0,0,0,0, lasa1_status,lasa2_status,lasa3_status,lasa4_status,lasa5_status,0,0,0,0,0))

    output_file.write('//Detector coordinates w.r.t CS002 LBA center\n')
    output_file.write('//Det_no.    X_Cord(m)    Y_Cord(m)    Z_Cord(m)    UTC_time    (10*nsecs)  int_ADC_count  gain(ADC/muon)    Particle_Density(/m2), baseline(ADC), rms(ADC), corrected_peak(ADC), corrected_threshold(ADC), avg_mean(ADC), avg_sigma(ADC)\n')

    for i in np.arange(len(detectors)):
        print i,detectors[i].x_cord,detectors[i].y_cord,detectors[i].z_cord,detectors[i].gps,detectors[i].cal_time,detectors[i].density
        output_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.4f}\t{7}\t{8:.4f}\t{9:.2f}\t{10:.2f}\t{11:.2f}\t{12:.2f}\t{13:.2f}\t{14:.2f}\t{15}\n'.format(i+1,detectors[i].x_cord,detectors[i].y_cord,detectors[i].z_cord,int(detectors[i].gps),int(detectors[i].cal_time),detectors[i].trace_int_counts,detectors[i].gain,detectors[i].density,detectors[i].trace_mean,detectors[i].trace_rms,detectors[i].corrected_peak,detectors[i].corrected_threshold,detectors[i].sec_mean,detectors[i].sec_sigma, detectors[i].trig))
    
          
    output_file.close()
    


#plot.plot_traces(detectors)

#first_hit=process.find_event_time(info)
#first_hit=process.find_event_time(info)
#plot.plot_traces(info)
#except:
#print 'no data'


