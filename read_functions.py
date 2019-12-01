import numpy as np
from datetime import datetime
import os
import ROOT
import LORAparameters as LORA

nTraceV1=4000
nDetV1=20
dtV1=2.5
timesV1=np.arange(0,4000*dtV1,dtV1)


def find_tag(timestamp,ns_timestamp,data_dir):
    #print '\n________________________\n'
    #timestamp = LORA4times[t][0]
    #ns_timestamp = LORA4times[t][1]
    
    dt_object = datetime.fromtimestamp(timestamp)
    year=dt_object.year
    month=dt_object.month
    day=dt_object.day
    hour=dt_object.hour
    minute=dt_object.minute
    second=dt_object.second
    
    #print year,month,day,hour
    #print int(timestamp),int(ns_timestamp)
    
    filestr1= str(year)+str(month).zfill(2)+str(day-1).zfill(2)
    filestr2= str(year)+str(month).zfill(2)+str(day).zfill(2)
    #filestr3= str(year)+str(month).zfill(2)+str(day-1).zfill(2)
    
    filelist=[]
    
    for file in os.listdir(data_dir):
        if filestr1 in file or filestr2 in file:
            if '.log' in file:
                #print(os.path.join(data_dir, file))
                filelist.append(os.path.join(data_dir, file))


    found=0
    filetag='no'
    
    for i in np.arange(len(filelist)):
        #print filelist[i]
        #for i in np.arange(1):
        file=open(filelist[i])
        stamps=np.genfromtxt(file,skip_header=10,usecols=(2,3))
        file.close()
        if len(stamps>0):
            if len(stamps[(stamps.T[0]==timestamp)*(stamps.T[1]==ns_timestamp)])>0:
                
                filetag=filelist[i].split('raw/')[1].split('.')[0]
                found=1

    if found==1:
        return filetag
    else:
        return 'no_match'

def getTime(det, entry):
    det.GetEntry(entry)
    ymd=det.GetLeaf('YMD').GetValue()
    gps=det.GetLeaf('GPS_time_stamp').GetValue()
    ctd=det.GetLeaf('CTD').GetValue()
    nsec=det.GetLeaf('nsec').GetValue()
    return ymd,gps,ctd,nsec

def getTimeSec(lasa, entry):
    lasa.GetEntry(entry)
    gps=lasa.GetLeaf('GPS_time_stamp').GetValue()
    return gps

def find_entry_number(lora_utc,lora_nsec,tree_event):
    event=-1
    
    det1=tree_event.GetBranch('Det1')
    det5=tree_event.GetBranch('Det5')
    det9=tree_event.GetBranch('Det9')
    det13=tree_event.GetBranch('Det13')
    det17=tree_event.GetBranch('Det17')
    
    det1.GetLeaf('GPS_time_stamp')
    nE= det1.GetEntries()
    #times=np.zeros([nLasa,nE])
    trigger_check=0
    diff_best=1e10
    for e in np.arange(nE):
        
        ymd1,gps1,ctd1,nsec1= getTime(det1,e)
        ymd2,gps2,ctd2,nsec2= getTime(det5,e)
        ymd3,gps3,ctd3,nsec3= getTime(det9,e)
        ymd4,gps4,ctd4,nsec4= getTime(det13,e)
        ymd5,gps5,ctd5,nsec5= getTime(det17,e)
        times=[gps1,gps2,gps3,gps4,gps5]
        times_ns=[nsec1,nsec2,nsec3,nsec4,nsec5]
        
        if lora_utc in times:
            diff=np.max(np.abs(lora_nsec-np.asarray(times_ns)[np.asarray(times_ns)>1]))
            if diff<10000:  # w/in 10 us to avoid mis-triggers
                if diff<diff_best:
                    diff_best=diff
                    trigger_check=1
                    event=e



    return event


def getDataV1(det, entry):
    
    det.GetEntry(entry)
    
    detector=det.GetLeaf('detector').GetValue()
    ymd=det.GetLeaf('YMD').GetValue()
    gps=det.GetLeaf('GPS_time_stamp').GetValue()
    ctd=det.GetLeaf('CTD').GetValue()
    nsec=det.GetLeaf('nsec').GetValue()
    trigg_condition=det.GetLeaf('Trigg_condition').GetValue()
    try:
        trigg_pattern=det.GetLeaf('Trigg_pattern').GetValue()
    except:
        trigg_pattern=-1
    total_counts=det.GetLeaf('Total_counts').GetValue()
    pulse_height=det.GetLeaf('Pulse_height').GetValue()
    pulse_width=det.GetLeaf('Pulse_width').GetValue()
    counts=det.GetLeaf('counts')
    hold=np.zeros([nTraceV1])
    for i in np.arange(nTraceV1):
        hold[i]=counts.GetValue(i)

    info={'det':detector,'ymd':ymd,'gps':gps,'ctd':ctd,'nsec':nsec,'trigg_condition':trigg_condition,'trigg_pattern':trigg_pattern,'total_counts':total_counts,'pulse_height':pulse_height,'pulse_width':pulse_width,'counts':hold}
    return info



def return_root(filename,utc,nsec,data_dir):

    log_file=open(data_dir+filename+'.log','r')
    root_file=ROOT.TFile.Open(data_dir+filename+'.root')
    
    tree_sec = root_file.Get("Tree_sec")
    tree_event = root_file.Get("Tree_event")
    tree_log = root_file.Get("Tree_log")
    tree_noise = root_file.Get("Tree_noise")

    event_index=find_entry_number(utc,nsec,tree_event)
    all_info=[]
    for i in np.arange(nDetV1):
        detname='Det'+str(1+i)
        det=tree_event.GetBranch(detname)
        info=getDataV1(det,event_index)
        all_info.append(info)

    return all_info















def find_sec_number(lora_utc,tree_sec,i):
    event=-1
    
    lasa1=tree_sec.GetBranch('Lasa'+str(i+1))


    lasa1.GetLeaf('GPS_time_stamp')
    nE= lasa1.GetEntries()
    #times=np.zeros([nLasa,nE])
    diff_best=1e10
    event=-1
    index_found=0
    
    for e in np.arange(nE):
        
        gps1= getTimeSec(lasa1,e)
        if gps1==lora_utc and index_found==0:
            event=e
            index_found=1
        if gps1>lora_utc and index_found==0:
            event=e
            index_found=1

    return event













def getSecV1(det, entry):
    
    det.GetEntry(entry)
    
    lasa=det.GetLeaf('Lasa').GetValue()
    
    YMD=det.GetLeaf('YMD').GetValue()

    GPS_time_stamp=det.GetLeaf('GPS_time_stamp').GetValue()
    sync=det.GetLeaf('sync').GetValue()
    CTP=det.GetLeaf('CTP').GetValue()
    quant=det.GetLeaf('quant').GetValue()
    Channel_1_Thres_count_high=det.GetLeaf('Channel_1_Thres_count_high').GetValue()
    Channel_1_Thres_count_low=det.GetLeaf('Channel_1_Thres_count_low').GetValue()
    Channel_2_Thres_count_high=det.GetLeaf('Channel_2_Thres_count_high').GetValue()
    Channel_2_Thres_count_low=det.GetLeaf('Channel_2_Thres_count_low').GetValue()
    Satellite_info=det.GetLeaf('Satellite_info').GetValue()

    info={'lasa':lasa,'YMD':YMD,'GPS_time_stamp':GPS_time_stamp,'sync':sync,'CTP':CTP,'quant':quant}
    return info


def return_second_data(filename,utc,nsec,data_dir):

    log_file=open(data_dir+filename+'.log','r')
    root_file=ROOT.TFile.Open(data_dir+filename+'.root')
    tree_sec = root_file.Get("Tree_sec")
    tree_event = root_file.Get("Tree_event")
    tree_log = root_file.Get("Tree_log")
    tree_noise = root_file.Get("Tree_noise")

    entry=np.zeros([LORA.nLASA])
    
    for i in np.arange(LORA.nLASA):
        entry[i]=find_sec_number(utc,tree_sec,i)
    
    all_info=[]
    all_info1=[]
    all_info2=[]

    
    for i in np.arange(LORA.nLASA):
        lasaname='Lasa'+str(1+i)
        det=tree_sec.GetBranch(lasaname)
        all_info.append(getSecV1(det, int(entry[i])))
        all_info1.append(getSecV1(det, int(entry[i]+1)))
        all_info2.append(getSecV1(det, int(entry[i]+2)))

    return all_info,all_info1,all_info2


def getLogV1(det, entry):
    
    det.GetEntry(entry)
    
    
    YMD=det.GetLeaf('YMD').GetValue()
    GPS_time_stamp=det.GetLeaf('Time_stamp').GetValue()
    Threshold_low=det.GetLeaf('Channel_thres_low').GetValue()
    info={'threshold':Threshold_low}
    return info

def getNoiseV1(det, entry):
    
    det.GetEntry(entry)
    
    sigma=det.GetLeaf('Sigma').GetValue()
    mean=det.GetLeaf('Mean').GetValue()
    info={'mean':mean,'sigma':sigma}
    return info

def find_noise_number(lora_utc,tree_noise,d):
    
    event=-1
    
    det1=tree_noise.GetBranch('Det'+str(d))
    
    det1.GetLeaf('GPS_time_stamp')
    nE= det1.GetEntries()
    #times=np.zeros([nLasa,nE])
    
    event=-1
    index_found=0
    
    for e in np.arange(nE):
        
        det1.GetEntry(e)
        gps1=det1.GetLeaf('GPS_time_stamp').GetValue()
        if gps1>lora_utc and index_found==0:
            event=e
            index_found=1

    return event

def find_log_number(lora_utc,tree_log,d):
    
    event=-1
    
    det1=tree_log.GetBranch('Det'+str(d))
    
    det1.GetLeaf('Time_stamp')
    nE= det1.GetEntries()
    #times=np.zeros([nLasa,nE])
  
    event=-1
    index_found=0
    
    for e in np.arange(nE):
        
        det1.GetEntry(e)
        gps1=det1.GetLeaf('Time_stamp').GetValue()
        if gps1>lora_utc and index_found==0:
            event=e
            index_found=1

    return event



def return_log_data(filename,utc,nsec,data_dir):
    
    log_file=open(data_dir+filename+'.log','r')
    root_file=ROOT.TFile.Open(data_dir+filename+'.root')
    tree_sec = root_file.Get("Tree_sec")
    tree_event = root_file.Get("Tree_event")
    tree_log = root_file.Get("Tree_log")
    tree_noise = root_file.Get("Tree_noise")
    
    entry1=find_log_number(utc,tree_log,1)
    entry2=find_log_number(utc,tree_log,5)
    entry3=find_log_number(utc,tree_log,9)
    entry4=find_log_number(utc,tree_log,13)
    entry5=find_log_number(utc,tree_log,17)
    
    all_info=[]
    
    
    for i in np.arange(LORA.nLORA):
        detname='Det'+str(1+i)
        det=tree_log.GetBranch(detname)
        if i==0 or i==1 or i==2 or i==3:
            all_info.append(getLogV1(det, entry1))
        if i==4 or i==5 or i==6 or i==7:
            all_info.append(getLogV1(det, entry2))
        if i==8 or i==9 or i==10 or i==11:
            all_info.append(getLogV1(det, entry3))
        if i==12 or i==13 or i==14 or i==15:
            all_info.append(getLogV1(det, entry4))
        if i==16 or i==17 or i==18 or i==19:
            all_info.append(getLogV1(det, entry5))
    
    return all_info






def return_noise_data(filename,utc,nsec,data_dir):
    
    log_file=open(data_dir+filename+'.log','r')
    root_file=ROOT.TFile.Open(data_dir+filename+'.root')
    tree_sec = root_file.Get("Tree_sec")
    tree_event = root_file.Get("Tree_event")
    tree_log = root_file.Get("Tree_log")
    tree_noise = root_file.Get("Tree_noise")
    
    entry1=find_noise_number(utc,tree_noise,1)
    entry2=find_noise_number(utc,tree_noise,5)
    entry3=find_noise_number(utc,tree_noise,9)
    entry4=find_noise_number(utc,tree_noise,13)
    entry5=find_noise_number(utc,tree_noise,17)
    
    all_info=[]
    
    
    for i in np.arange(LORA.nLORA):
        detname='Det'+str(1+i)
        det=tree_noise.GetBranch(detname)
        if i==0 or i==1 or i==2 or i==3:
            all_info.append(getNoiseV1(det, entry1))
        if i==4 or i==5 or i==6 or i==7:
            all_info.append(getNoiseV1(det, entry2))
        if i==8 or i==9 or i==10 or i==11:
            all_info.append(getNoiseV1(det, entry3))
        if i==12 or i==13 or i==14 or i==15:
            all_info.append(getNoiseV1(det, entry4))
        if i==16 or i==17 or i==18 or i==19:
            all_info.append(getNoiseV1(det, entry5))
    
    return all_info



def log_file(filename,data_dir):
    filepath=data_dir+filename+'.log'
    #log_file=open(data_dir+filename+'.log','r')
    lasa1_status=0
    lasa2_status=0
    lasa3_status=0
    lasa4_status=0
    lasa5_status=0
    LOFAR_trig=0

    with open(filepath,'r') as fp:
        line = fp.readline()
        cnt = 1
        while line:
            #print("Line {}: {}".format(cnt, line.strip()))
            if 'LOFAR trigger settings' in line:
                LOFAR_trig=int(line.strip().split(':')[1])
            
            if 'CS003:' in line:
                lasa1_status=int(line.strip().split(':')[1])
            if 'CS004:' in line:
                lasa2_status=int(line.strip().split(':')[1])
            if 'CS005:' in line:
                lasa3_status=int(line.strip().split(':')[1])
            if 'CS006:' in line:
                lasa4_status=int(line.strip().split(':')[1])
            if 'CS007:' in line:
                lasa5_status=int(line.strip().split(':')[1])
            
            line = fp.readline()
            cnt += 1


    '''
    print 'LOFAR trigger: ',LOFAR_trig
    print 'lasa 1: ',lasa1_status
    print 'lasa 2: ',lasa2_status
    print 'lasa 3: ',lasa3_status
    print 'lasa 4: ',lasa4_status
    print 'lasa 5: ',lasa5_status
    '''
    info={'LOFAR_trig':LOFAR_trig,'lasa1_status':lasa1_status,'lasa2_status':lasa2_status,'lasa3_status':lasa3_status,'lasa4_status':lasa4_status,'lasa5_status':lasa5_status}
    return info
