import numpy as np
from datetime import datetime
import os
import ROOT


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
        print i
        detname='Det'+str(1+i)
        det=tree_event.GetBranch(detname)
        info=getDataV1(det,event_index)
        all_info.append(info)

    return all_info
