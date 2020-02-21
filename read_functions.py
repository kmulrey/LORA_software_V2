import numpy as np
from datetime import datetime
import os
import ROOT
import LORAparameters as LORA
import os.path
from datetime import date

nTraceV1=4000
nDetV1=20
dtV1=2.5
timesV1=np.arange(0,4000*dtV1,dtV1)

avg_threshold=np.asarray([ 27.24867982,  27.08347496,  21.3469044,   22.45092564,  16.96505311,
16.62127466,  17.87311077,  17.90886191,  20.3131563,   20.87103187,
20.18339909,  22.10192716,  22.96687405,  23.04840668,  13.57128983,
13.92050076,  27.7231563,   15.69400607,  23.10057663,  24.20940819,])


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
    print 'day:    ',day
    print 'month:    ',month
    print 'year:    ',year

    filestr1= str(year)+str(month).zfill(2)+str(day).zfill(2)
    filestr3='-1'
    filestr4='-1'
    
    if day>1:
        filestr2= str(year)+str(month).zfill(2)+str(day-1).zfill(2)
        filestr4= str(year)+str(month).zfill(2)+str(day-2).zfill(2)
    
    if (month==10 or month==5 or month==7 or month==12) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(30).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(29).zfill(2)
    
    if (month==2 or month==4 or month==6 or month==8 or month==9 or month==11) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(31).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(30).zfill(2)
    
    if (month==3) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(28).zfill(2)
        filestr3= str(year)+str(month-1).zfill(2)+str(29).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(27).zfill(2)
    
    if (month==1) and day==1:
        filestr2= str(year-1)+str(12).zfill(2)+str(31).zfill(2)
        filestr4= str(year-1)+str(12).zfill(2)+str(30).zfill(2)

    #filestr3= str(year)+str(month).zfill(2)+str(day-1).zfill(2)
    
    filelist=[]
    
    for file in os.listdir(data_dir):
        if filestr1 in file or filestr2 in file or filestr3 in file or filestr4 in file:
            if '.log' in file:
                #print(os.path.join(data_dir, file))
                filelist.append(os.path.join(data_dir, file))


    found=0
    filetag='no'
    
    for i in np.arange(len(filelist)):
        #print filelist[i]
        #for i in np.arange(1):
        file=open(filelist[i])
        stamps=[]
        try:
            stamps=np.genfromtxt(file,skip_header=10,usecols=(2,3))
        except:
            print 'empty log'
        file.close()

        if len(stamps)>0:
            if len(stamps[(stamps.T[0]==timestamp)*(stamps.T[1]==ns_timestamp)])>0:
                
                filetag=filelist[i].split('raw/')[1].split('.')[0]
                found=1


    #print 'looking for this file: ',data_dir+filetag+'.root'

    if found==1 and os.path.exists(data_dir+filetag+'.root'):
        return filetag
    else:
        return 'no_match'

def find_tag_exception(timestamp,ns_timestamp,data_dir):
    #print '\n________________________\n'
    #timestamp = LORA4times[t][0]
    #ns_timestamp = LORA4times[t][1]
    
    print 'no standard log file avaliable, trying to look at .root files'
    
    dt_object = datetime.fromtimestamp(timestamp)
    year=dt_object.year
    month=dt_object.month
    day=dt_object.day
    hour=dt_object.hour
    minute=dt_object.minute
    second=dt_object.second
    
    #print year,month,day,hour
    #print int(timestamp),int(ns_timestamp)
    filestr1= str(year)+str(month).zfill(2)+str(day).zfill(2)
    filestr3='-1'
    filestr4='-1'
    
    if day>1:
        filestr2= str(year)+str(month).zfill(2)+str(day-1).zfill(2)
        filestr4= str(year)+str(month).zfill(2)+str(day-2).zfill(2)

    if (month==10 or month==5 or month==7 or month==12) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(30).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(29).zfill(2)

    if (month==2 or month==4 or month==6 or month==8 or month==9 or month==11) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(31).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(30).zfill(2)

    if (month==3) and day==1:
        filestr2= str(year)+str(month-1).zfill(2)+str(28).zfill(2)
        filestr3= str(year)+str(month-1).zfill(2)+str(29).zfill(2)
        filestr4= str(year)+str(month-1).zfill(2)+str(27).zfill(2)

    if (month==1) and day==1:
        filestr2= str(year-1)+str(12).zfill(2)+str(31).zfill(2)
        filestr4= str(year-1)+str(12).zfill(2)+str(30).zfill(2)




    filelist=[]
    
    for file in os.listdir(data_dir):
        if filestr1 in file or filestr2 in file or filestr3 in file or filestr4 in file:
            if '.root' in file:
                #print(os.path.join(data_dir, file))
                filelist.append(os.path.join(data_dir, file))
    print filelist
    
    found=0
    filetag='no'
    
    for i in np.arange(len(filelist)):
        root_file=ROOT.TFile.Open(filelist[i])
        try:
            tree_event = root_file.Get("Tree_event")
        
        
        
            event_index=-1
            event_index=find_entry_number(timestamp,ns_timestamp,tree_event)
            print 'did we find the index?   {0}'.format(event_index)
        
            if event_index>-1:
                found=1
                filetag=filelist[i].split('raw/')[1].split('.')[0]
                break
        except:
            print 'can\'t open the event leaf'
    
    if found==1:
        print 'returning file tag'
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

    #log_file=open(data_dir+filename+'.log','r')
    root_file=ROOT.TFile.Open(data_dir+filename+'.root')
    print 'reading root file: {0}'.format(data_dir+filename+'.root')
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

    #log_file=open(data_dir+filename+'.log','r')
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


def getLogV1(det,d, entry):
    nE= det.GetEntries()
    det.GetEntry(entry)
    
    
    YMD=det.GetLeaf('YMD').GetValue()
    GPS_time_stamp=det.GetLeaf('Time_stamp').GetValue()
    Threshold_low=det.GetLeaf('Channel_thres_low').GetValue()
    
    print '_________________________________________'
    print 'finding threshold:  {0}'.format(Threshold_low)
    if Threshold_low==0.0:
        print 'issue with Threshold_low {0}'.format(nE)
        thresh_avg=0
        thresh_count=0
        for i in np.arange(nE):
            det.GetEntry(i)
            thesh_temp=det.GetLeaf('Channel_thres_low').GetValue()
            if thesh_temp>0.0:
                thresh_avg=thresh_avg+thesh_temp
                thresh_count=thresh_count+1
        if thresh_count>0:
            Threshold_low=thresh_avg/(1.0*thresh_count)
        else:
            Threshold_low=avg_threshold[d]
        
    
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
    
    #log_file=open(data_dir+filename+'.log','r')
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
            all_info.append(getLogV1(det,i, entry1))
        if i==4 or i==5 or i==6 or i==7:
            all_info.append(getLogV1(det,i, entry2))
        if i==8 or i==9 or i==10 or i==11:
            all_info.append(getLogV1(det,i, entry3))
        if i==12 or i==13 or i==14 or i==15:
            all_info.append(getLogV1(det,i, entry4))
        if i==16 or i==17 or i==18 or i==19:
            all_info.append(getLogV1(det,i, entry5))
    
    return all_info






def return_noise_data(filename,utc,nsec,data_dir):
    
    #log_file=open(data_dir+filename+'.log','r')
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
    lasa1_status=-1
    lasa2_status=-1
    lasa3_status=-1
    lasa4_status=-1
    lasa5_status=-1
    LOFAR_trig=-1
    try:
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

    except:
        print 'can\'t find log file'
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






def return_event_V2(event_id,event_GPS, event_ns,event_data):

    Station=event_data['Station']
    Detector=event_data['Detector']
    Channel_Passed_Threshold=event_data['Channel_Passed_Threshold']
    Trigg_Threshold=event_data['Trigg_Threshold']
    Charge_Corrected=event_data['Charge_Corrected']
    Peak_Height_Corrected=event_data['Peak_Height_Corrected']
    Peak_Height_Raw=event_data['Peak_Height_Raw']
    Waveform_Raw=event_data['Waveform_Raw']
    Event_Id=event_data['Event_Id']
    Run_Id=event_data['Run_Id']
    GPS_Time_Stamp=event_data['GPS_Time_Stamp']
    CTD=event_data['CTD']
    nsec_Online=event_data['nsec_Online']
    HiSparc_Trigg_Pattern=event_data['HiSparc_Trigg_Pattern']
    HiSparc_Trigg_Condition=event_data['HiSparc_Trigg_Condition']


    dets=Detector[Event_Id==event_id]
    counts=Waveform_Raw[Event_Id==event_id]
    pulse_height=Peak_Height_Corrected[Event_Id==event_id]
    total_counts=Charge_Corrected[Event_Id==event_id]
    trigger_pattern=HiSparc_Trigg_Pattern[Event_Id==event_id]
    trigger_condition=HiSparc_Trigg_Condition[Event_Id==event_id]
    nsec=nsec_Online[Event_Id==event_id]
    ctd=CTD[Event_Id==event_id]
    gps=GPS_Time_Stamp[Event_Id==event_id]
    thresh=Trigg_Threshold[Event_Id==event_id]
    #print counts.shape
    all_info=[]
    log_all_info=[]
    

    for d in np.arange(LORA.nLORA):
        if (d+1) in dets:
            ind=np.where(dets==(d+1))[0][0]
            #print d+1, dets[ind], total_counts[ind]
            timestamp = date.fromtimestamp(gps[ind])
            ymd=int(str(timestamp.year).zfill(2)+str(timestamp.month).zfill(2)+str(timestamp.day).zfill(2))
            info={'det':dets[ind],'ymd':ymd,'gps':gps[ind],'ctd':ctd[ind],'nsec':nsec[ind],'trigg_condition':trigger_condition[ind],'trigg_pattern':trigger_pattern[ind],'total_counts':total_counts[ind],'pulse_height':pulse_height[ind],'pulse_width':0,'counts':counts[ind]}
            log_info={'threshold':thresh[ind]}
            all_info.append(info)
            log_all_info.append(log_info)
        else:
            info={'det':d+1,'ymd':0,'gps':0,'ctd':0,'nsec':0,'trigg_condition':0,'trigg_pattern':0,'total_counts':0,'pulse_height':0,'pulse_width':0,'counts':np.zeros([4000])}
            log_info={'threshold':0}

            all_info.append(info)
            log_all_info.append(log_info)



    return all_info,log_all_info



def return_second_data_V2(event_id,event_GPS, event_ns,osm_data_hisparc,osm_data_aera):


    Station_H=osm_data_hisparc['Station']
    Master_Or_Slave_H=osm_data_hisparc['Master_Or_Slave']
    
    
    GPS_Time_Stamp_H=osm_data_hisparc['GPS_Time_Stamp']
    Sync_Error_H=osm_data_hisparc['Sync_Error']
    Quant_Error_H=osm_data_hisparc['Quant_Error']
    CTP_H=osm_data_hisparc['CTP']

    all_info=[]
    all_info1=[]
    all_info2=[]
    '''
    for t in np.arange(3):
        #for i in np.arange(1):
        for i in np.arange(LORA.nLASA):
            lasa=i+1
            #master
            #print 'getting OSM for {0}'.format(i+1)
            gpsM= GPS_Time_Stamp[(GPS_Time_Stamp==(event_GPS+t))*(Station==i+1)*(Master_Or_Slave==0)]
            if len(gpsM>0):
                syncM= Sync_Error[(GPS_Time_Stamp==(event_GPS+t))*(Station==i+1)*(Master_Or_Slave==0)]
                quantM= Quant_Error[(GPS_Time_Stamp==(event_GPS+t))*(Station==i+1)*(Master_Or_Slave==0)]
                ctpM= CTP[(GPS_Time_Stamp==(event_GPS+t))*(Station==i+1)*(Master_Or_Slave==0)]
                #print gpsM
                timestamp = date.fromtimestamp(gpsM)
                ymdM=int(str(timestamp.year).zfill(2)+str(timestamp.month).zfill(2)+str(timestamp.day).zfill(2))
                #slave
                gpsS= GPS_Time_Stamp[(GPS_Time_Stamp==(event_GPS+t))*(Station==i+1)*(Master_Or_Slave==1)]
                syncS= Sync_Error[(GPS_Time_Stamp==(event_GPS+t))*(Station==i+1)*(Master_Or_Slave==1)]
                quantS= Quant_Error[(GPS_Time_Stamp==(event_GPS+t))*(Station==i+1)*(Master_Or_Slave==1)]
                ctpS= CTP[(GPS_Time_Stamp==(event_GPS+t))*(Station==i+1)*(Master_Or_Slave==1)]
                timestamp = date.fromtimestamp(gpsM)
                ymdS=int(str(timestamp.year).zfill(2)+str(timestamp.month).zfill(2)+str(timestamp.day).zfill(2))
        
        
                #print gpsM,syncM,ctpM
                
                info={'lasa':lasa,'YMD_M':ymdM,'GPS_time_stamp_M':gpsM,'sync_M':syncM,'CTP_M':ctpM,'quant_M':quantM,'YMD_S':ymdS,'GPS_time_stamp_S':gpsS,'sync_S':syncS,'CTP_S':ctpS,'quant_S':quantS}
            else:
                info={'lasa':lasa,'YMD_M':np.asarray([0]),'GPS_time_stamp_M':np.asarray([0]),'sync_M':np.asarray([0]),'CTP_M':np.asarray([0]),'quant_M':np.asarray([0]),'YMD_S':np.asarray([0]),'GPS_time_stamp_S':np.asarray([0]),'sync_S':np.asarray([0]),'CTP_S':np.asarray([0]),'quant_S':np.asarray([0])}
            
            if t==0:
                all_info.append(info)
            if t==1:
                all_info1.append(info)
            if t==2:
                all_info2.append(info)
    '''
    return all_info,all_info1,all_info2
