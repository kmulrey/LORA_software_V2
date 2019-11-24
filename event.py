import numpy as np
import LORAparameters as LORA
import detector as det
from ROOT import TH1F
from ROOT import TF1
from ROOT import TH2F
from ROOT import TF2


class Event:
    def __init__(self, name,type):
        self.name = name  # LOFAR event ID
        self.type = type  # V1 or V2
    
    theta=0
    elevation=0
    phi=0
    fit_theta=0
    fit_elevation=0
    fit_phi=0
    fit_theta_err=0
    fit_phi_err=0
    x_core=0
    y_core=0
    z_core=0


def find_counts(detector):
    
    # find background
    background=detector.counts[0:int(LORA.BG_No_Bin)]
    background_mean=np.average(background)
    background_rms=np.std(background)
    detector.trace_rms=background_rms
    detector.trace_mean=background_mean

    if background_rms<10.0:
        corrected=detector.counts-background_mean
        peak=np.max(corrected)
        max_bin=np.argmax(corrected)
        if peak<LORA.Max_ADC_Count:
            BIN_S=int(max_bin-detector.B_min) # start integration
            BIN_E=int(max_bin+(int(LORA.Sig_Time_Window/2.5))) # end integration

            total_count=np.sum(corrected[BIN_S:BIN_E])-np.sum(corrected[0:(int(LORA.Sig_Time_Window/2.5))])
            detector.trace_int_counts=total_count
    #print '{0:.2f}  -> {1:.2f}  :   {2:0.2f}'.format(total_count, detector.total_counts,total_count/detector.total_counts)



def retrive_sat_signal(detector):
    test=0


def get_arrival_time(detector):
    cut=LORA.Det_Thres*detector.trace_rms+detector.trace_mean
    flag=0
    for i in np.arange(LORA.nTrace):
        if detector.counts[i]>cut and flag==0:
            if i<400:
                continue
            else:
                detector.threshold_time=i*2.5*10  # unit of 0.1 ns
                flag=1

def get_event_timestamp(detector,lasa):
    
    if lasa.sec_flag!=1 and lasa.CTP[1]>0:
        detector.event_time_stamp=10*int((lasa.sync[0]+lasa.quant[1]+(1.0*detector.ctd/lasa.CTP[1])*(1000000000.0-lasa.quant[1]+lasa.quant[2])))


def cal_event_timestamp(detectors,lasa):
    #print '_____________________________________'
    #print 'lasa number: {0}'.format(lasa.number)
    lasa_ind=int(lasa.number-1)
    #print 'lasa index: {0}'.format(lasa_ind)
    #print 'trigger condition: {0}'.format(detectors[4*lasa_ind].trigg_condition)
    trigg_cond=detectors[4*lasa_ind].trigg_condition
    thresh_times=np.asarray([detectors[4*int(lasa.number-1)].threshold_time,detectors[4*int(lasa.number-1)+1].threshold_time,detectors[4*int(lasa.number-1)+2].threshold_time,detectors[4*int(lasa.number-1)+3].threshold_time])
    thresh_use=1.0*thresh_times[thresh_times!=0]
    args=np.argsort(thresh_use)
    #print int(args[detectors[0].trigg_condition]-1)
    if len(thresh_use)>1:
        #print 'sorted: {0}'.format(np.sort(thresh_use))
        #print thresh_use
        #print args
        #print 'use index: {0}'.format(args[int(trigg_cond)-1])
        trigg_time=thresh_use[args[int(trigg_cond)-1]]
        #print thresh_times,trigg_time

        for i in np.arange(4):
            if detectors[4*int(lasa_ind)+i].threshold_time>0:
            
                detectors[4*int(lasa_ind)+i].cal_time=detectors[4*int(lasa_ind)+i].threshold_time-trigg_time+detectors[lasa_ind*4].event_time_stamp+det.cable_delay[4*int(lasa_ind)+i]
                detectors[4*int(lasa_ind)+i].final_event_time=detectors[4*int(lasa_ind)+i].cal_time
                # maybe this has to be corrected for wrap-around seconds

            #print '{0}   {1}  {2}  {3}  {4}  {5}'.format(4*int(lasa_ind)+i,int(detectors[4*int(lasa_ind)+i].cal_time),detectors[4*int(lasa_ind)+i].threshold_time, trigg_time,detectors[4*int(lasa_ind)+i].event_time_stamp,det.cable_delay[4*int(lasa_ind)+i])




def do_arrival_time_diff(detectors):
    event_times=np.zeros([LORA.nDetA])
    event_weight=np.zeros([LORA.nDetA])

    ind_0=1000
    time_min=1e10
    for i in np.arange(LORA.nDetA):
        event_times[i]=detectors[i].final_event_time
        event_weight[i]=detectors[i].trace_int_counts
        #print i, event_times[i],event_weight[i]
        if event_times[i]<time_min and event_times[i]!=0:
            time_min=event_times[i]
            ind_0=i
    #print ind_0, event_times[ind_0]

    for i in np.arange(LORA.nDetA):
        if(event_times[i]>0 and event_weight[i]>0):
            detectors[i].cdt=(event_times[i]-event_times[ind_0])*0.1*(1.e-9*LORA.vel_light)


def do_arrival_direction(detectors,event):

    P=0
    Q=0
    R=0
    S=0
    W=0
    T1=0
    S1=0
    S2=0
    S3=0
    S4=0
    S5=0
    S6=0
    
    counter=0
    
    for i in np.arange(LORA.nDetA):
        if detectors[i].cdt>=0:
            counter=counter+1
            S=S+detectors[i].x_cord**2
            W=W+detectors[i].y_cord**2
            
            T1=T1+detectors[i].x_cord
            S1=S1+detectors[i].x_cord*detectors[i].cdt
            S2=S2+detectors[i].x_cord*detectors[i].y_cord
            S3=S3+detectors[i].y_cord*detectors[i].cdt
            S4=S4+detectors[i].y_cord
            S5=S5+detectors[i].cdt
            S6=S6+1

    P=(S1*S2)/S
    Q=(T1*S2)/S
    R=-((S2*S2)/S)+W

    t0=(T1*S1*R-R*S*S5+T1*P*S2-T1*S2*S3-S*S4*P+S*S4*S3)/(T1*Q*S2-T1*S2*S4+R*T1*T1-S*S4*Q+S*S4*S4-R*S*S6)
    m=(P-t0*Q-S3+t0*S4)/R
    l=(-S1/S)-((P*S2)/(R*S))+((t0*Q*S2)/(R*S))+((S2*S3)/(R*S))-((t0*S2*S4)/(R*S))+((t0*T1)/S)
    n=np.sqrt(1.0-(l*l+m*m))

    theta=(np.arcsin(np.sqrt(l*l+m*m)))*(180.0/np.pi)    #Zenith in degrees (from vertical direction +Z)
    phi=(np.arccos(m/np.sqrt(l*l+m*m)))*(180.0/np.pi)    #in degrees (Eastward from North)
    if l<0:
        phi=360.0-phi



    event.theta=theta
    event.phi=phi
    event.elevation=90.0-theta
    print 'theta: {0:.2f}   phi: {1:.2f}    el: {2:.2f}'.format(event.theta,event.phi,event.elevation)


def do_COM_core(detectors,event):

    print 'COM core'
    x=0
    Event_Size=0
    SumX=0
    SumY=0
    
    for i in np.arange(LORA.nDetA):
    
        x=detectors[i].trace_int_counts/detectors[i].gain/(LORA.Det_Area*np.cos(event.theta*np.pi/180.0))
        
        if(detectors[i].final_event_time>0 and x>LORA.Density_Cut):
        
            Event_Size=Event_Size+detectors[i].trace_int_counts/detectors[i].gain
            SumX=SumX+detectors[i].trace_int_counts/detectors[i].gain*detectors[i].x_cord
            SumY=SumY+detectors[i].trace_int_counts/detectors[i].gain*detectors[i].y_cord
        

        
    
    event.x_core=SumX/Event_Size
    event.y_core=SumY/Event_Size
    event.z_core=0     #we assume that all the LORA detectors are at z=0

    print 'core: ({0:.2f}, {1:.2f}, {2:.2f})'.format(event.x_core,event.y_core,event.z_core)


def find_density(detectors,event):

    
    for i in np.arange(LORA.nDetA):
        detectors[i].density=detectors[i].trace_int_counts/detectors[i].gain/(LORA.Det_Area*np.cos(event.theta*np.pi/180.0))
        detectors[i].err_density=np.power(detectors[i].density,0.5) ;    #Assuming possionian error
        if detectors[i].density<=LORA.Density_Cut:
            detectors[i].density=0
            detectors[i].err_density=0



def func_plane(x,par):    #//Plane function to fit the shower front: (Parameters: t0,l,m)
    return par[0]-par[1]*x[0]-par[2]*x[1]


def fit_arrival_direction(detectors,event):
    print 'fit arrival direction'
        
    plane=TH2F("Shower plane","",1000,-500,500,1000,-500,500)

    for i in np.arange(LORA.nDetA):
        if detectors[i].cdt>=0:
            plane.SetBinContent(plane.GetXaxis().FindBin(detectors[i].x_cord),plane.GetYaxis().FindBin(detectors[i].y_cord),detectors[i].cdt)

    fit_plane= TF2("fit_plane",func_plane,-500,500,-500,500,3)
    fit_plane.SetParameter(0,100) ;
    fit_plane.SetParameter(1,0.2) ;
    fit_plane.SetParameter(2,0.2) ;
    fit_plane.SetParLimits(1,-1,1) ;
    fit_plane.SetParLimits(2,-1,1) ;
    plane.Fit(fit_plane,"Q0") ;
    t0=fit_plane.GetParameter(0)
    l=fit_plane.GetParameter(1)
    m=fit_plane.GetParameter(2)
    n=np.sqrt(1.0-(l*l+m*m))
    err_l=fit_plane.GetParError(1)
    err_m=fit_plane.GetParError(2)
    theta=(np.arcsin(np.sqrt(l*l+m*m)))*(180.0/np.pi)#//Zenith in degrees (from vertical direction +Z)
    phi=(np.arccos(m/np.sqrt(l*l+m*m)))*(180.0/np.pi)# ;    //in degrees (Eastward from North)
    if(l<0):
        phi=360.0-phi

    err_theta=np.sqrt(np.power(l*err_l/(n*n),2)+np.power(m*err_m/(n*n),2))/np.tan(theta*np.pi/180.0)
    err_theta=err_theta*(180.0/np.pi)
    df_f2=np.power((l*l)/(m*(l*l+m*m)),2)*np.power(err_m,2)+np.power(l/(l*l+m*m),2)*np.power(err_l,2)
    err_phi=np.sqrt(df_f2/np.power(np.tan(phi*np.pi/180.0),2))
    err_phi=err_phi*(180.0/np.pi)
    if phi>360.0:
        phi=phi=360

    print 'fit    theta: {0:.2f}  theta_err: {2:.2f}    phi: {1:.2f}   phi_err: {3:.2f} '.format(theta,phi,err_theta,err_phi)
    event.fit_theta=theta
    event.fit_phi=phi
    event.fit_elevation=90.0-theta
    event.fit_theta_err=err_theta
    event.fit_phi_err=err_phi
