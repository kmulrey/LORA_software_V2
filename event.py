import numpy as np
import LORAparameters as LORA
import detector as det
from ROOT import TH1F
from ROOT import TF1
from ROOT import TH2F
from ROOT import TF2
from ROOT import TMath
from ROOT import TVirtualFitter
import math

atm_flag=0
atm_l=11
f_int=np.zeros([atm_l])
f_logN_Ref=np.zeros([atm_l])
err1=np.zeros([atm_l])
f_X0=np.zeros([atm_l])
err2=np.zeros([atm_l])
f_lamb=np.zeros([atm_l])
err3=np.zeros([atm_l])

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
    x_core_err=0
    y_core_err=0
    z_core=0
    UTC_min=0
    nsec_min=0
    energy=0
    energy_err=0
    Rm=0
    fit_elevation_err=0
    fit_phi_err=0
    Ne=0
    Ne_err=0
    CorCoef_xy=0
    Ne_RefA=0
    NeErr_RefA=0
    Energy_RefA=0
    EnergyErr_RefA=0




def read_attenuation():
    try:
        data=np.genfromtxt(LORA.atm_file)
        return data
    except:
        print 'problem reading atm file'
        return 0

def find_counts(detector):

    if detector.number<=20:
        # find background
        background=detector.counts[0:int(LORA.BG_No_Bin)]
        background_mean=np.average(background)
        background_rms=np.std(background)
        detector.trace_rms=background_rms
        detector.trace_mean=background_mean
        detector.peak=np.max(detector.counts)
    
        if(background_mean<100 and background_mean>0):
            detector.corrected_threshold=detector.threshold-background_mean
        
            print 'finding threshold from real background: {0}  {1}  {2}'.format(detector.corrected_threshold,detector.threshold,background_mean)
        
        else:
            detector.corrected_threshold=detector.threshold-detector.sec_mean
            print 'finding threshold from second background: {0}  {1}  {2}'.format(detector.corrected_threshold,detector.threshold,detector.sec_mean)
        
        if detector.corrected_threshold<0:
            detector.corrected_threshold=detector.threshold-detector.sec_mean
        
        if detector.corrected_threshold<0:
            print '~*~***~*~*~*~**~'
            print 'what the heck is going on with this?????'
            print detector.threshold,detector.sec_mean,background_mean
        
        if background_rms<10.0:
            corrected=detector.counts-background_mean
            peak=np.max(corrected)
            max_bin=np.argmax(corrected)
            if peak<LORA.Max_ADC_Count:
                BIN_S=int(max_bin-detector.B_min) # start integration
                BIN_E=int(max_bin+(int(LORA.Sig_Time_Window/2.5))) # end integration
                total_count=np.sum(corrected[BIN_S:BIN_E])-np.sum(corrected[0:(int(LORA.Sig_Time_Window/2.5))])
            else:
                total_count=0
                BIN_S=int(max_bin-detector.B_min) # start integration
                BIN_E=int(max_bin+(int(LORA.Sig_Time_Window/2.5))) # end integration
                total_count=np.sum(corrected[BIN_S:BIN_E])-np.sum(corrected[0:(int(LORA.Sig_Time_Window/2.5))])

            if total_count>0:
                detector.trace_int_counts=total_count

            detector.corrected_peak=peak

    else:
        print 'counts new version : {0}'.format(detector.number)

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
    print '_______event timestamp______'
    print lasa.CTP
    print lasa.sec_flag
    if lasa.sec_flag!=1 and lasa.CTP[1]>0:
        detector.event_time_stamp=10*(lasa.sync[0]+lasa.quant[1]+(1.0*detector.ctd/lasa.CTP[1])*(1000000000.0-lasa.quant[1]+lasa.quant[2]))
        print 'doing real time stamp '
        print detector.event_time_stamp
    else:
        print 'doing est. time stamp '

        detector.event_time_stamp=10*detector.nsec#10*(1.0*detector.ctd/200000000.0 )*(1000000000.0)
        print detector.event_time_stamp

def get_event_timestamp_V2(detector,lasa):
    
    #print detector.number, detector.number%2
    if detector.number%2==1:
        if lasa.sec_flag!=1:
            detector.event_time_stamp=10*(lasa.sync_M[0]+lasa.quant_M[1]+(1.0*detector.ctd/lasa.CTP_M[1])*(1000000000.0-lasa.quant_M[1]+lasa.quant_M[2]))
            #print detector.event_time_stamp
        else:
            #print 'flagged event'
            detector.event_time_stamp=10*detector.nsec
            #print detector.event_time_stamp

        
    if detector.number%2==0:
        if lasa.sec_flag!=1:
            detector.event_time_stamp=10*(lasa.sync_S[0]+lasa.quant_S[1]+(1.0*detector.ctd/lasa.CTP_S[1])*(1000000000.0-lasa.quant_S[1]+lasa.quant_S[2]))
            #print detector.event_time_stamp
        else:
            #print 'flagged event'
            detector.event_time_stamp=10*detector.nsec
            #print detector.event_time_stamp


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
        #print 'trigger condition: ',trigg_cond
        #print 'sorted: {0}'.format(np.sort(thresh_use))
        #print thresh_use
        #print args
        #print 'use index: {0}'.format(args[int(trigg_cond)-1])
        try:
            trigg_time=thresh_use[args[int(trigg_cond)-1]]
        except: #this is in here becasue in a few case the length of thresh_use was less than trigg. condition
            trigg_time=np.sort(thresh_use)[len(thresh_use)-1]
        #print thresh_times,trigg_time

        for i in np.arange(4):
            if detectors[4*int(lasa_ind)+i].threshold_time>0:
            
                detectors[4*int(lasa_ind)+i].cal_time=detectors[4*int(lasa_ind)+i].threshold_time-trigg_time+detectors[lasa_ind*4].event_time_stamp+det.cable_delay[4*int(lasa_ind)+i]
                detectors[4*int(lasa_ind)+i].final_event_time=detectors[4*int(lasa_ind)+i].cal_time
                #print detectors[4*int(lasa_ind)+i].final_event_time

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
    
    
    
    print 'direction params:  {0} {1}  {2}  {3}'.format(t0,m,l,n)

    if l*l+m*m<1:

        theta=(np.arcsin(np.sqrt(l*l+m*m)))*(180.0/np.pi)    #Zenith in degrees (from vertical direction +Z)
        phi=(np.arccos(m/np.sqrt(l*l+m*m)))*(180.0/np.pi)    #in degrees (Eastward from North)
    else:
        theta=0.0
        phi=(np.arccos(m/np.sqrt(l*l+m*m)))*(180.0/np.pi)    #in degrees (Eastward from North)

        
    if l<0:
        phi=360.0-phi



    event.theta=theta
    event.phi=phi
    event.elevation=90.0-theta
    print 'theta: {0:.2f}   phi: {1:.2f}    el: {2:.2f}'.format(event.theta,event.phi,event.elevation)


def do_COM_core(detectors,event):

    #print 'COM core'
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
    #print 'fit arrival direction'
        
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

    
    if math.isnan(theta)==True or math.isnan(phi)==True:
        event.fit_theta=event.theta
        event.fit_phi=event.phi
        event.fit_elevation=90.0-theta
        event.fit_theta_err=err_theta
        event.fit_elevation_err=err_theta
        event.fit_phi_err=err_phi
        print 'fitting doesn\'t work'
    else:
        event.fit_theta=theta
        event.fit_phi=phi
        event.fit_elevation=90.0-theta
        event.fit_theta_err=err_theta
        event.fit_elevation_err=err_theta
        event.fit_phi_err=err_phi

    print 'fit    theta: {0:.2f}  phi: {1:.2f} '.format(theta,phi,err_theta,err_phi)


def theta_phi(theta,phi,psi,x0,y0,z0):

    #1st ROTATION: counterclockwise from X-axis(N) for angle 'phi' about Z-axis.
    x1= x0*np.cos(phi)+y0*np.sin(phi)
    y1=-x0*np.sin(phi)+y0*np.cos(phi)
    z1= z0
    #-------------xxx----------------
    #2nd ROTATION: clockwise from Z-axis for angle 'theta' about Y-axis(W).
    x2= x1*np.cos(theta)-z1*np.sin(theta) ;
    y2= y1
    z2= x1*np.sin(theta)+z1*np.cos(theta) ;
    #-------------xxx----------------
    #3rd ROTATION: counterclockwise from X-axis(N) for angle 'psi' about Z-axis.
    x3= x2*np.cos(psi)+y2*np.sin(psi)
    y3=-x2*np.sin(psi)+y2*np.cos(psi)
    z3= z2
    #-------------xxx----------------
    return x3,y3,z3


def back_theta_phi(theta,phi,psi,x0,y0,z0):
    #1st ROTATION:clockwise from X-axis(N) for angle 'psi' about Z-axis.
    x1= x0*np.cos(psi)-y0*np.sin(psi)
    y1= x0*np.sin(psi)+y0*np.cos(psi)
    z1= z0
    #--------------xxx------------------
    #2nd ROTATION: counterclockwise from Z-axis for angle 'theta' about Y-axis(W).
    x2= x1*np.cos(theta)+z1*np.sin(theta)
    y2= y1
    z2= -x1*np.sin(theta)+z1*np.cos(theta)
    #--------------xxx-------------------
    #3rd ROTATION:clockwise from X-axis(N) for angle 'phi' about Z-axis.
    x3= x2*np.cos(phi)-y2*np.sin(phi) ;
    y3= x2*np.sin(phi)+y2*np.cos(phi) ;
    z3= z2 ;
    #--------------xxx------------------
    return x3,y3,z3



def func_nkg_show(x,par):    #NKG function for fitting 5 parameters: (x_core,y_core,Ne,rM,s)

    return (par[0]/np.power(par[1],2))*(TMath.Gamma(4.5-par[2]))/(2*np.pi*(TMath.Gamma(par[2]))*(TMath.Gamma(4.5-2*par[2])))*np.power(np.sqrt(np.power(par[3]-x[0],2)+np.power(par[4]-x[1],2))/par[1],par[2]-2)*np.power(1.0+np.sqrt(np.power(par[3]-x[0],2)+np.power(par[4]-x[1],2))/par[1],par[2]-4.5)

def func_nkg_ind_show(x,par):    #//NKG function for fitting 3 parameters: (Ne,rM,s)

    return (par[0]/np.power(par[1],2))*(TMath.Gamma(4.5-par[2]))/(2*np.pi*(TMath.Gamma(par[2]))*(TMath.Gamma(4.5-2*par[2])))*np.power(x[0]/par[1],par[2]-2)*np.power(1.0+x[0]/par[1],par[2]-4.5)


def TF2_Fit_NKG(h2,N0,r0,s0,x0,y0,x_min,x_max,y_min,y_max):
    
    fitf2= TF2('fitf2',func_nkg_show,x_min,x_max,y_min,y_max,5)
    fitf2.SetParameter(0,N0) #        Setting starting values (Ne)
    fitf2.SetParameter(1,r0) #            ,,    ,,     (r_M)
    fitf2.SetParameter(2,s0) #            ,,    ,,     (s)
    fitf2.SetParameter(3,x0) #            ,,    ,,     (x_core)
    fitf2.SetParameter(4,y0) #            ,,    ,,     (y_core)
    
    #//fitf2->SetParLimits(1,10,300) ;        //Setting r_M limits
    fitf2.SetParLimits(1,r0,r0)#        //Setting r_M limits
    #//fitf2->SetParLimits(2,1,2.5) ;        //Setting Age limits (keeping fixed)
    fitf2.SetParLimits(2,s0,s0)#        //Setting Age limits (Keeping fixed)
    h2.Fit(fitf2,"Q0")
    #h2SetXTitle("X (m)") ;
    #h2->GetXaxis()->CenterTitle() ;
    #h2->SetYTitle("Y (m)") ;
    #h2->GetYaxis()->CenterTitle() ;
    #//h2->Draw("HIST") ;
    #//fitf2->DrawCopy("surf1 same");
    #gStyle->SetOptFit(1111) ;
    Ne_fit=fitf2.GetParameter(0)
    rM_fit=fitf2.GetParameter(1)
    s_fit=fitf2.GetParameter(2)
    x_core_fit=fitf2.GetParameter(3)
    y_core_fit=fitf2.GetParameter(4)
    x_core_fit_err=fitf2.GetParError(3)
    y_core_fit_err=fitf2.GetParError(4)
    
    
    
    fitter=TVirtualFitter.GetFitter()
    
  
    
    corMatrix=np.zeros([5,5])
    
    sig_i=0
    sig_j=0
    cov_ij=0
    
    for i in np.arange(4):

        sig_i=np.sqrt(fitter.GetCovarianceMatrixElement(i,i))
        for j in np.arange(4):
            sig_j=np.sqrt(fitter.GetCovarianceMatrixElement(j,j))
            cov_ij=fitter.GetCovarianceMatrixElement(i,j)
            corMatrix[i][j]=cov_ij/(sig_i*sig_j)
            if(sig_i==0 or sig_j==0):
                corMatrix[i][j]=0
       

    corCoef_xy=corMatrix[1][2] #Correlation coeff. between the x & y core position
    return Ne_fit,rM_fit,s_fit,x_core_fit,y_core_fit,x_core_fit_err,y_core_fit_err, corCoef_xy




def TF1_Fit_NKG(h1,N0,r0,s0,R_Min,R_Max):

    fit_R_Min=0
    if(R_Min<5):
        fit_R_Min=R_Min-2
    else:
        fit_R_Min=R_Min-5

    if(fit_R_Min<0):
        fit_R_Min=0.1


    fitf1= TF1("fitf1",func_nkg_ind_show,fit_R_Min,R_Max+5,3) ;

    fitf1.SetParameter(0,N0)        #//Setting starting values (Ne)
    fitf1.SetParameter(1,r0)        #//Setting starting values (r_M)
    fitf1.SetParameter(2,s0)       #//Setting starting values (s)
    fitf1.SetParLimits(1,1,300)        #//Setting r_M limits
    #fitf1->SetParLimits(1,r0,r0) ;        //Setting r_M limits
    # fitf1->SetParLimits(2,1,2.5) ;        //Setting Age limits (keeping fixed)
    fitf1.SetParLimits(2,s0,s0)       # //Setting Age limits (keeping fixed)
    fitf1.SetParName(0,"N_{e}")
    fitf1.SetParName(1,"r_{M}")
    fitf1.SetParName(2,"s")

    #//h1->Fit("fitf1","V+","",fit_R_Min-5,R_Max+5) ;
    h1.Fit("fitf1","Q","",fit_R_Min,R_Max+5)
    #h1->SetYTitle("Charged particle density (m^{-2} )") ;
    #h1->GetYaxis()->CenterTitle() ;
    #h1->SetXTitle("Distance from shower axis (m)") ;
    #h1->GetXaxis()->CenterTitle() ;
    #h1->SetMarkerColor(4);
    #h1->SetMarkerStyle(20);
    #h1->SetMarkerSize(1);
    #h1->SetFillColor(0);
    ##//h1->Draw("HIST") ;
    #//gStyle->SetOptTitle(kFALSE) ;
    #gStyle->SetOptStat(10) ;
    #gStyle->SetStatColor(0) ;
    #gStyle->SetTitleFillColor(0) ;
    #gStyle->SetOptFit(0111) ;
    Ne=fitf1.GetParameter(0)
    Ne_err=fitf1.GetParError(0)
    rM=fitf1.GetParameter(1)
    s=fitf1.GetParameter(2)
    
    return Ne,rM,s,Ne_err




def fit_NKG(detectors,event):
    print 'fit core position'
    
    # directions for fitting
    theta=(event.fit_theta)*np.pi/180.0
    phi=(90.0+360.0-event.fit_theta)*np.pi/180.0  #Northward (anticlockwise) from East (X-axis)
    psi=2*np.pi-phi
    
    x_shower=np.zeros([LORA.nDetA])
    y_shower=np.zeros([LORA.nDetA])
    x_pos=np.zeros([LORA.nDetA])
    y_pos=np.zeros([LORA.nDetA])
    z_shower=np.zeros([LORA.nDetA])

    temp_det_array=[]
    temp_den_array=[]
    temp_x_pos=[]
    temp_y_pos=[]

    # find core position in shower plane for first guess for fit. 4 densest detectors are used
    
    for i in np.arange(LORA.nDetA):
        x_shower[i],y_shower[i],z_shower[i]=theta_phi(theta,phi,psi,detectors[i].x_cord,detectors[i].y_cord,detectors[i].z_cord)
        x_pos[i]=detectors[i].x_cord
        y_pos[i]=detectors[i].y_cord
        if detectors[i].density>=LORA.Density_Cut and detectors[i].density<=LORA.Density_Cut_High:
            temp_den_array.append(detectors[i].density)
            temp_det_array.append(i)
            temp_x_pos.append(x_shower[i])
            temp_y_pos.append(x_shower[i])

    ind=np.argsort(np.asarray(temp_den_array))[::-1]

    temp_den_array=np.asarray(temp_den_array)
    temp_det_array=np.asarray(temp_det_array)
    temp_x_array=np.asarray(temp_x_pos)
    temp_y_array=np.asarray(temp_y_pos)

    temp_total_den=temp_den_array[ind[0]]+temp_den_array[ind[1]]+temp_den_array[ind[2]]+temp_den_array[ind[3]]
    x_core=(x_shower[temp_det_array[ind[0]]]*temp_den_array[ind[0]]+x_shower[temp_det_array[ind[1]]]*temp_den_array[ind[1]]+x_shower[temp_det_array[ind[2]]]*temp_den_array[ind[2]]+x_shower[temp_det_array[ind[3]]]*temp_den_array[ind[3]])/temp_total_den
    y_core=(y_shower[temp_det_array[ind[0]]]*temp_den_array[ind[0]]+y_shower[temp_det_array[ind[1]]]*temp_den_array[ind[1]]+y_shower[temp_det_array[ind[2]]]*temp_den_array[ind[2]]+y_shower[temp_det_array[ind[3]]]*temp_den_array[ind[3]])/temp_total_den

    # find min/max detector positions

    nBinsX=int(((np.max(x_pos)+10)-(np.min(x_pos)-10))/LORA.Bin_Size_X)
    nBinsY=int(((np.max(y_pos)+10)-(np.min(y_pos)-10))/LORA.Bin_Size_Y)
    x_min=(np.min(x_pos)-10)
    x_max=(np.max(x_pos)+10)
    y_min=(np.min(y_pos)-10)
    y_max=(np.max(y_pos)+10)

    nBinsX_shower=int(((np.max(x_shower)+10)-(np.min(x_shower)-10))/LORA.Bin_Size_X)
    nBinsY_shower=int(((np.max(y_shower)+10)-(np.min(y_shower)-10))/LORA.Bin_Size_Y)
    x_min_shower=(np.min(x_shower)-150)
    x_max_shower=(np.max(x_shower)+150)
    y_min_shower=(np.min(y_shower)-150)
    y_max_shower=(np.max(y_shower)+150)



    evt_display = TH2F('Event','Event',nBinsX,x_min,x_max,nBinsY,y_min,y_max) # ground plane
    evt_display_shower = TH2F('Event_shower','Event_shower',nBinsX_shower,x_min_shower,x_max_shower,nBinsY_shower,y_min_shower,y_max_shower) # shower plane
    lat_den_show=TH1F('Lateral_density','Lateral_density',LORA.No_Bin_R,LORA.min_R,LORA.max_R)
    time_display = TH2F('time','time',nBinsX,x_min,x_max,nBinsY,y_min,y_max)    #ground plane

    rho=np.zeros([LORA.nDetA])
    rho_shower=np.zeros([LORA.nDetA])
    
    rho_err=np.zeros([LORA.nDetA])
    nDet_triggered=0
    f_den=0
    shower_size=0

    for i in np.arange(LORA.nDetA):
        rho[i]=detectors[i].density
        rho_err[i]=detectors[i].err_density
        if rho[i]>=LORA.Density_Cut and rho[i]<=LORA.Density_Cut_High:
            nDet_triggered=nDet_triggered+1
            evt_display.SetBinContent(evt_display.GetXaxis().FindBin(x_pos[i]),evt_display.GetYaxis().FindBin(y_pos[i]),rho[i])
            evt_display_shower.SetBinContent(evt_display_shower.GetXaxis().FindBin(x_shower[i]),evt_display_shower.GetYaxis().FindBin(y_shower[i]),rho[i])
            f_den=f_den+(1.0/np.power(LORA.rM,2))*(TMath.Gamma(4.5-LORA.Age))/(2*np.pi*(TMath.Gamma(1.0))*(TMath.Gamma(4.5-2*LORA.Age)))*np.power(np.sqrt(np.power(x_core-x_shower[i],2)+np.power(y_core-y_shower[i],2))/LORA.rM,LORA.Age-2)*np.power(1.0+np.sqrt(np.power(x_core-x_shower[i],2)+np.power(y_core-y_shower[i],2))/LORA.rM,(LORA.Age-4.5))

            shower_size=shower_size+rho[i]







    # first iteration of fit -> doing Ne, xcore, ycore
    Ne_fit,rM_fit,s_fit,x_core_fit,y_core_fit,x_core_fit_err,y_core_fit_err, corr_coef_xy=TF2_Fit_NKG(evt_display_shower,shower_size/f_den,LORA.rM,LORA.Age,x_core,y_core,x_min_shower,x_max_shower,y_min_shower,y_max_shower)
    print 'x:   {0:.2f}   ->  {1:.2f}'.format(x_core, x_core_fit)
    print 'y:   {0:.2f}   ->  {1:.2f}'.format(y_core, y_core_fit)

    #print 'Ne:  {0:.2f}   ->  {1:.2f}'.format(shower_size/f_den, Ne_fit)







    R_Max=0
    R_Min=1000000
    radius_show=np.zeros([LORA.nDetA])
    radius_bin_show=np.zeros([LORA.nDetA])

    for i in np.arange(LORA.nDetA):
        if rho[i]>=LORA.Density_Cut and rho[i]<=LORA.Density_Cut_High:
            radius_show[i]=np.abs(np.sqrt(np.power(x_core_fit-x_shower[i],2)+np.power(y_core_fit-y_shower[i],2)))
            radius_bin_show[i]=lat_den_show.FindBin(radius_show[i])
            if  radius_show[i]>R_Max:
                R_Max=radius_show[i]
            if radius_show[i]<R_Min:
                R_Min=radius_show[i]
            lat_den_show.SetBinContent(int(radius_bin_show[i]),rho[i])
            lat_den_show.SetBinError(int(radius_bin_show[i]),rho_err[i])



    Ne_fit,rM_fit,s_fit,Ne_fit_er=TF1_Fit_NKG(lat_den_show,Ne_fit,LORA.rM,LORA.Age,R_Min,R_Max)

    print 'Ne:  {0:.2f},   Rm:  {1:.2f}, xCore:  {2:.2f},  yCore:  {3:.2f}'.format(Ne_fit, rM_fit,x_core_fit, y_core_fit)



    # iterate 2d, 1d fits
    for k in np.arange(4):
    
        del lat_den_show
        lat_den_show=TH1F('Lateral_density','Lateral_density',LORA.No_Bin_R,LORA.min_R,LORA.max_R)
        Ne_fit,rM_fit,s_fit,x_core_fit,y_core_fit,x_core_fit_err,y_core_fit_err,corr_coef_xy=TF2_Fit_NKG(evt_display_shower,Ne_fit,rM_fit,s_fit,x_core_fit,y_core_fit,x_min_shower,x_max_shower,y_min_shower,y_max_shower)
        
        R_Max=0
        R_Min=1000000
  
    
        for i in np.arange(LORA.nDetA):
            if rho[i]>=LORA.Density_Cut and rho[i]<=LORA.Density_Cut_High:
                radius_show[i]=np.abs(np.sqrt(np.power(x_core_fit-x_shower[i],2)+np.power(y_core_fit-y_shower[i],2)))
                radius_bin_show[i]=lat_den_show.FindBin(radius_show[i])
            if  radius_show[i]>R_Max:
                R_Max=radius_show[i]
            if radius_show[i]<R_Min:
                R_Min=radius_show[i]
            lat_den_show.SetBinContent(int(radius_bin_show[i]),rho[i])
            lat_den_show.SetBinError(int(radius_bin_show[i]),rho_err[i])

        Ne_fit,rM_fit,s_fit,Ne_fit_err=TF1_Fit_NKG(lat_den_show,Ne_fit,rM_fit,s_fit,R_Min,R_Max)

        print 'Ne:  {0:.2f},   Rm:  {1:.2f}, xCore:  {2:.2f},  yCore:  {3:.2f}'.format(Ne_fit, rM_fit,x_core_fit,y_core_fit)

    

    # do atm correction ##################################

    atm_data=read_attenuation()
    size_theta=np.zeros([30])
    if len(atm_data)!=11:
        print 'no atm corrections'
    else:
        f_no=len(atm_data)
        f_int=atm_data.T[0]
        f_logN_Ref=atm_data.T[1]
        err1=atm_data.T[2]
        f_X0=atm_data.T[3]
        err2=atm_data.T[4]
        f_lamb=atm_data.T[5]
        err3=atm_data.T[6]
        Lambda0=0
        err_Lambda0=0
        for k in np.arange(len(f_int)):
            size_theta[k]=f_logN_Ref[k]-(f_X0[k]/f_lamb[k])*(1/np.cos(theta)-1/np.cos(np.pi*LORA.Ref_angle/180.0))*0.4342944819 #log10(size) for attenuation curve kk at zenith angle 'theta'
            #print size_theta[k]
            
            if np.log10(Ne_fit) >= size_theta[1]: ## typo?   1->k?
                Lambda0=f_lamb[k]    #  //Extropolation
                err_Lambda0=err3[k]
                print 'Extrapolation: k={0}  s_theta={1}  s2={2}  lamb={3}'.format(k,size_theta[k],np.log10(Ne_fit),Lambda0)
                break
            
            elif np.log10(Ne_fit) >= size_theta[k]:
                Lambda0=(f_lamb[k]*(np.log10(Ne_fit)-size_theta[k-1])+f_lamb[k-1]*(size_theta[k]-np.log10(Ne_fit)))/(size_theta[k]-size_theta[k-1]) #  //Interpolation
                err_Lambda0=(err3[k]*(np.log10(Ne_fit)-size_theta[k-1])+err3[k-1]*(size_theta[k]-np.log10(Ne_fit)))/(size_theta[k]-size_theta[k-1]) #  //Interpolation
                print 'Interpolation: k={0}  s_theta={1}  s2={2}  lamb={3}'.format(k,size_theta[k],np.log10(Ne_fit),Lambda0)
                break
    ####

    size_theta[f_no-1]=f_logN_Ref[f_no-1]-(f_X0[f_no-1]/f_lamb[f_no-1])*(1/np.cos(theta)-1/np.cos(np.pi*LORA.Ref_angle/180))*0.4342944819

    if np.log10(Ne_fit) >= size_theta[f_no-1]:
        Lambda0=f_lamb[f_no-1] #; //Extropolation
        err_Lambda0=err3[f_no-1]
        print 'Extrapolation: k={0}  s_theta={1}  s2={2}  lamb={3}'.format(k,size_theta[k],np.log10(Ne_fit),Lambda0)


    log_size_Ref=np.log10(Ne_fit)+(LORA.X0/Lambda0)*(1/np.cos(theta)-1/np.cos(np.pi*LORA.Ref_angle/180.0))*0.4342944819# ; //log10()
    size_Ref=np.power(10,log_size_Ref)
    err_size_Ref=np.sqrt(np.power(Ne_fit_err/Ne_fit,2)+np.power(np.log(10)*LORA.X0*(1/np.cos(theta)-1/np.cos(np.pi*LORA.Ref_angle/180.0))*0.4342944819,2)*np.power(err_Lambda0/np.power(Lambda0,2),2))*size_Ref

    #energy_Ref=pow(size_Ref,par_b)*pow(10,par_a)*pow(10,-6) ; //Energy(PeV) at Ref_angle: Formula from KASCADE simulation (2008)
    #err_energy_Ref=sqrt(pow(log(10)*err_a,2)+pow(log(size_Ref)*err_b,2)+pow(par_b*err_size_Ref/size_Ref,2))*energy_Ref ;    //error on energy at Ref_angle
    energy_Ref=np.power(size_Ref,LORA.par_b)*np.power(10.0,LORA.par_a)*np.power(10.0,-6.0) #; #//Energy(PeV) at Ref_angle: Formula from KASCADE simulation (2008)
    err_energy_Ref=np.sqrt(np.power(np.log(10)*LORA.err_a,2)+np.power(np.log10(size_Ref)*LORA.err_b,2)+np.power(LORA.par_b*err_size_Ref/size_Ref,2))*energy_Ref ##;    //error on energy at
    
    



    #undo core rotation
    X3,Y3,Z3=back_theta_phi(theta,phi,psi,x_core_fit,y_core_fit,0)
    l,m,n=back_theta_phi(theta,phi,psi,0,0,1)

    x_core_ground=l*(-Z3/n)+X3 ;
    y_core_ground=m*(-Z3/n)+Y3 ;
    print 'x_core={0:.2f}   y_core={1:.2f}'.format(x_core_ground,y_core_ground)

    energy=np.power(Ne_fit,LORA.par_b)*np.power(10.0,LORA.par_a)*np.power(10.0,-6.0)# ; //Energy (PeV): Formula from KASCADE simulation (2008)
    err_energy=np.sqrt(np.power(np.log(10)*LORA.err_a,2)+np.power(np.log10(Ne_fit)*LORA.err_b,2)+np.power(LORA.par_b*Ne_fit_err/Ne_fit,2))*energy #;    //error on energy

    
    #find first time stamp
    min_utc=1e10
    min_nsec=1e10
    for i in np.arange(LORA.nDetA):
        if detectors[i].gps>1 and detectors[i].cal_time>1:
            if detectors[i].gps<min_utc:
                min_utc=detectors[i].gps
            if detectors[i].cal_time<min_nsec:
                min_nsec=detectors[i].cal_time

    # assign event values
    event.x_core=x_core_ground
    event.y_core=y_core_ground
    event.x_core_err=x_core_fit_err
    event.y_core_err=y_core_fit_err
    event.z_core=0
    event.UTC_min=min_utc
    event.nsec_min=min_nsec
    event.energy= energy
    event.energy_err=err_energy
    event.Rm=rM_fit
    event.Ne=Ne_fit
    event.Ne_err=Ne_fit_err
    event.CorCoef_xy=corr_coef_xy
    event.Ne_RefA=size_Ref
    event.NeErr_RefA=err_size_Ref
    event.Energy_RefA=energy_Ref
    event.EnergyErr_RefA=err_energy_Ref

