import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import LORAparameters as LORA
import process_functions as process


plt.ion()


def plot_traces(detectors):

    time=np.arange(0,len(detectors[0].counts)*LORA.dtV1,LORA.dtV1)
    
    
    peaks, peak_args=process.find_peaks(detectors)
    
    first_peak=peak_args[np.argmax(peaks)]
    
    fig = plt.figure(figsize=(5,6))

    ax1 = fig.add_subplot(5,1,1)
    ax2 = fig.add_subplot(5,1,2)
    ax3 = fig.add_subplot(5,1,3)
    ax4 = fig.add_subplot(5,1,4)
    ax5 = fig.add_subplot(5,1,5)
    
    ax1.set_xticks([], [])
    ax2.set_xticks([], [])
    ax3.set_xticks([], [])
    ax4.set_xticks([], [])

    #x=time[int(detectors[0].threshold_time/2.5)]
    #y= detectors[0].counts[int(detectors[0].threshold_time/2.5)],
    for i in np.arange(LORA.nDet):
        ax1.plot(time,detectors[0+i].counts)
        #ax1.plot(time[int(detectors[0+i].threshold_time/2.5)],detectors[0+i].counts[int(detectors[0+i].threshold_time/2.5)],'o')
        ax2.plot(time,detectors[4+i].counts)
        ax3.plot(time,detectors[8+i].counts)
        ax4.plot(time,detectors[12+i].counts)
        ax5.plot(time,detectors[16+i].counts)
        
    

    
    #ax1.plot(x,y,'o')

    
    ax1.set_ylabel('ADC')
    ax2.set_ylabel('ADC')
    ax3.set_ylabel('ADC')
    ax4.set_ylabel('ADC')
    ax5.set_ylabel('ADC')

    ax5.set_xlabel('time (ns)')


    pos_window=1500
    neg_window=500
    '''
    ax1.set_xlim(time[int(first_peak)]-neg_window,time[int(first_peak)]+pos_window)
    ax2.set_xlim(time[int(first_peak)]-neg_window,time[int(first_peak)]+pos_window)
    ax3.set_xlim(time[int(first_peak)]-neg_window,time[int(first_peak)]+pos_window)
    ax4.set_xlim(time[int(first_peak)]-neg_window,time[int(first_peak)]+pos_window)
    ax5.set_xlim(time[int(first_peak)]-neg_window,time[int(first_peak)]+pos_window)
    '''

    #ax2.axis([-0.6,0.6,1.75,3])
    #    ax3.axis([-0.6,0.6,1.75,3])
    #ax4.axis([-0.6,0.6,1.75,3])



    fig.tight_layout()

    plt.show()
    raw_input()
    plt.close()
    
