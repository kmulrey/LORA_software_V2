import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

nDet=4
dtV1=2.5

plt.ion()


def plot_traces(info):

    time=np.arange(0,len(info[0]['counts'])*dtV1,dtV1)
    
    fig = plt.figure()

    ax1 = fig.add_subplot(5,1,1)
    ax2 = fig.add_subplot(5,1,2)
    ax3 = fig.add_subplot(5,1,3)
    ax4 = fig.add_subplot(5,1,4)
    ax5 = fig.add_subplot(5,1,5)

    for i in np.arange(nDet):
        ax1.plot(time,info[0+i]['counts'])
        ax2.plot(time,info[4+i]['counts'])
        ax3.plot(time,info[8+i]['counts'])
        ax4.plot(time,info[12+i]['counts'])
        ax5.plot(time,info[16+i]['counts'])

    
    '''
    ax1.axis([-0.6,0.6,1.75,3])
    ax2.axis([-0.6,0.6,1.75,3])
    ax3.axis([-0.6,0.6,1.75,3])
    ax4.axis([-0.6,0.6,1.75,3])

    ax1.set_ylabel('log(q/Me)')
    ax2.set_ylabel('log(q/Me)')
    ax3.set_ylabel('log(q/Me)')
    ax4.set_ylabel('log(q/Me)')
    ax4.set_xlabel('x (m)')
        '''


    plt.show()
    raw_input()
    plt.close()
