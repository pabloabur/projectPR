# VAN KEULEN 1981
# Recording obtained from a pair of connected RC and MN
# RC recorded extracellularly and MN intracellularly. 
# RC excited through step currents with 0.5ms duration and
# through glutamate injection.
# Dorsal roots cut in the anaesthesized cat spinal cord.

# ROSS ET AL. 1975
# Same seeting, but with decerebrated cat with dorsal and ventral roots cut.
# Current chosen to be just above threshlod. Here RC bursts are less regularly
# and repetitively than in van keulen.

# TODO problems: 
# Cannot simulate in glutamate injection on RC described on the experiment and 
# extracellular recording

import time


import matplotlib.pyplot as plt
import numpy as np

from Configuration import Configuration
from MotorUnitPoolOpt import MotorUnitPool
from InterneuronPoolOpt import InterneuronPool
from NeuralTract import NeuralTract
from SynapsesFactory import SynapsesFactory

def simulator():

    conf = Configuration('confuchiyama.rmto')

    pools = dict()
    pools[0] = MotorUnitPool(conf, 'SOL')
    pools[1] = InterneuronPool(conf, 'RC', 'ext')

    Syn = SynapsesFactory(conf, pools)

    t = np.arange(0.0, conf.simDuration_ms, conf.timeStep_ms)

    RC_mV = np.zeros_like(t)
    MN_mV = np.zeros_like(t)

    tic = time.clock()
    for i in xrange(0, len(t)):
        # t=step*i [ms]
        # Current in nA, adjusted to reproduce approximately the experimental 3mV
        if t[i]>10 and t[i]<510:
            for j in xrange(len(pools[0].unit)):
                pools[0].iInjected[2*j+1] = 13
        else:
            for j in xrange(len(pools[0].unit)):
                pools[0].iInjected[2*j+1] = 0
        pools[0].atualizeMotorUnitPool(t[i]) # MN pool
        pools[2].atualizePool(t[i]) # RC synaptic Noise
        pools[1].atualizeInterneuronPool(t[i]) # RC pool
        RC_mV[i] = pools[1].unit[0].v_mV[0] 
        MN_mV[i] = pools[0].unit[0].v_mV[1] 
    toc = time.clock()
    print str(toc - tic) + ' seconds'

    pools[0].listSpikes()
    pools[1].listSpikes()
 
    plt.figure()
    plt.plot(pools[1].poolSomaSpikes[:, 0],
             pools[1].poolSomaSpikes[:, 1]+1, '.')
    
    
    plt.figure()
    plt.plot(pools[0].poolSomaSpikes[:, 0],
             pools[0].poolSomaSpikes[:, 1]+1, '.')

    plt.figure()
    plt.plot(t, pools[0].Muscle.force, '-')

    plt.figure()
    plt.plot(t, RC_mV, '-')

    plt.figure()
    plt.plot(t, MN_mV, '-')
    
if __name__ == '__main__':

    np.__config__.show()
    simulator()
    plt.show()
