# CLEVELAND ET AL. 1981
# Experiment in decerebrate cats with ipsilateral dorsal and ventral roots cut.
# Renshaw cell recorded extracellularly and aMN impaled, depolarized
# intracellularly, and then subjected to antidromic conditioning stimulation of
# various frequencies. Stimulus pulses were 200mus wide and supramaximal for
# a-axons.The period of discharge adaptation on the onset of stimulus must have
# come to an end.

import time


import matplotlib.pyplot as plt
import numpy as np

from Configuration import Configuration
from MotorUnitPool import MotorUnitPool
from InterneuronPool import InterneuronPool
from NeuralTract import NeuralTract
from SynapsesFactory import SynapsesFactory

def simulator():

    conf = Configuration('confnlinIO.rmto')

    pools = dict()
    pools[0] = MotorUnitPool(conf, 'SOL')
    pools[1] = InterneuronPool(conf, 'RC', 'ext')

    for i in xrange(0,len(pools[0].unit)):
        pools[0].unit[i].createStimulus()

    Syn = SynapsesFactory(conf, pools)

    t = np.arange(0.0, conf.simDuration_ms, conf.timeStep_ms)

    RC_mV = np.zeros_like(t)
    MN_mV = np.zeros_like(t)

    tic = time.clock()
    for i in xrange(0, len(t)):
        pools[0].atualizeMotorUnitPool(t[i]) # MN pool
        pools[2].atualizePool(t[i]) # RC synaptic Noise
        pools[1].atualizeInterneuronPool(t[i]) # RC pool
        RC_mV[i] = pools[1].v_mV[0] 
        MN_mV[i] = pools[0].v_mV[1] 
    toc = time.clock()
    print str(toc - tic) + ' seconds'

    pools[0].listSpikes()

    plt.figure()
    plt.plot(t, pools[0].unit[0].nerveStimulus_mA)
    
    plt.figure()
    plt.plot(t, MN_mV, '-')
    
    #plt.figure()
    #plt.plot(pools[0].poolSomaSpikes[:, 0],
    #    pools[0].poolSomaSpikes[:, 1]+1, '.')

    plt.figure()
    plt.plot(t, RC_mV, '-')

    #plt.figure()
    #plt.plot(t, pools[0].Muscle.force, '-')

if __name__ == '__main__':

    np.__config__.show()
    simulator()
    plt.show()
