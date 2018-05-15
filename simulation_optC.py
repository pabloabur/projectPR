'''
Created on Sep 30, 2015

@author: root
'''


import cProfile
import profile
import time
import sys

import matplotlib.pyplot as plt
import numpy as np

from Configuration import Configuration
from CythonMotorUnitPoolOpt import MotorUnitPool
from InterneuronPoolOpt import InterneuronPool
from NeuralTract import NeuralTract
from SynapsesFactory import SynapsesFactory

def simulator():

    conf = Configuration('confTest.rmto')

    pools = dict()
    pools[0] = MotorUnitPool(conf, 'SOL')
    #pools.append(InterneuronPool(conf, 'RC'))

    Syn = SynapsesFactory(conf, pools)
    del Syn
    
    t = np.arange(0.0, conf.simDuration_ms, conf.timeStep_ms)

    dendV = np.zeros_like(t)
    somaV = np.zeros_like(t)
    internodeV = np.zeros_like(t)
    nodeV = np.zeros_like(t)

    tic = time.time()
    for i in xrange(0, len(t)):
        for j in xrange(len(pools[0].unit)):
            pools[0].iInjected[j] = 15
        pools[0].atualizeMotorUnitPool(t[i])
        dendV[i] = pools[0].unit[2].v_mV[0]
        somaV[i] = pools[0].unit[2].v_mV[1] 
    toc = time.time()
    print str(toc - tic) + ' seconds'

    pools[0].listSpikes()
   
    '''
    #plt.figure()
    #plt.plot(pools[1].poolTerminalSpikes[:, 0],
    #         pools[1].poolTerminalSpikes[:, 1]+1, '.')
    
    
    plt.figure()
    plt.plot(pools[0].poolSomaSpikes[:, 0],
             pools[0].poolSomaSpikes[:, 1]+1, '.')

    plt.figure()
    plt.plot(t, pools[0].Muscle.force, '-')

    plt.figure()
    plt.plot(t, dendV, '-')

    plt.figure()
    plt.plot(t, somaV, '-')

    plt.show()
    '''
    
if __name__ == '__main__':
    #cProfile.run('simulator()', sort = 'tottime')
    #np.__config__.show()
    
    simulator()
