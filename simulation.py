'''
Created on Sep 30, 2015
'''


import cProfile
import profile
import time
import sys


import matplotlib.pyplot as plt
import numpy as np

from Configuration import Configuration
from MotorUnitPool import MotorUnitPool
from InterneuronPool import InterneuronPool
from NeuralTract import NeuralTract
from SynapsesFactory import SynapsesFactory
from jointAnkleForceTask import jointAnkleForceTask

def simulator():

    conf = Configuration('confTest.rmto')

    #idx = np.where(conf.confArray['f0']=='MUnumber_SOL-S')[0][0]
    #conf.confArray['f1'][idx] = nMN

    pools = dict()
    pools[0] = MotorUnitPool(conf, 'SOL')
    pools[1] = NeuralTract(conf, 'CMExt')
    #pools[2] = InterneuronPool(conf, 'RC', 'ext')

    #ankle = jointAnkleForceTask(conf, pools)
    Syn = SynapsesFactory(conf, pools)
    del Syn

    #conf.simDuration_ms = duration
    t = np.arange(0.0, conf.simDuration_ms, conf.timeStep_ms)

    tic = time.time()
    for i in xrange(0, len(t)):
        pools[1].atualizePool(t[i])
        pools[0].atualizeMotorUnitPool(t[i])
        #pools[3].atualizePool(t[i]) # RC synaptic Noise
        #pools[2].atualizeInterneuronPool(t[i]) # RC pool
    toc = time.time()
    print str(toc - tic) + ' seconds'

    pools[0].listSpikes()
    pools[1].listSpikes()
    '''
    plt.figure()
    plt.plot(pools[1].poolTerminalSpikes[:, 0],
             pools[1].poolTerminalSpikes[:, 1]+1, '.')
    
    
    plt.figure()
    plt.plot(pools[0].poolSomaSpikes[:, 0],
             pools[0].poolSomaSpikes[:, 1]+1, '.')

    plt.figure()
    plt.plot(t, pools[0].Muscle.force, '-')

    plt.figure()
    plt.plot(t, dendV, '-')

    plt.figure()
    plt.plot(t, somaV, '-')
    '''
if __name__ == '__main__':

    #cProfile.run('simulator()', sort = 'tottime')
    
    #np.__config__.show()
    
    
    simulator()
    '''
    plt.show()
    '''
