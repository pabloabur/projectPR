#!/usr/bin/env /home/renato/intel/intelpython27/bin/python
import sys
sys.path.insert(0, '..')
import time
import matplotlib.pyplot as plt
#matplotlib inline  
import numpy as np

from Configuration import Configuration
from MotorUnitPoolMPI import MotorUnitPool
from InterneuronPool import InterneuronPool
from NeuralTract import NeuralTract
from SynapsesFactory import SynapsesFactory

conf = Configuration('confTest.rmto')

# Time vector for the simulation
t = np.arange(0.0, conf.simDuration_ms, conf.timeStep_ms)

membPotential = np.zeros_like(t, dtype = 'd')
dendV = np.zeros_like(t)
somaV = np.zeros_like(t)

pools = dict()
pools[0] = MotorUnitPool(conf, 'SOL')
Syn = SynapsesFactory(conf, pools)

tic = time.clock()
for i in xrange(0, len(t)-1):
    for j in xrange(len(pools[0].unit)):
        pools[0].unit[j].iInjected[1] = 10
    pools[0].atualizeMotorUnitPool(t[i]) # MN pool
    dendV[i] = pools[0].unit[2].v_mV[0]
    somaV[i] = pools[0].unit[2].v_mV[1]
toc = time.clock()
print str(toc - tic) + ' seconds', pools[0].rank
#print str(toc - tic) + ' seconds'

pools[0].listSpikes()

#if pools[0].rank == 1:
if True:
    plt.figure()
    plt.plot(pools[0].poolSomaSpikes[:, 0],
                pools[0].poolSomaSpikes[:, 1]+1, '.')
    #plt.figure()
    #plt.plot(pools[0].poolTerminalSpikes[:, 0],
    #            pools[0].poolTerminalSpikes[:, 1]+1, '.')         
    
    plt.figure()
    plt.plot(t, pools[0].Muscle.force, '-')

    #print 'M = ' + str(np.mean(pools[0].Muscle.force[int(1000/conf.timeStep_ms):-1]))
    #print 'SD = ' + str(np.std(pools[0].Muscle.force[int(1000/conf.timeStep_ms):-1]))

    plt.figure()
    plt.plot(t, dendV, '-')

    plt.figure()
    plt.plot(t, somaV, '-')

    plt.show()
