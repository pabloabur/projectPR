import matplotlib.pyplot as plt
import numpy as np

from Configuration import Configuration
from MotorUnitPool import MotorUnitPool
from InterneuronPool import InterneuronPool
from NeuralTract import NeuralTract
from SynapsesFactory import SynapsesFactory

def simulator1():

    conf = Configuration('confGH.rmto')

    pools = dict()
    pools[0] = MotorUnitPool(conf, 'SOL')
    pools[1] = InterneuronPool(conf, 'RC', 'ext')
    pools[2] = NeuralTract(conf, 'CMExt')

    Syn = SynapsesFactory(conf, pools)

    t = np.arange(0.0, conf.simDuration_ms, conf.timeStep_ms)

    RC_mV = np.zeros_like(t)
    MN_mV = np.zeros_like(t)

    for i in xrange(0, len(t)):
        pools[2].atualizePool(t[i]) # Neural tract
        pools[0].atualizeMotorUnitPool(t[i]) # MN pool
        pools[3].atualizePool(t[i]) # RC synaptic Noise
        pools[1].atualizeInterneuronPool(t[i]) # RC pool
        RC_mV[i] = pools[1].unit[0].v_mV[0] 
        MN_mV[i] = pools[0].unit[0].v_mV[1] 

    pools[0].listSpikes()
    pools[1].listSpikes()
    pools[2].listSpikes()

    # Calculate feedback gain
    ## Calculate number of spikes
    numRecruitedUnits = int(max(pools[0].poolSomaSpikes[:, 1]))
    numberSpikes = [None] * numRecruitedUnits
    firingRate = [None] * numRecruitedUnits
    for i in xrange(numRecruitedUnits):
        # Count number of spikes for each recruited unit
        numberSpikes[i] = pools[0].poolSomaSpikes[:,1].tolist().count(i)
        firingRate[i] = numberSpikes[i]/conf.simDuration_ms
    ## Calcultate average MN activity
    activity = np.mean(firingRate) * numRecruitedUnits

    plt.figure()
    plt.plot(pools[0].poolSomaSpikes[:, 0],
        pools[0].poolSomaSpikes[:, 1]+1, '.')

    plt.figure()
    plt.plot(t, MN_mV)

    print '{} recruited units, with average activity of {}'.format(numRecruitedUnits,activity)

def simulator2():

    conf = Configuration('confGH.rmto')

    pools = dict()
    pools[0] = MotorUnitPool(conf, 'SOL')
    pools[1] = NeuralTract(conf, 'CMExt')

    Syn = SynapsesFactory(conf, pools)

    t = np.arange(0.0, conf.simDuration_ms, conf.timeStep_ms)

    MN_mV = np.zeros_like(t)

    for i in xrange(0, len(t)):
        pools[1].atualizePool(t[i]) # Neural tract
        pools[0].atualizeMotorUnitPool(t[i]) # MN pool
        MN_mV[i] = pools[0].unit[0].v_mV[1] 

    pools[0].listSpikes()
    pools[1].listSpikes()

    # Calculate feedback gain
    ## Calculate number of spikes
    numRecruitedUnits = int(max(pools[0].poolSomaSpikes[:, 1]))
    numberSpikes = [None] * numRecruitedUnits
    firingRate = [None] * numRecruitedUnits
    for i in xrange(numRecruitedUnits):
        # Count number of spikes for each recruited unit
        numberSpikes[i] = pools[0].poolSomaSpikes[:,1].tolist().count(i)
        firingRate[i] = numberSpikes[i]/conf.simDuration_ms
    ## Calcultate average MN activity
    activity = np.mean(firingRate) * numRecruitedUnits

    plt.figure()
    plt.plot(pools[0].poolSomaSpikes[:, 0],
        pools[0].poolSomaSpikes[:, 1]+1, '.')

    plt.figure()
    plt.plot(t, MN_mV)

    print '{} recruited units, with average activity of {}'.format(numRecruitedUnits,activity)

if __name__ == '__main__':

    np.__config__.show()
    simulator1()
    simulator2()
    plt.show()

