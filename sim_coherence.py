from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

from Configuration import Configuration
from MotorUnitPool import MotorUnitPool
from InterneuronPool import InterneuronPool
from NeuralTract import NeuralTract
from SynapsesFactory import SynapsesFactory

def simulator():

    conf = Configuration('confMNRC.rmto')

    pools = dict()
    pools[0] = MotorUnitPool(conf, 'SOL')
    pools[1] = InterneuronPool(conf, 'RC', 'ext')

    Syn = SynapsesFactory(conf, pools)

    t = np.arange(0.0, conf.simDuration_ms, conf.timeStep_ms)

    RC_mV = np.zeros_like(t)
    MN_mV = np.zeros_like(t)

    for i in xrange(0, len(t)):
        # t=step*i [ms]
        # Current in nA, adjusted to reproduce approximately the experimental 3mV
        pools[0].iInjected[1] = 13
        pools[0].atualizeMotorUnitPool(t[i]) # MN pool
        pools[2].atualizePool(t[i]) # RC synaptic Noise
        pools[1].atualizeInterneuronPool(t[i]) # RC pool
        RC_mV[i] = pools[1].unit[0].v_mV[0] 
        MN_mV[i] = pools[0].unit[0].v_mV[1] 

    fs=1/(conf.timeStep_ms*1e-3)
    window = 1000
    f, Cxy = signal.coherence(MN_mV, RC_mV, fs, nperseg=window)

    plt.figure()
    plt.plot(f, Cxy)
    plt.xlabel('f (Hz)')
    plt.ylabel('Coherence')
    plt.xlim((0, 500))
    
if __name__ == '__main__':

    np.__config__.show()
    simulator()
    plt.show()
