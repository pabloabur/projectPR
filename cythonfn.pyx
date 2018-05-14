#cython: boundscheck=False

import numpy as np
cimport numpy as np

def runge_kutta(float t, double[:] x, float timeStep, float timeStepByTwo,
        float timeStepBySix, unsigned int MUnumber, unit, double[:] iIonic,
        G, double[:] iInjected,
        double[:] EqCurrent_nA, double[:] capacitanceInv):
    k1 = dVdt(t, x, MUnumber, unit, iIonic, G, iInjected, EqCurrent_nA, capacitanceInv)
    k2 = dVdt(t + timeStepByTwo, x + timeStepByTwo * k1, MUnumber, unit, iIonic, G, iInjected, EqCurrent_nA, capacitanceInv)
    k3 = dVdt(t + timeStepByTwo, x + timeStepByTwo * k2, MUnumber, unit, iIonic, G, iInjected, EqCurrent_nA, capacitanceInv)
    k4 = dVdt(t + timeStep, x + timeStep * k3, MUnumber, unit, iIonic, G, iInjected, EqCurrent_nA, capacitanceInv)
    
    #xout = x + timeStepBySix * (np.add(
    #                                   np.add(
    #                                         np.add(k1, k2, order = 'C'), np.add(k2, k3, order='C')
    #                                         ),
    #                                   np.add(k3, k4, order='C'), order='C')
    #                           )
    xout = x + timeStepBySix * k1 + k2 + k2 + k3 + k3 + k4
    return xout

'''np.ndarray[np.float, ndim=2] G'''

def dVdt(float t, double[:] V, unsigned int MUnumber, unit, double[:] iIonic,
        G, double[:] iInjected,
        double[:] EqCurrent_nA, double[:] capacitanceInv):
    #k = 0
    cdef unsigned int i, j
    cdef int[:] compNumber = np.empty(len(unit), dtype=np.int32)

    for i in xrange(len(unit)):
        compNumber[i] = unit[i].compNumber

    for i in xrange(MUnumber):
        for j in xrange(compNumber[i]):
            iIonic[i*compNumber[0]+j] = computeCurrent(t,
                                        V[i*compNumber[0]+j],
                                        unit[i].compartment[j])
            #k += 1
          
    res = G.dot(V)
    Vout = (iIonic + res + iInjected + EqCurrent_nA) * capacitanceInv
    return Vout

def computeCurrent(self, t, V_mV, compartment):
    I = 0.0

    if comparment.SynapsesIn[0].numberOfIncomingSynapses:
        I += computeSynapticCurrent(t, V_mV, comparment.SynapsesIn[0])
    if comparment.SynapsesIn[1].numberOfIncomingSynapses:
        I += computeSynapticCurrent(t, V_mV, comparment.SynapsesIn[1])
    for i in xrange(0, self.numberChannels):
        I += self.Channels[i].computeChannelCurrent(t, V_mV)
    
    return I

def computeSynapticCurrent(self, t, V_mV, synapse):
    if len(self.tEndOfPulse) == 0:
        self.tBeginOfPulse = np.ones_like(self.gmax_muS,
                                          dtype=float) * float("-inf")
        self.tEndOfPulse = np.ones_like(self.gmax_muS,
                                        dtype=float) * float("-inf")
        self.tLastPulse = np.ones_like(self.gmax_muS,
                                       dtype=float) * float("-inf")
        self.conductanceState = np.zeros_like(self.gmax_muS,
                                              dtype=int)
        self.ri = np.zeros_like(self.gmax_muS, dtype=float)
        self.ti = np.zeros_like(self.gmax_muS, dtype=float)
        self.dynamicGmax = np.zeros_like(self.gmax_muS, dtype=float)
        self.synContrib = self.gmax_muS / self.gMaxTot_muS
        self.computeCurrent = self.computeCurrent2
    
    return self.computeConductance(t) * (self.EqPot_mV - V_mV)
