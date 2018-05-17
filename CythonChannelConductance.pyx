import numpy as np
import math
from CythonPulseConductanceState import PulseConductanceState

cdef double compCondKs(double V_mV, double gmax, list state, double EqPot):
        return gmax * (state[0].value * state[0].value) * (EqPot - V_mV)


cdef double compCondKsaxon(double V_mV, double gmax, list state, double EqPot):
        return gmax * state[0].value * (EqPot - V_mV)

cdef double compCondNa(double V_mV, double gmax, list state, double EqPot):
        return gmax * (state[0].value * state[0].value * state[0].value) * state[1].value * (EqPot - V_mV)

cdef double compCondNap(double V_mV, double gmax, list state, double EqPot):
        return gmax * (state[0].value * state[0].value * state[0].value) * (EqPot - V_mV)

cdef double compCondKf(double V_mV, double gmax, list state, double EqPot):
        return gmax * (state[0].value * state[0].value * state[0].value * state[0].value) * (EqPot - V_mV)

cdef double compCondH(double V_mV, double gmax, list state, double EqPot):
        return gmax * state[0].value * (EqPot - V_mV)

cdef class ChannelConductance:
    cdef str kind, neuronKind, compKind, stateType, pool
    cdef double compArea, EqPot_mV, gmax_muS
    cdef unsigned int index, lenStates, i
    cdef list condState

    def __init__(self, str kind, conf, double compArea, str pool, str neuronKind, str compKind, unsigned int index):
        ## string with the type of the ionic channel. For now it 
        ## can be *Na* (Sodium), *Ks* (slow Potassium), *Kf* (fast Potassium) or 
        ## *Ca* (Calcium).
        # TODO maybe this is not necessary
        #self.kind = <str>kind
        
        ## Equilibrium Potential of the ionic channel, mV.
        self.EqPot_mV = <double>conf.parameterSet('EqPot_' + kind + '@' + compKind, pool, index)
        ## Maximal conductance, in \f$\mu\f$S, of the ionic channel. 
        self.gmax_muS = compArea * <double>conf.parameterSet('gmax_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool, index)
                        
        ## String with type of dynamics of the states. For now it accepts the string pulse.
        self.stateType = conf.parameterSet('StateType', pool, index)
        
        if self.stateType == 'pulse':
            ConductanceState = PulseConductanceState

        ## List of ConductanceState objects, representing each state of the ionic channel.
        self.condState = []
        
        if self.kind == 'Kf':
            self.condState.append(ConductanceState('n', conf, pool, neuronKind, compKind, index))
#            beta = float(conf.parameterSet('beta_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool, index))
#            alpha = float(conf.parameterSet('alpha_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool,index))
#            PulseDuration = float(conf.parameterSet('PulseDur_' + kind, pool, index))
#            self.condState.append(ConductanceState('n', beta, alpha, PulseDuration, conf.timeStep_ms))
            ## Function that computes the conductance dynamics.
            self.compCond = compCondKf
        if self.kind == 'Ks':
#            beta = float(conf.parameterSet('beta_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool, index))
#            alpha = float(conf.parameterSet('alpha_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool,index))
#            PulseDuration = float(conf.parameterSet('PulseDur_' + kind, pool, index))
#            self.condState.append(ConductanceState('q', beta, alpha, PulseDuration, conf.timeStep_ms))
            self.condState.append(ConductanceState('q', conf, pool, neuronKind, compKind, index))
            self.compCond = compCondKs
        if self.kind == 'Na':
#            beta = float(conf.parameterSet('beta_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool, index))
#            alpha = float(conf.parameterSet('alpha_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool,index))
#            PulseDuration = float(conf.parameterSet('PulseDur_' + kind, pool, index))
#            self.condState.append(ConductanceState('m', beta, alpha, PulseDuration, conf.timeStep_ms))
#            self.condState.append(ConductanceState('h', beta, alpha, PulseDuration, conf.timeStep_ms))
            self.condState.append(ConductanceState('m', conf, pool, neuronKind, compKind, index))
            self.condState.append(ConductanceState('h', conf, pool, neuronKind, compKind, index))
            self.compCond = compCondNa
        if self.kind == 'Ca':
            pass  # to be implemented
        if self.kind == 'Nap':
#            beta = float(conf.parameterSet('beta_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool, index))
#            alpha = float(conf.parameterSet('alpha_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool,index))
#            PulseDuration = float(conf.parameterSet('PulseDur_' + kind, pool, index))
#            self.condState.append(ConductanceState('mp', beta, alpha, PulseDuration, conf.timeStep_ms))
            self.condState.append(ConductanceState('mp', conf, pool, neuronKind, compKind, index))
            self.compCond = compCondNap
        if self.kind == 'KsAxon':
#            beta = float(conf.parameterSet('beta_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool, index))
#            alpha = float(conf.parameterSet('alpha_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool,index))
#            PulseDuration = float(conf.parameterSet('PulseDur_' + kind, pool, index))
#            self.condState.append(ConductanceState('s', beta, alpha, PulseDuration, conf.timeStep_ms))
            self.condState.append(ConductanceState('s', conf, pool, neuronKind, compKind, index))
            self.compCond = compCondKsaxon
        if self.kind == 'H':
#            beta = float(conf.parameterSet('beta_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool, index))
#            alpha = float(conf.parameterSet('alpha_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool,index))
#            PulseDuration = float(conf.parameterSet('PulseDur_' + kind, pool, index))
#            self.condState.append(ConductanceState('qh', beta, alpha, PulseDuration, conf.timeStep_ms))
            self.condState.append(ConductanceState('qh', conf, pool, neuronKind, compKind, index))
            self.compCond = compCondH            
        
        ## Integer with the number of states in the ionic channel.    
        self.lenStates = len(self.condState)          
    
    cdef double computeCurrent(self, double t, double V_mV): 
        '''
        Computes the current genrated by the ionic Channel
        
        - Inputs:
            + **t**: instant in ms.
            + **V_mV**: membrane potential of the compartment in mV.
        
        - Outputs:
            + Ionic current, in nA
        '''        
         
        for i in xrange(0, self.lenStates): 
            self.condState[i].computeStateValue(t)        
                          
        return self.compCond(V_mV, self.gmax_muS, self.condState, self.EqPot_mV)

    cdef reset(self):
        for i in xrange(self.lenStates):
            self.condState[i].reset()
