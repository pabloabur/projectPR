import math
import numpy as np

cdef class PulseConductanceState:
    '''
    Implements the Destexhe pulse approximation of the solution of 
    the states of the Hodgkin-Huxley neuron model.
    '''
    # TODO careful manipulation of strings
    cdef str kind, actType, neuronKind, compKind, pool
    #cdef public double Beta, Alpha, PulseDuration, TimeStep, AlphaExp, BetaExp, endOfPulse_ms
    cdef double alpha_ms1, beta_ms1, PulseDur_ms, AlphaExp, BetaExp, endOfPulse_ms
    cdef double value
    cdef bint state
    cdef unsigned int index
   

    def __init__(self, str kind, conf, str pool, str neuronKind, str compKind, unsigned int index):
    #def __init__(self, kind, Beta, Alpha, PulseDuration, TimeStep):
        self.kind = kind
        self.beta_ms1 = <double>conf.parameterSet('beta_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool, index)
        self.alpha_ms1 = <double>conf.parameterSet('alpha_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool,index)
        self.PulseDur_ms = <double>conf.parameterSet('PulseDur_' + kind, pool, index)
        #self.beta_ms1 = Beta
        #self.alpha_ms1 = Alpha
        #self.PulseDur_ms = PulseDuration
        self.value = <double>0
        #self.value = 0.0
        self.state = False

        self.AlphaExp = math.exp(-self.alpha_ms1 * conf.timeStep_ms)
        self.BetaExp = math.exp(-self.beta_ms1 * conf.timeStep_ms)
        #self.AlphaExp = math.exp(-self.alpha_ms1 * TimeStep)
        #self.BetaExp = math.exp(-self.beta_ms1 * TimeStep)

        self.endOfPulse_ms = self.PulseDur_ms

        if (self.kind == 'm'):
            self.actType = 'activation'
        if (self.kind == 'h'):
            self.actType = 'inactivation'
        if (self.kind == 'n'):
            self.actType = 'activation'
        if (self.kind == 'q'):
            self.actType = 'activation'
        if (self.kind == 'mp'):
            self.actType = 'activation'
        if (self.kind == 's'):
            self.actType = 'activation'
        if (self.kind == 'qh'):
            self.actType = 'inactivation'

        if (self.actType == 'activation'):
            self.computeStateValue = self.computeStateValueActivation            
        else:
            self.computeStateValue = self.computeStateValueInactivation

    cdef changeState(self, double t):
        '''
        Void function that modify the current situation (true/false)
        of the state.

        - Inputs:
            + **t**: current instant, in ms.
        '''

        self.state = not self.state
        self.endOfPulse_ms = self.PulseDur_ms + t

    cdef computeStateValueActivation(self, double t):
        '''
        Compute the state value by using the approximation of Destexhe (1997) to
        compute the Hodgkin-Huxley states of *activation* type.

        - Input:
            + **t**: current instant, in ms.

        The value of the state \f$v\f$ is computed according to the following
        equation before and after the pulse:

        \f{equation}{
            v(t) = v_0\exp[-\beta(t-t_0)]
        \f} 

        and according to the following equation during the pulse:

        \f{equation}{
            v(t) = 1 + (v_0 - 1)\exp[-\alpha(t-t_0)]
        \f} 
        where \f$t_0\f$ is the time at which the pulse changed
        the value (on to off or off to on) and \f$v_0\f$ is value
        of the state at that time.
        '''

        if not self.state:
            self.value *= self.BetaExp
        else:
            if t > self.endOfPulse_ms:
                self.changeState(t)
                self.value *= self.BetaExp                 
            else: 
                self.value = (self.value - 1) * self.AlphaExp + 1                
        
    cdef computeStateValueInactivation(self, double t):
        '''
        Compute the state value by using the approximation of Destexhe (1997) to
        compute the Hodgkin-Huxley states of *inactivation* type.

        - Input:
            + **t**: current instant, in ms.

        The value of the state \f$v\f$ is computed according to the following
        equation before and after the pulse:

        \f{equation}{
            v(t) = v_0\exp[-\beta(t-t_0)]
        \f} 

        and according to the following equation during the pulse:

        \f{equation}{
            v(t) = 1 + (v_0 - 1)\exp[-\alpha(t-t_0)]
        \f} 
        where \f$t_0\f$ is the time at which the pulse changed
        the value (on to off or off to on) and \f$v_0\f$ is value
        of the state at that time.
        '''

        if not self.state:
            self.value = (self.value - 1) * self.AlphaExp + 1
        else:
            if t > self.endOfPulse_ms:
                self.changeState(t)
                self.value = (self.value - 1) * self.AlphaExp + 1
            else: self.value *= self.BetaExp  

    cdef reset(self):
        self.value = <double>0
        self.endOfPulse_ms = self.PulseDur_ms

        
