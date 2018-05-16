'''
    Neuromuscular simulator in Python.
    Copyright (C) 2017  Renato Naville Watanabe

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: renato.watanabe@usp.br
'''


import math
#from numba import jit
import numpy as np


#@jit
def compValOn(v0, alpha, beta, t, t0):
    '''
    Time course of the state during the pulse for the *inactivation* states
    and before and after the pulse for the *activation* states.

    The value of the state \f$v\f$ is computed according to the following
    equation:

    \f{equation}{
        v(t) = v_0\exp[-\beta(t-t_0)]
    \f} 
    where \f$t_0\f$ is the time at which the pulse changed
    the value (on to off or off to on) and \f$v_0\f$ is value
    of the state at that time.
    '''
    return v0 * np.exp(beta * (t0 - t))

#@jit
def compValOff(v0, alpha, beta, t, t0):
    '''
    Time course of the state during the pulse for the *activation* states
    and before and after the pulse for the *inactivation* states.

    The value of the state \f$v\f$ is computed according to the following
    equation:

    \f{equation}{
        v(t) = 1 + (v_0 - 1)\exp[-\alpha(t-t_0)]
    \f} 
    where \f$t_0\f$ is the time at which the pulse changed
    the value (on to off or off to on) and \f$v_0\f$ is value
    of the state at that time.
    '''
    return 1.0 + (v0 - 1.0)  *  np.exp(alpha * (t0 - t))

cdef class PulseConductanceState:
    '''
    Implements the Destexhe pulse approximation of the solution of 
    the states of the Hodgkin-Huxley neuron model.
    '''
    # TODO careful manipulation of strings
    cdef public str kind, actType
    cdef public double endOfPulse_ms
    cdef public double value
    cdef public bint state
    cdef public double alpha_ms1, beta_ms1, PulseDur_ms, AlphaExp, BetaExp
    '''
    Implements the Destexhe pulse approximation of the solution of 
    the states of the Hodgkin-Huxley neuron model.
    '''

   
    def __init__(self, kind, conf, pool, neuronKind, compKind, index):
        '''
        Initializes the pulse conductance state.

        Variables:
            + **kind**: string with type of the state (m, h, n, q).
            
            + **conf**:  an instance of the Configuration class with the functions to correctly parameterize the model. See the Configuration class.
            
            + **pool**: string with the pool that this state belongs.
            
            + **neuronKind**: string with the type of the motor unit. It used for 
            motoneurons. It can be *S* (slow), *FR* (fast and resistant), and *FF* 
            (fast and fatigable). 

            + **compKind**: The kind of compartment that the Channel belongs. 
            For now, it can be *soma*, *dendrite*, *node* or *internode*.

            + **index**: the index of the unit that this state belongs.                    
        '''
        self.kind = kind
        self.value = <double>0
        
        
        
       
        
        self.state = False
        
        self.beta_ms1 = <double>conf.parameterSet('beta_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool, index)
        self.alpha_ms1 = <double>conf.parameterSet('alpha_' + kind + ':' + pool + '-' + neuronKind + '@' + compKind, pool,index)
        self.PulseDur_ms = <double>conf.parameterSet('PulseDur_' + kind, pool, index)

        self.AlphaExp = math.exp(-self.alpha_ms1 * conf.timeStep_ms)
        self.BetaExp = math.exp(-self.beta_ms1 * conf.timeStep_ms)

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

    #@profile
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
        
    
    #@profile
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

    def reset(self):
        '''

        '''
        self.value = float(0)
        self.endOfPulse_ms = self.PulseDur_ms

        
#import math
#import numpy as np
#
#def compValOn(v0, alpha, beta, t, t0):
#    '''
#    Time course of the state during the pulse for the *inactivation* states
#    and before and after the pulse for the *activation* states.
#
#    The value of the state \f$v\f$ is computed according to the following
#    equation:
#
#    \f{equation}{
#        v(t) = v_0\exp[-\beta(t-t_0)]
#    \f} 
#    where \f$t_0\f$ is the time at which the pulse changed
#    the value (on to off or off to on) and \f$v_0\f$ is value
#    of the state at that time.
#    '''
#    return v0 * np.exp(beta * (t0 - t))
#
#def compValOff(v0, alpha, beta, t, t0):
#    '''
#    Time course of the state during the pulse for the *activation* states
#    and before and after the pulse for the *inactivation* states.
#
#    The value of the state \f$v\f$ is computed according to the following
#    equation:
#
#    \f{equation}{
#        v(t) = 1 + (v_0 - 1)\exp[-\alpha(t-t_0)]
#    \f} 
#    where \f$t_0\f$ is the time at which the pulse changed
#    the value (on to off or off to on) and \f$v_0\f$ is value
#    of the state at that time.
#    '''
#    return 1.0 + (v0 - 1.0)  *  np.exp(alpha * (t0 - t))
#
#cdef class PulseConductanceState:
#    '''
#    Implements the Destexhe pulse approximation of the solution of 
#    the states of the Hodgkin-Huxley neuron model.
#    '''
#    # TODO careful manipulation of strings
#    cdef public str kind, actType
#    cdef public double Beta, Alpha, PulseDuration, TimeStep, AlphaExp, BetaExp, endOfPulse_ms
#    cdef public double value
#    cdef public bint state
#   
#    # TODO import Configuration instead
#    def __init__(self, kind, Beta, Alpha, PulseDuration, TimeStep):
#        '''
#        Initializes the pulse conductance state.
#
#        Variables:
#            + **kind**: string with type of the state (m, h, n, q).
#            
#        '''
#        self.kind = kind
#        self.beta_ms1 = Beta
#        self.alpha_ms1 = Alpha
#        self.PulseDur_ms = PulseDuration
#        self.value = 0.0
#        self.state = False
#        self.state = not self.state
#
#        self.AlphaExp = math.exp(-self.alpha_ms1 * TimeStep)
#        self.BetaExp = math.exp(-self.beta_ms1 * TimeStep)
#
#        self.endOfPulse_ms = self.PulseDur_ms
#
#        if (self.kind == 'm'):
#            self.actType = 'activation'
#        if (self.kind == 'h'):
#            self.actType = 'inactivation'
#        if (self.kind == 'n'):
#            self.actType = 'activation'
#        if (self.kind == 'q'):
#            self.actType = 'activation'
#        if (self.kind == 'mp'):
#            self.actType = 'activation'
#        if (self.kind == 's'):
#            self.actType = 'activation'
#        if (self.kind == 'qh'):
#            self.actType = 'inactivation'
#
#        if (self.actType == 'activation'):
#            self.computeStateValue = self.computeStateValueActivation            
#        else:
#            self.computeStateValue = self.computeStateValueInactivation
#
#
#    def changeState(self, t):
#        '''
#        Void function that modify the current situation (true/false)
#        of the state.
#
#        - Inputs:
#            + **t**: current instant, in ms.
#        '''
#
#        self.state = not self.state
#        self.endOfPulse_ms = self.PulseDur_ms + t
#
#    #@profile
#    def computeStateValueActivation(self, t):
#        '''
#        Compute the state value by using the approximation of Destexhe (1997) to
#        compute the Hodgkin-Huxley states of *activation* type.
#
#        - Input:
#            + **t**: current instant, in ms.
#
#        The value of the state \f$v\f$ is computed according to the following
#        equation before and after the pulse:
#
#        \f{equation}{
#            v(t) = v_0\exp[-\beta(t-t_0)]
#        \f} 
#
#        and according to the following equation during the pulse:
#
#        \f{equation}{
#            v(t) = 1 + (v_0 - 1)\exp[-\alpha(t-t_0)]
#        \f} 
#        where \f$t_0\f$ is the time at which the pulse changed
#        the value (on to off or off to on) and \f$v_0\f$ is value
#        of the state at that time.
#        '''
#
#        if not self.state:
#            self.value *= self.BetaExp
#        else:
#            if t > self.endOfPulse_ms:
#                self.changeState(t)
#                self.value *= self.BetaExp                 
#            else: 
#                self.value = (self.value - 1) * self.AlphaExp + 1                
#        
#    
#    #@profile
#    def computeStateValueInactivation(self, t):
#        '''
#        Compute the state value by using the approximation of Destexhe (1997) to
#        compute the Hodgkin-Huxley states of *inactivation* type.
#
#        - Input:
#            + **t**: current instant, in ms.
#
#        The value of the state \f$v\f$ is computed according to the following
#        equation before and after the pulse:
#
#        \f{equation}{
#            v(t) = v_0\exp[-\beta(t-t_0)]
#        \f} 
#
#        and according to the following equation during the pulse:
#
#        \f{equation}{
#            v(t) = 1 + (v_0 - 1)\exp[-\alpha(t-t_0)]
#        \f} 
#        where \f$t_0\f$ is the time at which the pulse changed
#        the value (on to off or off to on) and \f$v_0\f$ is value
#        of the state at that time.
#        '''
#
#        if not self.state:
#            self.value = (self.value - 1) * self.AlphaExp + 1
#        else:
#            if t > self.endOfPulse_ms:
#                self.changeState(t)
#                self.value = (self.value - 1) * self.AlphaExp + 1
#            else: self.value *= self.BetaExp  
#
#    def reset(self):
#        '''
#
#        '''
#        self.value = float(0)
#        self.endOfPulse_ms = self.PulseDur_ms
#
#        
