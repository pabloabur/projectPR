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


import numpy as np

from NeuralTract import NeuralTract
from SynapticNoise import SynapticNoise

class SynapsesFactory(object):
    '''
    Class to build all the synapses in the system.
    '''


    def __init__(self, conf, pools):
        '''
        Constructor

        - Inputs:
            + **conf**: Configuration object with the simulation parameters.

            + **pools**: list of all the pools in the system.

        
        '''
        ## Total number of synapses in the system.
        self.numberOfSynapses = 0
        #pools.append(NeuralTract(conf, 'NoiseRC'))

        # The 'Out' suffix is used here to identify the presynaptic neurons.
        for poolOut in xrange(len(pools)):
            for unitOut in xrange(len(pools[poolOut].unit)):
                # Get all the synapses that unitOut from poolOut can make.
                # They may or may not be connected later, depending on simulation settings
                pools[poolOut].unit[unitOut].SynapsesOut = conf.determineSynapses(pools[poolOut].pool + '-' + 
                                                                                  pools[poolOut].unit[unitOut].kind)
                # The 'In' suffix, on the other hand, is for the postsynaptic neurons.
                # The loop will prepare all possible sinapses from pools[poolOut].unit[unitOut]
                for synapseIn in xrange(len(pools[poolOut].unit[unitOut].SynapsesOut)):
                    presynapticPool = pools[poolOut].unit[unitOut].SynapsesOut[synapseIn][0]
                    presynapticType = pools[poolOut].unit[unitOut].SynapsesOut[synapseIn][1]
                    presynapticCompartment = pools[poolOut].unit[unitOut].SynapsesOut[synapseIn][2]
                    presynapticKind = pools[poolOut].unit[unitOut].SynapsesOut[synapseIn][3]
                    conn = float(conf.parameterSet('Con:' + pools[poolOut].pool + '-' 
                                                   + pools[poolOut].unit[unitOut].kind + '>'
                                                   + presynapticPool
                                                   + '-'
                                                   + presynapticType
                                                   + '@'
                                                   + presynapticCompartment
                                                   + '|'
                                                   + presynapticKind,
                                                   '', 0)) / 100.0
                    gmax = float(conf.parameterSet('gmax:' + pools[poolOut].pool + '-'
                                                   + pools[poolOut].unit[unitOut].kind + '>'
                                                   + presynapticPool
                                                   + '-'
                                                   + presynapticType
                                                   + '@'
                                                   + presynapticCompartment
                                                   + '|'
                                                   + presynapticKind,
                                                   '', 0))
                    delay = float(conf.parameterSet('delay:' + pools[poolOut].pool + '-'
                                                    + pools[poolOut].unit[unitOut].kind + '>'
                                                    + presynapticPool
                                                    + '-'
                                                    + presynapticType
                                                    + '@'
                                                    + presynapticCompartment
                                                    + '|'
                                                    + presynapticKind,
                                                    '', 0))
                    declineFactor = float(conf.parameterSet('dec:' + pools[poolOut].pool + '-'
                                                    + pools[poolOut].unit[unitOut].kind + '>'
                                                    + presynapticPool
                                                    + '-'
                                                    + presynapticType
                                                    + '@'
                                                    + presynapticCompartment
                                                    + '|'
                                                    + presynapticKind,
                                                    '', 0))
                    dyn = conf.parameterSet('dyn:' + pools[poolOut].pool + '-'
                                            + pools[poolOut].unit[unitOut].kind + '>'
                                            + presynapticPool
                                            + '-' + presynapticType
                                            + '@' + presynapticCompartment
                                            + '|' + presynapticKind,
                                            '', 0)
                    if dyn != 'None':
                        var = float(conf.parameterSet('var:' + pools[poolOut].pool + '-'
                                                      + pools[poolOut].unit[unitOut].kind + '>'
                                                      + presynapticPool
                                                      + '-' + presynapticType
                                                      + '@' + presynapticCompartment
                                                      + '|' + presynapticKind,
                                                      '', 0))
                        tau = float(conf.parameterSet('tau:' + pools[poolOut].pool + '-'
                                                      + pools[poolOut].unit[unitOut].kind + '>'
                                                      + presynapticPool
                                                      + '-' + presynapticType
                                                      + '@' + presynapticCompartment
                                                      + '|' + presynapticKind,
                                                      '', 0))
                    else:
                        var = 0
                        tau = 100000
                    # This loop will determine which neurons will actually receive the synapse just built.
                    for poolIn in xrange(len(pools)):
                        if pools[poolIn].pool in presynapticPool:
                            for unitIn in xrange(len(pools[poolIn].unit)):
                                for compartmentIn in xrange(len(pools[poolIn].unit[unitIn].compartment)):
                                    if presynapticPool == pools[poolIn].pool and presynapticType == pools[poolIn].unit[unitIn].kind and presynapticCompartment == pools[poolIn].unit[unitIn].compartment[compartmentIn].kind:
                                        if np.isfinite(declineFactor):
                                            neuronsDistance = np.abs(pools[poolIn].unit[unitIn].position_mm
                                                    - pools[poolOut].unit[unitOut].position_mm)
                                            weight = 1
                                            Pconn = conn*declineFactor / (declineFactor + neuronsDistance**2)
                                        else:
                                            weight = 1
                                            Pconn = conn
                                        if np.random.uniform(0.0, 1.0) <= Pconn:
                                            for synapse in xrange(len(pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn)): 
                                                if pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn[synapse].kind == presynapticKind:
                                                    pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn[synapse].addConductance(gmax*weight, delay, dyn, var, tau)
                                                    pools[poolOut].unit[unitOut].transmitSpikesThroughSynapses.append(pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn[synapse])                                                            
                                                    pools[poolOut].unit[unitOut].indicesOfSynapsesOnTarget.append(len(pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn[synapse].gmax_muS) - 1)
                                                    self.numberOfSynapses += 1

        print 'All the ' + str(self.numberOfSynapses) +  ' synapses were built'


        ## Total number of synaptic noises in the system.
        self.numberOfSynapticNoise = 0

        NoiseSynapsesOut = conf.determineSynapses('Noise')
        for synapseIn in xrange(len(NoiseSynapsesOut)):
            pools[len(pools)] = SynapticNoise(conf, NoiseSynapsesOut[synapseIn][0])
            poolOut = len(pools) - 1
            gmax = float(conf.parameterSet('gmax:Noise>' + NoiseSynapsesOut[synapseIn][0] 
                                           + '-' + NoiseSynapsesOut[synapseIn][1]
                                           + '@' + NoiseSynapsesOut[synapseIn][2] + '|'
                                           + NoiseSynapsesOut[synapseIn][3],
                                           '', 0))
            delay = float(conf.parameterSet('delay:Noise>' + NoiseSynapsesOut[synapseIn][0]
                                            + '-' + NoiseSynapsesOut[synapseIn][1]
                                            + '@' + NoiseSynapsesOut[synapseIn][2] + '|'
                                            + NoiseSynapsesOut[synapseIn][3],
                                            '', 0))
            declineFactor = float(conf.parameterSet('dec:Noise>' + NoiseSynapsesOut[synapseIn][0]
                                                    + '-' + NoiseSynapsesOut[synapseIn][1]
                                                    + '@' + NoiseSynapsesOut[synapseIn][2] + '|'
                                                    + NoiseSynapsesOut[synapseIn][3],
                                                    '', 0))
            dyn = conf.parameterSet('dyn:Noise>' + NoiseSynapsesOut[synapseIn][0] 
                                    + '-' + NoiseSynapsesOut[synapseIn][1]
                                    + '@' + NoiseSynapsesOut[synapseIn][2] + '|'
                                    + NoiseSynapsesOut[synapseIn][3],
                                    '', 0)
            if dyn != 'None':
                var = float(conf.parameterSet('var:Noise>' + NoiseSynapsesOut[synapseIn][0]
                                              + '-' + NoiseSynapsesOut[synapseIn][1]
                                              + '@' + NoiseSynapsesOut[synapseIn][2] + '|' 
                                              + NoiseSynapsesOut[synapseIn][3],
                                              '', 0))
                tau = float(conf.parameterSet('tau:Noise>' + NoiseSynapsesOut[synapseIn][0]
                                              + '-' + NoiseSynapsesOut[synapseIn][1]
                                              + '@' + NoiseSynapsesOut[synapseIn][2]
                                              + '|' + NoiseSynapsesOut[synapseIn][3],
                                              '', 0))
            else:
                var = 0
                tau = 10000
            for unitOut in xrange(len(pools[poolOut].unit)):
                for poolIn in xrange(len(pools)):
                    if NoiseSynapsesOut[synapseIn][0] == pools[poolIn].pool and pools[poolIn].kind != 'SN':
                        for unitIn in xrange(len(pools[poolIn].unit)):
                            for compartmentIn in xrange(len(pools[poolIn].unit[unitIn].compartment)):
                                if NoiseSynapsesOut[synapseIn][1] == pools[poolIn].unit[unitIn].kind and NoiseSynapsesOut[synapseIn][2] == pools[poolIn].unit[unitIn].compartment[compartmentIn].kind and pools[poolIn].unit[unitIn].index == pools[poolOut].unit[unitOut].index:
                                    for synapse in xrange(len(pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn)): 
                                        if pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn[synapse].kind == NoiseSynapsesOut[synapseIn][3]:
                                            if np.isfinite(declineFactor):
                                                neuronsDistance = np.abs(pools[poolIn].unit[unitIn].position_mm
                                                                            - pools[poolOut].unit[unitOut].position_mm)
                                                weight = declineFactor / (declineFactor + neuronsDistance**2)
                                            else:
                                                weight = 1
                                            pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn[synapse].addConductance(gmax*weight, delay, dyn, var, tau)
                                            pools[poolOut].unit[unitOut].transmitSpikesThroughSynapses.append(pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn[synapse])                                                            
                                            pools[poolOut].unit[unitOut].indicesOfSynapsesOnTarget.append(len(pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn[synapse].gmax_muS) - 1)
                                            self.numberOfSynapticNoise += 1
                  

        

        print 'All the ' + str(self.numberOfSynapticNoise) +  ' synaptic noises were built'                   
        
