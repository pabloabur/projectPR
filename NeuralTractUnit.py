'''
Created on Oct 26, 2015

@author: root
'''

from PointProcessGenerator import PointProcessGenerator
import math
import numpy as np
from multiprocessing import Pool
import functools



class NeuralTractUnit(object):
    '''
    Class that implements a neural tract unit. 
    It consists of a point process generator.
    '''

    
    def __init__(self, conf, pool, index):
        '''
        Constructor

        - Inputs:
            + **conf**: Configuration object with the simulation parameters.

            + **pool**: string with the name of the Neural tract.

            + **index**: integer corresponding to the neural tract unit identification.

        '''     
        # point process generator data
        ## Integer order of the Gamma distribution.     
        self.GammaOrder = int(conf.parameterSet('GammaOrder_' + pool, pool, index))
        
        ## A PointProcessGenerator object, corresponding the generator of
        ## spikes of the neural tract unit.   
        self.spikesGenerator = PointProcessGenerator(self.GammaOrder, index)  
        ## List of the spikes of the neural tract unit.       
        self.terminalSpikeTrain = self.spikesGenerator.points
        
        
        
         
        # Build synapses       
        ## 
        self.SynapsesOut = []
        self.transmitSpikesThroughSynapses = []
        self.indicesOfSynapsesOnTarget = []
        
        ## Integer corresponding to the neural tract unit identification.
        self.index = index
    
      
    def atualizeNeuralTractUnit(self, t, FR):
        '''

        - Inputs:
            + **t**:

            + **FR**:
        '''        
        
        self.spikesGenerator.atualizeGenerator(t, FR)
        if self.terminalSpikeTrain and abs(t - self.terminalSpikeTrain[-1][0]) < 1e-3: self.transmitSpikes(t)
        
        
    
    def transmitSpikes(self, t):
        '''

        - Inputs:
            + **t**:
        '''
        for i in xrange(len(self.indicesOfSynapsesOnTarget)):
            self.transmitSpikesThroughSynapses[i].receiveSpike(t, self.indicesOfSynapsesOnTarget[i])
        