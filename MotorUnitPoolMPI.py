'''
Created on Jul, 28 2015

@author: root
'''

import numpy as np
from MotorUnit import MotorUnit
from MuscularActivation1 import MuscularActivation
from MuscleNoHill import MuscleNoHill
from MuscleHill import MuscleHill
from scipy.sparse import lil_matrix

 

from mpi4py import MPI
import copy
import itertools

class MotorUnitPool(object):
    '''
    Class that implements a motor unit pool. Encompasses a set of motor
    units that controls a single  muscle.
    '''


    def __init__(self, conf, pool):
        '''
        Constructor

        - Inputs:
            + **conf**: Configuration object with the simulation parameters.

            + **pool**: string with Motor unit pool to which the motor unit belongs.
        '''

        ## Indicates that is Motor Unit pool.
        self.kind = 'MU'

        ## Configuration object with the simulation parameters.
        self.conf = conf
        ## String with Motor unit pool to which the motor unit belongs.
        self.pool = pool
        MUnumber_S = int(conf.parameterSet('MUnumber_' + pool + '-S', pool, 0))
        MUnumber_FR = int(conf.parameterSet('MUnumber_' + pool + '-FR', pool, 0))
        MUnumber_FF = int(conf.parameterSet('MUnumber_' + pool + '-FF', pool, 0))
        ## Number of motor units.
        self.MUnumber = MUnumber_S + MUnumber_FR + MUnumber_FF
        
        ## List of MotorUnit objects.
        self.unit = []
        
        
        for i in xrange(0, self.MUnumber): 
            if i < MUnumber_S:
                self.unit.append(MotorUnit(conf, pool, i, 'S'))
            elif i < MUnumber_S + MUnumber_FR:
                self.unit.append(MotorUnit(conf, pool, i, 'FR'))
            else:
                self.unit.append(MotorUnit(conf, pool, i, 'FF'))

        ## Vector with the instants of spikes in the soma compartment, in ms.            
        self.poolSomaSpikes = np.array([])    
        ## Vector with the instants of spikes in the terminal, in ms.
        self.poolTerminalSpikes = np.array([])
        
        #activation signal
        self.Activation = MuscularActivation(self.conf,self.pool, self.MUnumber,self.unit)
        
        #Force
        ## String indicating whther a Hill model is used or not. For now, it can be *No*.
        self.hillModel = conf.parameterSet('hillModel',pool, 0)
        if self.hillModel == 'No': 
            self.Muscle = MuscleNoHill(self.conf, self.pool, self.MUnumber, MUnumber_S, self.unit)
        else:
            self.Muscle = MuscleHill(self.conf, self.pool, self.MUnumber, MUnumber_S, self.unit)
        
        # MPI
        self.comm = MPI.COMM_WORLD
        self.size = self.comm.Get_size ()
        self.rank = self.comm.Get_rank ()
        self.procSize = self.MUnumber / (self.size - 1)
        self.unitSlice=[]
        self.aux=[]
        self.unitSlice = self.unit[(self.rank - 1) * self.procSize:self.rank * self.procSize]

        print 'size ' + str(self.size)
        print 'rank ' + str(self.rank)
        print 'procSize ' + str(self.procSize)
        
    def atualizeMotorUnitPool(self, t):
        '''
        Update all parts of the Motor Unit pool. It consists
        to update all motor units, the activation signal and
        the muscle force.

        - Inputs:
            + **t**: current instant, in ms.
        '''
        muscleSpike = []
        # Rank 0 process is responsible for force computations
        if self.rank != 0:
            for i in self.unitSlice: 
                i.atualizeMotorUnit(t)
                # Some terminalSpikeTrain values need to be empty lists
                if not i.terminalSpikeTrain:
                    muscleSpike.append([])
                else:
                    muscleSpike.append(i.terminalSpikeTrain[-1][0])
            # Each process send spikes to rank 0 process
            self.comm.send(muscleSpike, dest=0, tag=self.rank)

        if self.rank == 0:
            unit = []
            for ps in xrange(self.size - 1):
                unit.extend (self.comm.recv(source=ps + 1, tag=ps + 1))

            self.Activation.atualizeActivationSignal(t, unit)
            self.Muscle.atualizeForce(self.Activation.activation_Sat)

    def listSpikes(self):
        '''
        List the spikes that occurred in the soma and in
        the terminal of the different motor units.
        '''
        
        #if self.rank != 0:
        #    self.comm.send(, dest=0, tag=self.rank)

        # TODO this for each process, with appropriate length
        # Join results in mai
        for i in xrange(0,self.MUnumber):
            if i == 0:
                somaSpikeTrain = np.array(self.unit[i].somaSpikeTrain)
                terminalSpikeTrain = np.array(self.unit[i].terminalSpikeTrain)
            else:
                somaSpikeTrain = np.append(somaSpikeTrain, np.array(self.unit[i].somaSpikeTrain))
                terminalSpikeTrain = np.append(terminalSpikeTrain, np.array(self.unit[i].terminalSpikeTrain))
        self.poolSomaSpikes = somaSpikeTrain
        self.poolTerminalSpikes = terminalSpikeTrain
            
        self.poolSomaSpikes = np.reshape(self.poolSomaSpikes, (-1, 2))
        self.poolTerminalSpikes = np.reshape(self.poolTerminalSpikes, (-1, 2))
