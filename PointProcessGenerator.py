'''
    Neuromuscular simulator in Python.
    Copyright (C) 2016  Renato Naville Watanabe

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
import numpy.random
import math
from numba import jit






@jit
def  gammaPoint(GammaOrder, GammaOrderInv):
    '''
    Generates a number according to a Gamma Distribution with an integer order **GammaOrder**.

    - Inputs:
        + **GammaOrder**: integer order of the Gamma distribution.

        + **GammaOrderInv**:  inverse of the GammaOrder. This is necessary
                                for computational efficiency.

    - Outputs:
        + The number generated from the Gamma distribution.

    The number is generated according to:

    \f{equation}{
        \Gamma = -\frac{1}{\lambda}\ln(\limits\prod_{i=1}^{\lambda} U(0,1))
    \f}
    where \f$\lambda\f$ is the order of the Gamma distribution and U(a,b) is
    a uniform distribution from a to b.

    '''
    aux = np.prod(numpy.random.uniform(0.0, 1.0, size=GammaOrder))

    return - GammaOrderInv * np.log(aux)

class PointProcessGenerator(object):
    '''
    Generator of point processes.
    '''
    
    def __init__(self, GammaOrder, index):
        '''
        Constructor

        - Inputs:
            + **GammaOrder**: integer order of the Gamma distribution.

            + **index**: integer corresponding to the unit order in the pool.
        '''
        ## Integer order of the Gamma distribution. 
        ## Gamma order 1 is Poisson process and order 10 is a Gaussian process.
        self.GammaOrder = int(GammaOrder);
        ## Inverse of the GammaOrder. This is necessary
        ## for computational efficiency.
        self.GammaOrderInv = 1.0 / GammaOrder
        ## Integer corresponding to the unit order in the pool to which this
        ## generator is associated.
        self.index = index

        ## Auxiliary variable cummulating a value that indicates
        ## whether there will be a new spike or not.
        self.threshold = gammaPoint(self.GammaOrder, self.GammaOrderInv)

        
        ## List of spike instants of the generator.
        self.points = []

    def atualizeGenerator(self, t, firingRate):
        '''

        - Inputs:
            + **t**: current instant, in ms.

            + **firingRate**: instant firing rate, in spikes/s.
        '''
        self.threshold -= firingRate
        if self.threshold <= 0 and t != 0:
            self.points.append([t, self.index])
            self.threshold = gammaPoint(self.GammaOrder, self.GammaOrderInv)
            
