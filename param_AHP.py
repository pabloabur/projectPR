import matplotlib.pyplot as plt
import numpy as np
from inspyred import ec
from random import Random

from Configuration import Configuration
from InterneuronPoolOpt import InterneuronPool
from SynapsesFactory import SynapsesFactory

def generate_netparams(random, args):
    size = args.get('num_inputs')
    initialParams = [random.uniform(minParamValues[i], maxParamValues[i]) for i in range(size)]
    return initialParams

#gmax_Kf:RC_ext-@soma,1000,
#gmax_Ks:RC_ext-@soma,4000,
#def evaluate_netparams(candidates, args):
#    fitnessCandidates = []
#
#    for icand,cand in enumerate(candidates):
#        # modify network params based on this candidate params (genes)
#        tut2.netParams.connParams['S->M']['probability'] = cand[0]
#        tut2.netParams.connParams['S->M']['weight'] = cand[1]
#        tut2.netParams.connParams['S->M']['delay'] = cand[2]
#
#        # create network
#        sim.createSimulate(netParams=tut2.netParams, simConfig=tut2.simConfig)
#
#        # calculate firing rate
#        numSpikes = float(len(sim.simData['spkt']))
#        numCells = float(len(sim.net.cells))
#        duration = tut2.simConfig.duration/1000.0
#        netFiring = numSpikes/numCells/duration
#
#        # calculate fitness for this candidate
#        fitness = abs(targetFiring - netFiring)  # minimize absolute difference in firing rate
#
#        # add to list of fitness for each candidate
#        fitnessCandidates.append(fitness)
#
#        # print candidate parameters, firing rate, and fitness
#        print '\n CHILD/CANDIDATE %d: Network with prob:%.2f, weight:%.2f, delay:%.1f \n  firing rate: %.1f, FITNESS = %.2f \n'\
#        %(icand, cand[0], cand[1], cand[2], netFiring, fitness)
#
#
#    return fitnessCandidates

def simulator(numberRC, duration, current, newParametrization):

    conf = Configuration('confuchiyama.rmto')

    # Number of cells
    idx = np.where(conf.confArray['f0']=='Number_RC_ext')[0][0]
    conf.confArray['f1'][idx] = numberRC
    
    # Duration of simulation
    conf.simDuration_ms = duration
    
    # Parameters from java
    ## Threshold (makes a rheobase of 0.5 nA)
    idx = np.where(conf.confArray['f0']=='threshold:RC_ext-')[0][0]
    conf.confArray['f1'][idx] = 0.1724138
    conf.confArray['f2'][idx] = 0.1724138
        
    ## Connectivity
    idx = np.where(conf.confArray['f0']=='Con:RC_ext->SOL-S@soma|inhibitory')[0][0]
    conf.confArray['f1'][idx] = 100
    idx = np.where(conf.confArray['f0']=='Con:RC_ext->SOL-FR@soma|inhibitory')[0][0]
    conf.confArray['f1'][idx] = 100
    idx = np.where(conf.confArray['f0']=='Con:RC_ext->SOL-FF@soma|inhibitory')[0][0]
    conf.confArray['f1'][idx] = 100
    idx = np.where(conf.confArray['f0']=='Con:SOL-S>RC_ext-@soma|excitatory')[0][0]
    conf.confArray['f1'][idx] = 100
    idx = np.where(conf.confArray['f0']=='Con:SOL-FR>RC_ext-@soma|excitatory')[0][0]
    conf.confArray['f1'][idx] = 100
    idx = np.where(conf.confArray['f0']=='Con:SOL-FF>RC_ext-@soma|excitatory')[0][0]
    conf.confArray['f1'][idx] = 100

    ## Conductances
    idx = np.where(conf.confArray['f0']=='gmax:RC_ext->SOL-S@soma|inhibitory')[0][0]
    conf.confArray['f1'][idx] = 0.44
    idx = np.where(conf.confArray['f0']=='gmax:RC_ext->SOL-FR@soma|inhibitory')[0][0]
    conf.confArray['f1'][idx] = 0.3
    idx = np.where(conf.confArray['f0']=='gmax:RC_ext->SOL-FF@soma|inhibitory')[0][0]
    conf.confArray['f1'][idx] = 0.24
    idx = np.where(conf.confArray['f0']=='gmax:SOL-S>RC_ext-@soma|excitatory')[0][0]
    conf.confArray['f1'][idx] = 0.15
    idx = np.where(conf.confArray['f0']=='gmax:SOL-FR>RC_ext-@soma|excitatory')[0][0]
    conf.confArray['f1'][idx] = 0.17
    idx = np.where(conf.confArray['f0']=='gmax:SOL-FF>RC_ext-@soma|excitatory')[0][0]
    conf.confArray['f1'][idx] = 0.3

    ## Morphology
    idx = np.where(conf.confArray['f0']=='d@soma:RC_ext-')[0][0]
    conf.confArray['f1'][idx] = 64.77885
    conf.confArray['f2'][idx] = 64.77885
    idx = np.where(conf.confArray['f0']=='l@soma:RC_ext-')[0][0]
    conf.confArray['f1'][idx] = 285
    conf.confArray['f2'][idx] = 285
    idx = np.where(conf.confArray['f0']=='res@soma:RC_ext-')[0][0]
    conf.confArray['f1'][idx] = 8000
    conf.confArray['f2'][idx] = 8000
    
    pools = dict()
    pools[0] = InterneuronPool(conf, 'RC', 'ext')

    Syn = SynapsesFactory(conf, pools)

    t = np.arange(0.0, conf.simDuration_ms, conf.timeStep_ms)

    RC_mV = np.zeros_like(t)
    for i in xrange(0, len(t)):
        if t[i]>10 and t[i]<20:
            for j in xrange(len(pools[0].unit)):
                pools[0].iInjected[j] = current
        else:
            for j in xrange(len(pools[0].unit)):
                pools[0].iInjected[j] = 0
        #pools[1].atualizePool(t[i]) # RC synaptic Noise
        pools[0].atualizeInterneuronPool(t[i]) # RC pool
        RC_mV[i] = pools[0].unit[0].v_mV[0]

    pools[0].listSpikes()

    return RC_mV, t

nRC = 1
t = 200
newParams = False
i = 0.5
#i = 0.1 # for time constant
RCMembrane, t = simulator(nRC, t, i, newParams)

# Calculate AHP duration
i = 0
f = 0
for j in xrange(len(RCMembrane)):
    if RCMembrane[j] < 0:
        i = j
        break
for j in xrange(i, len(RCMembrane)):
    if np.isclose(RCMembrane[j], 0.0, atol=0.1):
        f = j
        break
AHPa = min(RCMembrane)
AHPt = t[f]-t[i]
print "AHP amplitude is " + str(AHPa)
print "AMP duration is " + str(AHPt)

plt.figure()
plt.plot(t, RCMembrane)
plt.title('RC Membrane Potential')
plt.show()
