import matplotlib.pyplot as plt
import numpy as np
from inspyred import ec
from random import Random

from Configuration import Configuration
from InterneuronPoolOpt import InterneuronPool
from SynapsesFactory import SynapsesFactory

def my_constraint_function(RCMembrane):
    """Return the number of constraints that candidate violates."""
    # In this case, we'll just say that the point has to lie 
    # within a circle centered at (0, 0) of radius 1.
    if min(RCMembrane) == -30 or min(RCMembrane) == 0 or max(RCMembrane) == 120:
        return 1
    else:
        return 0

def my_evaluator(candidates, args):
    # The fitness will be how far the point is from
    # the origin. (We're maximizing, in this case.)
    # Note that the constraint heavily punishes individuals
    # who go beyond the unit circle. Therefore, these
    # two functions combined focus the evolution toward
    # finding individual who lie ON the circle.
    fitnessCandidates = []

    for icand,cand in enumerate(candidates):
        # Make simulations
        RCMembrane, t = simulator(cand[0], cand[1], cand[2], cand[3], cand[4], cand[5], cand[6],
                                  cand[7], cand[8], cand[9], cand[10])

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

        # minimize absolute differences
        f1 = abs(targetAHPa - AHPa)
        f2 = abs(targetAHPt - AHPt)

        # add to list of fitness for each candidate
        if my_constraint_function(RCMembrane) > 0:
            fitnessCandidates.append(ec.emo.Pareto([10000, 10000]))
        else:
            fitnessCandidates.append(ec.emo.Pareto([f1, f2]))

        print "\n Child Candidate %d: gKf = %d, gKs = %d, gNa = %d\n"%(icand, cand[0], cand[1], cand[2])
        print "AHP amplitude: %.2f; AHP duration: %.2f; \n"%(AHPa, AHPt)
        print fitnessCandidates

    return fitnessCandidates

def generate_netparams(random, args):
    size = args.get('num_inputs')
    initialParams = [random.uniform(minParamValues[i], maxParamValues[i]) for i in range(size)]
    return initialParams

def simulator(gKf, gKs, gNa, alpham, alphan, alphah, alphaq, betam, betan, betah, betaq):
    conf = Configuration('confuchiyama.rmto')

    # Number of cells
    idx = np.where(conf.confArray['f0']=='Number_RC_ext')[0][0]
    conf.confArray['f1'][idx] = 1
    
    # Duration of simulation
    conf.simDuration_ms = 100
    
    # Parameters from java
    ## conductances
    idx = np.where(conf.confArray['f0']=='gmax_Kf:RC_ext-@soma')[0][0]
    conf.confArray['f1'][idx] = gKf
    idx = np.where(conf.confArray['f0']=='gmax_Ks:RC_ext-@soma')[0][0]
    conf.confArray['f1'][idx] = gKs
    idx = np.where(conf.confArray['f0']=='gmax_Na:RC_ext-@soma')[0][0]
    conf.confArray['f1'][idx] = gNa
    idx = np.where(conf.confArray['f0']=='alpha_m:RC_ext-@soma')[0][0]
    conf.confArray['f1'][idx] = alpham
    idx = np.where(conf.confArray['f0']=='alpha_n:RC_ext-@soma')[0][0]
    conf.confArray['f1'][idx] = alphan
    idx = np.where(conf.confArray['f0']=='alpha_h:RC_ext-@soma')[0][0]
    conf.confArray['f1'][idx] = alphah
    idx = np.where(conf.confArray['f0']=='alpha_q:RC_ext-@soma')[0][0]
    conf.confArray['f1'][idx] = alphaq
    idx = np.where(conf.confArray['f0']=='beta_m:RC_ext-@soma')[0][0]
    conf.confArray['f1'][idx] = betam
    idx = np.where(conf.confArray['f0']=='beta_n:RC_ext-@soma')[0][0]
    conf.confArray['f1'][idx] = betan
    idx = np.where(conf.confArray['f0']=='beta_h:RC_ext-@soma')[0][0]
    conf.confArray['f1'][idx] = betah
    idx = np.where(conf.confArray['f0']=='beta_q:RC_ext-@soma')[0][0]
    conf.confArray['f1'][idx] = betaq

    ## Threshold (makes a rheobase of 0.5 nA)
    idx = np.where(conf.confArray['f0']=='threshold:RC_ext-')[0][0]
    conf.confArray['f1'][idx] = 0.1724138
    conf.confArray['f2'][idx] = 0.1724138
        
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
                pools[0].iInjected[j] = 0.5
        else:
            for j in xrange(len(pools[0].unit)):
                pools[0].iInjected[j] = 0
        #pools[1].atualizePool(t[i]) # RC synaptic Noise
        pools[0].atualizeInterneuronPool(t[i]) # RC pool
        RC_mV[i] = pools[0].unit[0].v_mV[0]

    pools[0].listSpikes()

    return RC_mV, t

if __name__ == '__main__':
    # create random seed for evolutionary computation algorithm --> my_ec = ec.EvolutionaryComputation(rand)
    rand = Random()
    #rand.seed(1)

    # targets
    targetAHPt = 30
    targetAHPa = -2

    # min and max allowed value for each AHP parameter optimized:
    #                 gKs,      gKf,    gNa
    minParamValues = [100,      100,    100,    0.1,        0.1,       0.1,       0.1,       0.1,       0.1,       0.1,       0.1]
    maxParamValues = [100000,   100000, 100000,  20,         20,        20,        20,        20,        20,        20,       20 ]

    # instantiate evolutionary computation algorithm with random seed
    my_ec = ec.emo.PAES(rand)
    # tournament sampling of individuals from population (<num_selected> individuals
    # are chosen based on best fitness performance in tournament)
    my_ec.selector = ec.selectors.tournament_selection
    #toggle variators
    # biased coin flip to determine whether 'mom' or 'dad' element is passed to offspring design
    # gaussian mutation which makes use of bounder function as specified in
    # --> my_ec.evolve(...,bounder=ec.BOunder(minParamValues, maxParamValues),...)
    my_ec.variator = [ec.variators.uniform_crossover,   
                     ec.variators.gaussian_mutation]    
    # existing generation is replaced by offspring, with elitism (<num_elites> existing
    # individuals will survive if they have better fitness than offspring)
    my_ec.replacer = ec.replacers.generational_replacement    
    # termination dictated by number of evaluations that have been run
    my_ec.terminator = ec.terminators.evaluation_termination  
    #toggle observers
    my_ec.observer = [ec.observers.stats_observer,  # print evolutionary computation statistics
                      ec.observers.best_observer]   # print the best individual in the population to screen
    #call evolution iterator
    final_pop = my_ec.evolve(generator=generate_netparams,  # assign design parameter generator to iterator parameter generator
                          evaluator=my_evaluator,           # assign fitness function to iterator evaluator
                          pop_size=10,                      # each generation of parameter sets will consist of 10 individuals
                          maximize=False,                   # best fitness corresponds to minimum value
                          bounder=ec.Bounder(minParamValues, maxParamValues), # boundaries for parameter set ([probability, weight, delay])
                          max_evaluations=50,               # evolutionary algorithm termination at 50 evaluations
                          num_selected=10,                  # number of generated parameter sets to be selected for next generation
                          mutation_rate=0.8,                # rate of mutation
                          num_inputs=11,                    # len([probability, weight, delay])
                          num_elites=1)                     # 1 existing individual will survive to next generation if it has better fitness than an individual selected by the tournament selection

    # plot raster of top solutions
    #final_pop.sort(reverse=True)                            # sort final population so best fitness (minimum difference) is first in list
    #bestCand = final_pop[0].candidate                       # bestCand <-- individual @ start of list

    final_arc = my_ec.archive
    RCMembrane, t = simulator(final_arc[0]._candidate[0], final_arc[0]._candidate[1], final_arc[0]._candidate[2],
                              final_arc[0]._candidate[3], final_arc[0]._candidate[4], final_arc[0]._candidate[5],
                              final_arc[0]._candidate[6], final_arc[0]._candidate[7], final_arc[0]._candidate[8],
                              final_arc[0]._candidate[9], final_arc[0]._candidate[10])
    plt.figure()
    plt.plot(t, RCMembrane)
    plt.title('RC Membrane Potential')
    plt.show()
