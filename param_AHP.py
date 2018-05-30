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

def evaluate_netparams(candidates, args):
    fitnessCandidates = []

    for icand,cand in enumerate(candidates):
        # Make simulations
        RCMembrane, t = simulator(cand[0], cand[1])

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

        # calculate fitness for this candidate
        fitness = abs(targetAHPa - AHPa)  # minimize absolute difference in firing rate

        # add to list of fitness for each candidate
        fitnessCandidates.append(fitness)

        print "\n Child Candidate %d: gKs = %d, gKf = %d\n"%(icand, cand[0], cand[1])
        print "AHP amplitude: %.2f; AMP duration: %.2f; Fitness: %.2f;\n"%(AHPa, AHPt,fitness)

    return fitnessCandidates

def constrained_tournament_selection(random, population, args):
    num_selected = args.setdefault('num_selected', 1)
    constraint_func = args.setdefault('constraint_function', None)
    tournament_size = 2
    pop = list(population)
    selected = []
    for _ in range(num_selected):
        tournament = random.sample(pop, tournament_size)
        # If there is not a constraint function,
        # just do regular tournament selection.
        if constraint_func is None:
            selected.append(max(tournament))
        else:
            cons = [constraint_func(t.candidate) for t in tournament]
            # If no constraints are violated, just do 
            # regular tournament selection.
            if max(cons) == 0:
                selected.append(max(tournament))
            # Otherwise, choose the least violator 
            # (which may be a non-violator).
            else:
                selected.append(tournament[cons.index(min(cons))])
    return selected

def simulator(gKf, gKs):
    conf = Configuration('confuchiyama.rmto')

    # Number of cells
    idx = np.where(conf.confArray['f0']=='Number_RC_ext')[0][0]
    conf.confArray['f1'][idx] = 1
    
    # Duration of simulation
    conf.simDuration_ms = 100
    
    # Parameters from java
    ## Potassium conductances
    idx = np.where(conf.confArray['f0']=='gmax_Ks:RC_ext-@soma')[0][0]
    conf.confArray['f1'][idx] = gKs
    idx = np.where(conf.confArray['f0']=='gmax_Kf:RC_ext-@soma')[0][0]
    conf.confArray['f1'][idx] = gKf

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
    rand.seed(1)

    # target mean firing rate in Hz
    targetAHPt = 30
    targetAHPa = -2

    # min and max allowed value for each AHP parameter optimized:
    #                 gKs,      gKf,
    minParamValues = [500,      500]
    maxParamValues = [10000,   10000]

    # instantiate evolutionary computation algorithm with random seed
    my_ec = ec.EvolutionaryComputation(rand)
    # tournament sampling of individuals from population (<num_selected> individuals
    # are chosen based on best fitness performance in tournament)
    my_ec.selector = constrained_tournament_selection  
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
                          evaluator=evaluate_netparams,     # assign fitness function to iterator evaluator
                          pop_size=10,                      # each generation of parameter sets will consist of 10 individuals
                          maximize=False,                   # best fitness corresponds to minimum value
                          bounder=ec.Bounder(minParamValues, maxParamValues), # boundaries for parameter set ([probability, weight, delay])
                          max_evaluations=50,               # evolutionary algorithm termination at 50 evaluations
                          num_selected=10,                  # number of generated parameter sets to be selected for next generation
                          mutation_rate=0.1,                # rate of mutation
                          num_inputs=2,                     # len([probability, weight, delay])
                          num_elites=1)                     # 1 existing individual will survive to next generation if it has better fitness than an individual selected by the tournament selection

    # plot raster of top solutions
    final_pop.sort(reverse=True)                            # sort final population so best fitness (minimum difference) is first in list
    bestCand = final_pop[0].candidate                       # bestCand <-- individual @ start of list

    RCMembrane, t = simulator(bestCand[0], bestCand[1])
    plt.figure()
    plt.plot(t, RCMembrane)
    plt.title('RC Membrane Potential')
    plt.show()
