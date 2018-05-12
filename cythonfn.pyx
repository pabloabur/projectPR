import numpy as np

#pythran export runge_kutta(float, float64[], float, float, float, int)
def runge_kutta(t, x, timeStep, timeStepByTwo, timeStepBySix, MUnumber, unit, iIonic, G, iInjected, EqCurrent_nA, capacitanceInv):
    k1 = dVdt(t, x, MUnumber, unit, iIonic, G, iInjected, EqCurrent_nA, capacitanceInv)
    k2 = dVdt(t + timeStepByTwo, x + timeStepByTwo * k1, MUnumber, unit, iIonic, G, iInjected, EqCurrent_nA, capacitanceInv)
    k3 = dVdt(t + timeStepByTwo, x + timeStepByTwo * k2, MUnumber, unit, iIonic, G, iInjected, EqCurrent_nA, capacitanceInv)
    k4 = dVdt(t + timeStep, x + timeStep * k3, MUnumber, unit, iIonic, G, iInjected, EqCurrent_nA, capacitanceInv)
    
    return x + timeStepBySix * (np.add(
                                       np.add(
                                             np.add(k1, k2, order = 'C'), np.add(k2, k3, order='C')
                                             ),
                                       np.add(k3, k4, order='C'), order='C')
                               )

def dVdt(t, V, MUnumber, unit, iIonic, G, iInjected, EqCurrent_nA, capacitanceInv):
    #k = 0
    for i in xrange(MUnumber):
        for j in xrange(unit[i].compNumber):
            iIonic.itemset(i*unit[0].compNumber+j,
                                unit[i].compartment[j].computeCurrent(t,
                                                                        V.item(i*unit[0].compNumber+j)))
            #k += 1
          
    return (iIonic + G.dot(V) + iInjected + EqCurrent_nA) * capacitanceInv
