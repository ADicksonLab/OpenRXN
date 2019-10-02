"""
Derivative function builders take a list of sources and sinks
upon initialization.  Their self.deriv_func functions are used
by the integrator to move the system forward in time.

Rate constants given to DerivFuncBuilder are assumed to be 
in units of 1/s.  Units are stripped upon initialization.
"""

from openrxn import unit
import numpy as np

class DerivFuncBuilder(object):
    """
    Derivative functions are built with sources and sinks, 
    which have the form:
    
    (k, [quantities to multiply], n_per]

    k is a rate constant
    [quantities to multiply] is a list of Q indexes
    and n_per is the number of this species that are 
    produced (for sources) or consumed (for sinks) per reaction
    """
    def __init__(self,sources,sinks,sources_reservoir):
        #
        # sources and sinks are kept in the form (prefactor, inds)
        # where inds is a list of indexes of Q to multiply and prefactor is a scalar
        #
        self.terms = []
        for s in sources:
            pref = s[0].to(1/unit.sec).magnitude*s[2]
            if pref > 0:
                self.terms.append((pref,s[1]))
                
        for s in sinks:
            pref = s[0].to(1/unit.sec).magnitude*s[2]
            if pref > 0:
                self.terms.append((-pref,s[1]))

        self.sources_reservoir = [(s[0].to(unit.nm**3/unit.sec).magnitude*s[2], s[1]) for s in sources_reservoir]
        
    def deriv_func(self,Q,t):
        dqdt = 0
        for tup in self.terms:
            dqdt += tup[0]*np.prod(Q[tup[1]])
        for tup in self.sources_reservoir:
            dqdt += tup[0] * tup[1](t)
        return dqdt
