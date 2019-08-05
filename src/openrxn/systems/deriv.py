"""
Derivative function builders take a list of sources and sinks
upon initialization.  Their self.deriv_func functions are used
by the integrator to move the system forward in time.

Rate constants given to DerivFuncBuilder are assumed to be 
in units of 1/s.  Units are stripped upon initialization.
"""

from openrxn import unit

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
    def __init__(self,sources,sinks):
        self.sources = [(s[0].to(1/unit.sec).magnitude,s[1],s[2]) for s in sources]
        self.sinks = [(s[0].to(1/unit.sec).magnitude,s[1],s[2]) for s in sinks]
        
    def deriv_func(self,Q,t):
        dqdt = 0
        for tup in self.sources:
            tmp = 1
            for idx in tup[1]:
                tmp *= Q[idx]
            dqdt += tup[0]*tmp*tup[2]
        for tup in self.sinks:
            tmp = 1
            for idx in tup[1]:
                tmp *= Q[idx]
            dqdt -= tup[0]*tmp*tup[2]         
        return dqdt

