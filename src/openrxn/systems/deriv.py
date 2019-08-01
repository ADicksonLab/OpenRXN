"""
Derivative function builders take a list of sources and sinks
upon initialization.  Their self.deriv_func functions are used
by the integrator to move the system forward in time.

Rate constants given to DerivFuncBuilder are assumed to be 
in units of 1/s.  Units are stripped upon initialization.
"""

from openrxn import unit

class DerivFuncBuilder(object):

    def __init__(self,sources,sinks):
        self.sources = [(s[0].to(1/unit.sec).magnitude,s[1]) for s in sources]
        self.sinks = [(s[0].to(1/unit.sec).magnitude,s[1]) for s in sinks]
        
    def deriv_func(self,Q,t):
        dqdt = 0
        for tup in self.sources:
            tmp = 1
            for idx in tup[1]:
                tmp *= Q[idx]
            dqdt += tup[0]*tmp
        for tup in self.sinks:
            tmp = 1
            for idx in tup[1]:
                tmp *= Q[idx]
            dqdt -= tup[0]*tmp                
        return dqdt

