"""Gillespie systems have discrete values for species quantities
and are propagated forward in time using stochastic algorithms
like the Gillespie algorithm.

self.processes is a list of reactions formatted as follows:

(k,q_list,delta_list)

where q_list is a list of quantities that should be multiplied by
k to determine the reaction rate, and 

delta_list is a list of tuples formatted as: (species_idx, delta)
that govern the number count updates if that reaction is chosen. 
"""

from openrxn import unit
from openrxn.systems.state import State
from openrxn.systems.deriv import DerivFuncBuilder
from openrxn.systems.system import System
from openrxn.propagators import Gillespie

import numpy as np
import logging

EPSILON = 1e-8

class GillespieSystem(System):

    def __init__(self, *args, **kwargs):

        super().__init__(*args,**kwargs)
        self.processes = self._build_processes()

    def propagate(self,t_interval,**kwargs):
        """
        Interfaces with openrxn.propagators.Gillespie()
        """

        new_q, final_t = Gillespie(self.processes,t_interval,self.state.q_val)
        self.state.q_val = new_q

        return (new_q, final_t)

    def _build_processes(self):
        """
        Processes is a list with elements of format:
        (rate,q_list,delta_list)

        where rate is a constant (in 1/s) that is multiplied by a series of 
        quantities (specified by q_list) that are in units of number counts, 
        in order to determine the reaction flux (number of reactions occurring 
        per second)

        The delta list describes how the quantities should be updated if this
        reaction is chosen to proceed, and has elements that are tuples of the 
        format (index, delta).  Delta for e.g. is usually +1 or -1.
        """
        processes = []
        for c in self.model.compartments.values():
            # first add reactions 
            for r in c.reactions:
                # add forward reaction (losing reactants, gaining products)
                
                # qlist: quantities that must be multiplied together to
                # get the reaction flux
                q_list = []
                delta_list = []
                vol_fac = 1
                n_r = 0
                for j,x in enumerate(r.reactants):
                    q_list += [self.state.index[c.ID][x.ID]]*r.stoich_r[j]
                    delta_list.append((self.state.index[c.ID][x.ID],-r.stoich_r[j]))
                    n_r += r.stoich_r[j]
                if hasattr(c,'volume') and n_r - 1 > 0:
                    vol_fac = c.volume**(n_r-1)

                for j,x in enumerate(r.products):
                    delta_list.append((self.state.index[c.ID][x.ID],r.stoich_p[j]))

                processes.append((self._cast_rate(r.kf/vol_fac),q_list,delta_list))

                # add reverse reaction (gaining reactants, losing products)
                q_list = []
                delta_list = []
                vol_fac = 1
                n_r = 0
                for j,x in enumerate(r.reactants):
                    delta_list.append((self.state.index[c.ID][x.ID],r.stoich_r[j]))

                for j,x in enumerate(r.products):
                    delta_list.append((self.state.index[c.ID][x.ID],-r.stoich_p[j]))
                    q_list += [self.state.index[c.ID][x.ID]]*r.stoich_p[j]
                    n_r += r.stoich_p[j]
                if hasattr(c,'volume') and n_r - 1 > 0:
                    vol_fac = c.volume**(n_r-1)

                processes.append((self._cast_rate(r.kr/vol_fac),q_list,delta_list))

            # then, add diffusion processes
            for other_lab, conn in c.connections.items():
                for s in conn[1].species_rates.keys():
                    # add "out" diffusion process
                    # Note: volumes must be defined if diffusion processes are occurring
                    #
                    processes.append(self._cast_rate(conn[1].species_rates[s][0]/c.volume),
                                     [self.state.index[c.ID][s]],
                                     [(self.state.index[c.ID][s],-1), (self.state.index[other_lab][s],1)])
                                     
        return processes

    def set_q(self,idxs,Q):
        """Set the state.q_val array at the specified indexes
        to the value Q.

        idxs : list, int
        List of indexes to set.

        Q  : Quantity 
        Must be a unitless integer.
        """

        if hasattr(Q,'units'):
            if Q.units != unit.dimensionless:
                raise ValueError("Quantity values for Gillespie systems must be dimensionless")
            else:
                Q = Q.magnitude

        self.state.q_val[idxs] = Q

    def _cast_rate(self,rate):
        # checks rates to make sure they are in the correct units (1/s), and then
        # strips the units away

        rate.ito(1/unit.sec)
        
        return rate.magnitude
