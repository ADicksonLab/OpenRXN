"""ODE systems have continuous values for species quantities
and are propagated forward in time using ODE solvers like
scipy's solve_idp function."""

from openrxn import unit
from openrxn.systems.state import State
from openrxn.systems.deriv import DerivFuncBuilder
from openrxn.systems.system import System
from openrxn.compartments.compartment import Reservoir
from openrxn.connections import DivByVConnection

from scipy.integrate import solve_ivp
import numpy as np
import logging

EPSILON = 1e-8

class ODESystem(System):

    def __init__(self, *args, **kwargs):

        super().__init__(*args,**kwargs)
        self.NA = 6.022e23
        
        self.dqdt = self._build_dqdt()

    def set_q(self,idxs,Q):
        """Set the state.q_val array at the specified indexes
        to the value Q.

        idxs : list, int
        List of indexes to set.

        Q  : Quantity 
        If unitless, assumed to be number counts of species
        If mol/L is passed, it uses compartment values to
        convert to mol.
        """

        mult_volume = False
        mult_na = False
        if hasattr(Q,'units'):
            if Q.units == unit.mol/unit.L:
                mult_volume = True
                mult_na = True
            elif Q.units == unit.mol:
                mult_na = True
            elif Q.units != unit.dimensionless:
                raise ValueError("Quantity values should be either mol, mol/L or dimensionless")
        else:
            Q *= unit.dimensionless
 
        for i in idxs:
            if mult_volume:
                tmp = self.model.compartments[self.state.compartment[i]].volume * Q
            else:
                tmp = Q
            if mult_na:
                tmp *= self.NA

            # removing units for q_val array (cast in units of numbers of molecules)
            self.state.q_val[i] = tmp.magnitude
        
    def _build_dqdt(self):
        """Uses a model to build a list of derivative functions, with 
        indices that are consistent with the state vector.  self.dqdt[i], 
        calculates the rate of change of self.q_val[i], which holds the 
        amount of quantity i (a particular species in a particular compartment)
        in mol.

        Formulates each derivative equation as:
        
        dq[i]/dt =  sum_{j in sources}  k_j * prod_{k} q_k * n_per_ij
                  - sum_{j in sinks}    k_j * prod_{k} q_k * n_per_ij

        To construct the derivative function, it passes a list of 
        sources, where the elements of the list are:

        (k_j, [q_k0, q_k1...], n_per_ij)

        and a similar list of sinks.  The quantity n_per_ij is the number of 
        species i that are produced (or consumed) by reaction j.
        The rates passed are in units such that k_j * prod_{k} q_k has units 
        of s^-1.

        There is a special list of sources that come from Reservoir compartments,
        whose concentrations are not part of the state vector, but instead 
        vary deterministically as a function of time.  The elements of 
        source_reservoir lists are formatted as:

        (k_j, conc_func, n_per_ij)

        where conc_func is a link to the concentration function, which returns
        the concentration of the reservoir, given time as an input.

        Inputs:
        
        state : openrxn.systems.state.State object
        model : openrxn.model.FlatModel object
        """
        dqdt = []
        for i in range(self.state.size):
            # collect source and sink terms for this species in this compartment

            c = self.model.compartments[self.state.compartment[i]]
            s = self.state.species[i]

            sources = []
            sources_reservoir = []
            sinks = []
            
            # look through the reactions in this compartment for ones that
            # involve this species
            for r in c.reactions:
                if s in r.reactant_IDs:
                    s_idx = r.reactant_IDs.index(s)
                    if r.kf > 0:
                        # append forward reaction to sinks
                        q_list = []
                        n_r = 0
                        for j,x in enumerate(r.reactants):
                            q_list += [self.state.index[c.ID][x.ID]]*r.stoich_r[j]
                            n_r += r.stoich_r[j]
                        if n_r - 1 > 0 and c.volume is not None:
                            vol_fac = (self.NA*c.volume/unit.mol)**(n_r-1)
                            rate = r.kf/vol_fac
                        else:
                            rate = r.kf
                        sinks.append((rate, q_list, r.stoich_r[s_idx]))
                    if r.kr > 0:
                        # append reverse reaction to sources
                        q_list = []
                        n_p = 0
                        for j,x in enumerate(r.products):
                            q_list += [self.state.index[c.ID][x.ID]]*r.stoich_p[j]
                            n_p += r.stoich_p[j]
                        if n_p - 1 > 0 and c.volume is not None:
                            vol_fac = (self.NA*c.volume/unit.mol)**(n_p-1)
                            rate = r.kr/vol_fac
                        else:
                            rate = r.kr
                        sources.append((rate, q_list, r.stoich_r[s_idx]))
                    
                if s in r.product_IDs:
                    s_idx = r.product_IDs.index(s)
                    if r.kf > 0:
                        # append forward reaction to sources
                        q_list = []
                        n_r = 0
                        for j,x in enumerate(r.reactants):
                            q_list += [self.state.index[c.ID][x.ID]]*r.stoich_r[j]
                            n_r += r.stoich_r[j]
                        if n_r - 1 > 0 and c.volume is not None:
                            vol_fac = (self.NA*c.volume/unit.mol)**(n_r-1)
                            rate = r.kf/vol_fac
                        else:
                            rate = r.kf
                        sources.append((rate, q_list, r.stoich_p[s_idx]))

                    if r.kr > 0:
                        # append reverse reaction to sinks
                        q_list = []
                        n_p = 0
                        for j,x in enumerate(r.products):
                            q_list += [self.state.index[c.ID][x.ID]]*r.stoich_p[j]
                            n_p += r.stoich_p[j]
                        if n_p - 1 > 0 and c.volume is not None:
                            vol_fac = (self.NA*c.volume/unit.mol)**(n_p-1)
                            rate = r.kr/vol_fac
                        else:
                            rate = r.kr
                        sinks.append((rate, q_list, r.stoich_p[s_idx]))

            # look through connections for those that involve this species
            for other_lab, conn in c.connections.items():
                if s in conn[1].species_rates:
                                  
                    # add "out" diffusion process
                    if isinstance(conn,DivByVConnection):
                        sinks.append((conn[1].species_rates[s][0]/c.volume, [i], 1))
                    else:
                        sinks.append((conn[1].species_rates[s][0], [i], 1))

                    # add "in" diffusion process
                    if isinstance(self.model.compartments[other_lab],Reservoir):
                        sources_reservoir.append((conn[1].species_rates[s][1],
                                        self.model.compartments[other_lab].conc_funcs[s],
                                        1))
                    else:
                        if isinstance(conn,DivByVConnection):
                            sources.append((conn[1].species_rates[s][1]/conn[0].volume,
                                            [self.state.index[other_lab][s]],
                                            1))
                        else:
                            sources.append((conn[1].species_rates[s][1],
                                            [self.state.index[other_lab][s]],
                                            1))
                    
            dqdt.append(DerivFuncBuilder(sources, sinks, sources_reservoir))

        return dqdt
                            
    def _dQ_dt(self,t,Q):
        return np.array([builder.deriv_func(Q,t) for builder in self.dqdt])
    
    def propagate(self,t_interval,**kwargs):
        """For ODE systems, propagate directly calls the scipy
        solve_ivp function.  state.q_val is also updated."""

        result = solve_ivp(self._dQ_dt,t_interval,self.state.q_val,**kwargs)
        self.state.q_val = result.y[:,-1]
        
        return result
