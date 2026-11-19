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
from scipy.integrate._ivp.ivp import OdeResult
import numpy as np
import logging

EPSILON = 1e-8

class ODESystem(System):

    def __init__(self, *args, **kwargs):

        super().__init__(*args,**kwargs)
        
        self.NA = 6.022e23
        self.dqdt = self._build_dqdt()
        # C1: harvest the per-state terms above into vectorized, order-grouped
        # arrays for a fast right-hand side (self._dQ_dt).
        self._build_rhs()
        self._build_vicinities()

    def set_q(self,idxs,Q):
        """Set the state.q_val array at the specified indexes
        to the value Q.

        idxs : list, int
        List of indexes to set.

        Q  : number or Quantity
        - unitless or dimensionless: number of molecules
        - an amount (mol, mmol, ...): converted to molecules (x NA)
        - a concentration (mol/L, unit.molar, uM, ...): multiplied by
          each compartment's volume, then converted to molecules.  The
          compartment must have a 3D volume (e.g. L or nm**3).
        """

        if np.ndim(idxs) == 0:
            idxs = [idxs]

        if not hasattr(Q,'units') or Q.dimensionless:
            n = Q.to(unit.dimensionless).magnitude if hasattr(Q,'units') else Q
            for i in idxs:
                self.state.q_val[i] = n
        elif Q.check('[substance]'):
            n = Q.to(unit.mol).magnitude*self.NA
            for i in idxs:
                self.state.q_val[i] = n
        elif Q.check('[substance]/[volume]'):
            for i in idxs:
                c_ID = self.state.compartment[i]
                vol = self.model.compartments[c_ID].volume
                if vol is None or not hasattr(vol,'units') or not vol.check('[volume]'):
                    raise ValueError("Compartment {0} needs a 3D volume (e.g. L or nm**3) to set a "
                                     "concentration; it has volume {1}".format(c_ID, vol))
                self.state.q_val[i] = (vol*Q).to(unit.mol).magnitude*self.NA
        else:
            raise ValueError("Quantity values should be an amount (e.g. mol), a concentration "
                             "(e.g. mol/L, molar) or dimensionless; got {0}".format(Q.units))
        
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
        state = self.state

        dqdt = []
        for i in range(state.size):
            # collect source and sink terms for this species in this compartment

            c = self.model.compartments[self.state.compartment[i]]
            s = state.species[i]
            vol = (self.NA*c.volume/unit.mol) if c.volume is not None else None
            state_idx = state.index[c.ID]

            sources = []
            sources_reservoir = []
            sinks = []
            
            # look through the reactions in this compartment for ones that
            # involve this species

            for r in c.reactions:
                if r.vicinity is not None:
                    # handled by _build_vicinities
                    continue
                if s in r.reactant_IDs:
                    s_idx = r.reactant_IDs.index(s)
                    if r.kf > 0:
                        # append forward reaction to sinks
                        q_list = []
                        n_r = 0
                        for j,x in enumerate(r.reactants):
                            q_list.extend([state_idx[x.ID]]*r.stoich_r[j])
                            n_r += r.stoich_r[j]
                        if n_r - 1 > 0 and vol is not None:
                            vol_fac = vol**(n_r-1)
                            rate = r.kf/vol_fac
                        else:
                            rate = r.kf
                        sinks.append((rate, q_list, r.stoich_r[s_idx]))
                    if r.kr > 0:
                        # append reverse reaction to sources
                        q_list = []
                        n_p = 0
                        for j,x in enumerate(r.products):
                            q_list.extend([state_idx[x.ID]]*r.stoich_p[j])
                            n_p += r.stoich_p[j]
                        if n_p - 1 > 0 and vol is not None:
                            vol_fac = vol**(n_p-1)
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
                            q_list.extend([state_idx[x.ID]]*r.stoich_r[j])
                            n_r += r.stoich_r[j]
                        if n_r - 1 > 0 and vol is not None:
                            vol_fac = vol**(n_r-1)
                            rate = r.kf/vol_fac
                        else:
                            rate = r.kf
                        sources.append((rate, q_list, r.stoich_p[s_idx]))

                    if r.kr > 0:
                        # append reverse reaction to sinks
                        q_list = []
                        n_p = 0
                        for j,x in enumerate(r.products):
                            q_list.extend([state_idx[x.ID]]*r.stoich_p[j])
                            n_p += r.stoich_p[j]
                        if n_p - 1 > 0 and vol is not None:
                            vol_fac = vol**(n_p-1)
                            rate = r.kr/vol_fac
                        else:
                            rate = r.kr
                        sinks.append((rate, q_list, r.stoich_p[s_idx]))

            # look through connections for those that involve this species
            for other_lab, conn in c.connections.items():
                # conn is a tuple: (neighbor_compartment, connection_object)
                neighbor, connection = conn
                if s in connection.species_rates:

                    # add "out" diffusion process
                    if isinstance(connection,DivByVConnection):
                        sinks.append((connection.species_rates[s][0]/c.volume, [i], 1))
                    else:
                        sinks.append((connection.species_rates[s][0], [i], 1))

                    # add "in" diffusion process
                    if isinstance(neighbor,Reservoir):
                        sources_reservoir.append((connection.species_rates[s][1],
                                        neighbor.conc_funcs[s],
                                        1))
                    else:
                        if isinstance(connection,DivByVConnection):
                            sources.append((connection.species_rates[s][1]/neighbor.volume,
                                            [state.index[other_lab][s]],
                                            1))
                        else:
                            sources.append((connection.species_rates[s][1],
                                            [state.index[other_lab][s]],
                                            1))
                    
            dqdt.append(DerivFuncBuilder(sources, sinks, sources_reservoir))

        return dqdt
                            
    def _build_rhs(self):
        """Harvest the per-state terms from the DerivFuncBuilder list into
        flat, order-grouped arrays for a vectorized right-hand side.

        Each term is (prefactor, indices) contributing
        prefactor * prod(Q[indices]) to one state variable. 
        """
        groups = {}            # order -> {'targets': [], 'prefs': [], 'idx': []}
        reservoir_terms = []   # (target_state, coef, conc_func)
        for i, builder in enumerate(self.dqdt):
            for pref, inds in builder.terms:
                order = len(inds)
                g = groups.setdefault(order, {'targets': [], 'prefs': [], 'idx': []})
                g['targets'].append(i)
                g['prefs'].append(pref)
                g['idx'].append(list(inds))
            for coef, conc_func in builder.sources_reservoir:
                reservoir_terms.append((i, coef, conc_func))

        self._rhs_groups = []
        for order, g in groups.items():
            targets = np.asarray(g['targets'], dtype=np.intp)
            prefs = np.asarray(g['prefs'], dtype=float)
            idx = np.asarray(g['idx'], dtype=np.intp) if order > 0 else None
            self._rhs_groups.append((order, targets, prefs, idx))
        self._reservoir_terms = reservoir_terms

    def _build_vicinities(self):
        """Collect the reactions tagged with a Vicinity, one group per vicinity
        and compartment, for the quasi-steady-state rebinding rate law:

        L_v = (k_drift * q_L + sum_j f_j) / (k_drift + sum_i a_i)

        f_j are the fluxes of reaction directions that release the ligand into
        the vicinity, a_i * L_v the fluxes of directions that consume it
        (a_i = prefactor * product of the other reactants), and the bulk ligand
        exchanges with the vicinity at k_drift * (L_v - q_L).

        Each element of self._vicinities is (k_drift, lig, captures, releases),
        where captures and releases are lists of (prefactor, inds, deltas).
        """
        self._vicinities = []
        for c in self.model.compartments.values():
            state_idx = self.state.index[c.ID]
            vol = (self.NA*c.volume/unit.mol) if c.volume is not None else None
            groups = {}
            for r in c.reactions:
                if r.vicinity is None:
                    continue
                lig = state_idx[r.vicinity.ligand.ID]
                g = groups.setdefault(id(r.vicinity), (r.vicinity.k_drift.to(1/unit.sec).magnitude, lig, [], []))
                for k, side_in, st_in, side_out, st_out in [(r.kf, r.reactants, r.stoich_r, r.products, r.stoich_p),
                                                            (r.kr, r.products, r.stoich_p, r.reactants, r.stoich_r)]:
                    if not k > 0:
                        continue
                    inds = [state_idx[x.ID] for x, n in zip(side_in, st_in) for _ in range(n)]
                    rate = k/vol**(len(inds)-1) if len(inds) > 1 and vol is not None else k
                    deltas = {}
                    for x, n in zip(side_in, st_in):
                        deltas[state_idx[x.ID]] = deltas.get(state_idx[x.ID], 0) - n
                    for x, n in zip(side_out, st_out):
                        deltas[state_idx[x.ID]] = deltas.get(state_idx[x.ID], 0) + n
                    # the ligand is exchanged with the vicinity, not the bulk
                    capture = lig in inds
                    if capture:
                        inds.remove(lig)
                    deltas[lig] += 1 if capture else -1
                    term = (rate.to(1/unit.sec).magnitude, np.asarray(inds, dtype=np.intp),
                            [(i, d) for i, d in deltas.items() if d != 0])
                    (g[2] if capture else g[3]).append(term)
            self._vicinities += groups.values()

    def _dQ_dt(self, t, Q):
        """Vectorized right-hand side: computes dQ/dt for the whole state at
        once via order-grouped products and scatter-adds."""
        n = self.state.size
        dQ = np.zeros(n)
        for order, targets, prefs, idx in self._rhs_groups:
            if order == 0:
                vals = prefs
            else:
                vals = prefs * np.prod(Q[idx], axis=1)
            dQ += np.bincount(targets, weights=vals, minlength=n)
        for i, coef, conc_func in self._reservoir_terms:
            dQ[i] += coef * conc_func(t)
        for kd, lig, caps, rels in self._vicinities:
            a = [pref * np.prod(Q[inds]) for pref, inds, _ in caps]
            f = [pref * np.prod(Q[inds]) for pref, inds, _ in rels]
            den = kd + sum(a)
            Lv = (kd * Q[lig] + sum(f)) / den
            # = kd * (Lv - Q[lig])
            dQ[lig] += kd * (sum(f) - Q[lig] * sum(a)) / den
            for flux, (_, _, deltas) in zip([x * Lv for x in a] + f, caps + rels):
                for i, d in deltas:
                    dQ[i] += d * flux
        return dQ

    def _dQ_dt_original(self, t, Q):
        """Original per-state-variable RHS. Not used."""
        return np.array([builder.deriv_func(Q, t) for builder in self.dqdt])

    def _jac(self, t, Q):
        """Analytic Jacobian d(dQ/dt)/dQ, assembled from the same order-grouped
        terms as the RHS.

        For pref * prod(Q[inds]), calculate the derivative with respect to each species by multiplying pref 
        by all factors except the one being differentiated. If a species appears more than once, add the contributions from each occurrence.
        Order-0 (constant) and reservoir terms have zero derivative.
        """
        n = self.state.size
        J = np.zeros((n, n))
        for order, targets, prefs, idx in self._rhs_groups:
            if order == 0:
                continue
            for p in range(order):
                if order == 1:
                    prod_others = np.ones(len(prefs))
                else:
                    others = [q for q in range(order) if q != p]
                    prod_others = np.prod(Q[idx[:, others]], axis=1)
                np.add.at(J, (targets, idx[:, p]), prefs * prod_others)
        for kd, lig, caps, rels in self._vicinities:
            a = [self._mono(pref, inds, Q) for pref, inds, _ in caps]
            f = [self._mono(pref, inds, Q) for pref, inds, _ in rels]
            A, F = sum(x for x, _ in a), sum(x for x, _ in f)
            dA, dF = sum((g for _, g in a), np.zeros(n)), sum((g for _, g in f), np.zeros(n))
            den = kd + A
            Lv = (kd * Q[lig] + F) / den
            # quotient rule, including the derivative of the denominator
            dLv = (dF - Lv * dA) / den
            dLv[lig] += kd / den
            J[lig] += kd * (dF - (Q[lig] + (F - Q[lig] * A) / den) * dA) / den
            J[lig, lig] -= kd * A / den
            for (x, g), (_, _, deltas) in zip(a, caps):
                for i, d in deltas:
                    J[i] += d * (g * Lv + x * dLv)
            for (x, g), (_, _, deltas) in zip(f, rels):
                for i, d in deltas:
                    J[i] += d * g
        return J

    @staticmethod
    def _mono(pref, inds, Q):
        """ Value and gradient of the mass action term prefactor * prod(Q[inds]).
            Vicinity rates depend on L_v, so their gradients must be combined before updating the Jacobian.
        
        """
        g = np.zeros(Q.size)
        for p in range(len(inds)):
            g[inds[p]] += pref * np.prod(np.delete(Q[inds], p))
        return pref * np.prod(Q[inds]), g

    def propagate(self, t_interval, method='Radau', rtol=1e-6, atol=1e-6, jac=None, **kwargs):
        """For ODE systems, propagate calls scipy's solve_ivp and updates
        state.q_val.

        If t_eval is given, only the points inside t_interval are used, so
        run() can pass one t_eval to every segment.
        """
        if jac is None and method in ('Radau', 'BDF', 'LSODA'):
            jac = self._jac
        if jac is not None:
            kwargs['jac'] = jac          
        t_eval = kwargs.get('t_eval')
        if t_eval is not None:
            t_eval = np.asarray(t_eval, dtype=float)
            lo, hi = min(t_interval), max(t_interval)
            if np.any((t_eval < lo) | (t_eval > hi)):
                kwargs['t_eval'] = t_eval[(t_eval >= lo) & (t_eval <= hi)]
        result = solve_ivp(self._dQ_dt, t_interval, self.state.q_val,
                           method=method, rtol=rtol, atol=atol, **kwargs)
        # a copy, so a later set_q cannot overwrite the last point of this result
        self.state.q_val = result.y[:, -1].copy()
        return result

    def _combine_results(self, segments):
        """Join the solve_ivp results of run() segments into one OdeResult
        whose t and y cover the whole run.  A time point shared by two
        segments is kept once.  Dense output (sol) and events are not
        joined: they are set to None when there is more than one segment."""
        if len(segments) == 1:
            return segments[0]
        ts, ys = [], []
        last_t = None
        for seg in segments:
            t, y = np.asarray(seg.t), np.asarray(seg.y)
            if last_t is not None and t.size and t[0] == last_t:
                t, y = t[1:], y[:, 1:]
            if t.size:
                ts.append(t); ys.append(y); last_t = t[-1]
        out = OdeResult(t=np.concatenate(ts), y=np.concatenate(ys, axis=1),
                        sol=None, t_events=None, y_events=None,
                        nfev=sum(s.nfev for s in segments),
                        njev=sum(s.njev for s in segments),
                        nlu=sum(s.nlu for s in segments),
                        status=segments[-1].status, message=segments[-1].message,
                        success=all(s.success for s in segments))
        return out
