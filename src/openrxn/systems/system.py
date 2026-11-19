"""Systems are fully "realized" models, with a concentration
values specified for each compartment at a given time point
in system.state .
A system is moved forward in time with the system.propagate()
function.

In detail, all systems use a list of functions (system.dqdt), 
where system.dqdt[i] returns the rate of change of quantity i.

The system need not store the values of every species in every 
compartment.  Reporter functions can be attached to system objects
which can return the sum or average the concentrations in different
sets of compartments.  Reporters can also be configured to return 
ALL of the data and this is enabled by default.

After running, results (a.k.a. the reports from the reporters) are 
stored in system.results."""

from openrxn import unit
from openrxn.systems.state import State

import numpy as np
import logging

EPSILON = 1e-8

class System(object):

    def __init__(self, flatmodel, init_state=None, reporters=[]):
        """Systems must be initialized with FlatModel objects.
        initial states can be specified in the init_state argument,
        but care must be taken to ensure that this is compatible
        with the Model.

        It is strongly recommended to let the System initialize its
        own State, and then set its initial values using a set of
        selections and assignments.

        e.g. 
        s = System(my_flatmodel)
        species_a_bottom_layer = np.where(np.logical_and(
                                    s.state.z_pos < 1, s.state.species == a.ID))
        s.state.set_q(species_a_bottom_layer, 1 * ureg.mol)
        """

        self.model = flatmodel
        
        if init_state != None:
            self.state = init_state
        else:
            self.state = State(model=self.model)

        self.reporters = []
        self.reporters += reporters

    def add_reporter(self,reporter):
        self.reporters.append(reporter)

    def add_reporters(self,reporters):
        self.reporters += reporters

    def run(self,total_time,**kwargs):
        """
        Runs the system forward in time using the system-specific
        self.propagate function,

        Takes the total time of integration as an argument.  Other
        keyword arguments (kwargs) are passed to the self.propagate
        function.

        Note:  this function assumes that the initial state vector (q_val)
        has already been set, and that reporters have already been 
        defined and attached to the system.

        Every reporter records the initial state at t = 0 and then every
        r.freq.  With reporters the run is split into segments at the
        report times; the returned result covers the whole run from 0 to
        total_time (see self._combine_results).
        """

        report_freqs = [r.freq for r in self.reporters]

        # add endpoints as default, even if reporters aren't present
        checkpoints = [0,total_time]
        for freq in report_freqs:
            n = int(total_time/freq + EPSILON) + 1
            checkpoints += [freq*i for i in range(n)]

        checkpoints = list(set(checkpoints))
        checkpoints.sort()

        if report_freqs:
            # merge checkpoints that differ only by round-off, keeping 0 and total_time exact
            tol = EPSILON*min(report_freqs)
            merged = [0]
            for c in checkpoints:
                if tol < c - merged[-1] and c < total_time - tol:
                    merged.append(c)
            checkpoints = merged + [total_time]

        for r in self.reporters:
            r.report(checkpoints[0], self.state.q_val.copy())

        segments = []
        for i in range(len(checkpoints)-1):
            init_t = checkpoints[i]
            final_t = checkpoints[i+1]
            
            result = self.propagate((init_t,final_t),**kwargs)
            segments.append(result)
            if 'final_t' in result:
                checkpoints[i+1] = result['final_t']

            logging.info("Reached checkpoint: t = {0}".format(checkpoints[i+1]))
            
            for r in self.reporters:
                # check whether final_t is a multiple of r.freq (up to round-off)
                if abs(final_t/r.freq - round(final_t/r.freq)) < EPSILON:
                    r.report(checkpoints[i+1], self.state.q_val.copy())

        return self._combine_results(segments)

    def _combine_results(self, segments):
        """Merge the per-segment results of run() into one result.
        The base class returns the last segment (it holds the final state);
        ODESystem joins the trajectories."""
        return segments[-1]
        
    def propagate(self,t_interval,**kwargs):
        raise NotImplementedError
        
        
